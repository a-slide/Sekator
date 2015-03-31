#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

try:
    # Standard library imports
    from multiprocessing import Value, Array, Process, Queue, cpu_count
    from time import time
    import gzip
    from os import path
    import ConfigParser
    import optparse
    import sys
    import os
    from gzip import open as gopen

    # Local Package import
    from AdapterTrimmer import AdapterTrimmer
    from QualityTrimmer import QualityTrimmer
    from Conf_file import write_example_conf
    from Fastq import FastqSeq, FastqReader
    from Sample import Sample
    from ProgressBar import ProgressBar

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sekator (object):
    """
    Sekator is a multiprocessing fastq file adapter and quality trimmer
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "Sekator 0.2"
    USAGE = "Usage: %prog -c Conf.txt [-i -h]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)
        optparser.add_option('-c', dest="conf_file",
            help= "Path to the configuration file [Mandatory]")
        optparser.add_option('-i', dest="init_conf", action='store_true',
            help= "Generate an example configuration file and exit [Facultative]")

        # Parse arguments
        options, args = optparser.parse_args()

        return Sekator(options.conf_file, options.init_conf)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, conf_file=None, init_conf=None):
        """
        Initialization function, parse options in command line and configuration file and verify
        their values. All self.variables are initialized explicitly in init
        """

        # Create a example conf file if needed
        if init_conf:
            print("Create an example configuration file in the current folder")
            write_example_conf()
            sys.exit(0)

        print("Initialize Sekator")
        # Parse the configuration file and verify the values of variables
        try:

            #verify if conf file was
            assert conf_file, "A path to the configuration file is mandatory"
            self.conf = conf_file

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=False)
            cp.read(self.conf)

            # General section
            self.min_size = cp.getint("general", "min_size")
            self.n_thread = cpu_count() if cp.getboolean("general", "auto_thread") else cp.getint("general", "n_thread")
            self.write_report = cp.getboolean("general", "write_report")
            self.compress_output = cp.getboolean("general", "compress_output")

            # Quality Trimming section
            self.left_trim = cp.getboolean("quality", "left_trim")
            self.right_trim = cp.getboolean("quality", "right_trim")
            self.quality_trim = self.left_trim or self.right_trim
            if self.quality_trim:
                self.win_size = cp.getint("quality", "win_size")
                self.step = cp.getint("quality", "step")
                self.qual_cutdown = cp.getint("quality", "qual_cutdown")

            # Adapter Trimming section
            self.adapter_trim = cp.getboolean("adapter", "adapter_trim")
            if self.adapter_trim:
                self.min_match_len = cp.getfloat("adapter", "min_match_len")
                self.min_match_score = cp.getfloat("adapter", "min_match_score")
                self.ssw_match = cp.getint("adapter", "ssw_match")
                self.ssw_mismatch = cp.getint("adapter", "ssw_mismatch")
                self.ssw_gapO = cp.getint("adapter", "ssw_gapO")
                self.ssw_gapE = cp.getint("adapter", "ssw_gapE")

            # Samples are a special case, since the number of sections is variable
            # Iterate only on sections starting by "sample", create Sample objects
            # And store them in a list
            self.sample_list = []
            for sample in [i for i in cp.sections() if i.startswith("sample")]:
                self.sample_list.append (Sample (
                    name = cp.get(sample, "name"),
                    R1_path = cp.get(sample, "R1_path"),
                    R2_path = cp.get(sample, "R2_path"),
                    adapter_list = [] if not self.adapter_trim else cp.get(sample, "adapter_list").split(),
                    compress_output = self.compress_output))

            # Values are tested in a private function
            self._test_values()

        # Handle the many possible errors occurring during conf file parsing or variable test
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            sys.exit(1)
        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            sys.exit(1)

        print ("All configuration file parameters are valid")

        self.buffer_size = 20 # Buffer size for writing in fastq file

    def __str__(self):
        msg = "SEKATOR CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def __call__ (self):
        """ Main function of the script """

        start_time = time()

        for n, sample in enumerate (self.sample_list):

            print ("ANALYSING SAMPLE {} ({}/{})".format(sample.name, n+1, len(self.sample_list)))

            print ("\tVerify Fastq and count the number of reads")
            n_read1 = self._count_fastq (sample.R1_path)
            n_read2 = self._count_fastq (sample.R2_path)
            assert n_read1 == n_read2, "Fastq R1 and Fastq R2 files do not contain the same number of reads"
            self.progress_bar = ProgressBar(total_seq = n_read1, number_step = 10)

            # Init generic shared memory counters
            self.total = Value('i', 0)
            self.pass_qual = Value('i', 0)
            self.pass_adapt = Value('i', 0)
            self.total_pass = Value('i', 0)

            # Define Quality Trimmer Object and specific shared memory counters
            if self.quality_trim:
                self.quality_trimmer = QualityTrimmer(
                    qual_cutdown = self.qual_cutdown,
                    win_size = self.win_size,
                    step = self.step,
                    min_size = self.min_size,
                    left_trim = self.left_trim,
                    right_trim = self.right_trim)

                self.qual_total = Value('i', 0)
                self.qual_untrimmed = Value('i', 0)
                self.qual_trimmed = Value('i', 0)
                self.qual_fail = Value('i', 0)
                self.qual_base_trimmed = Value('i', 0)
                self.qual_mean_sum = Value('i', 0)

            # Define Adapter Trimmer Object and specific shared memory counters
            if self.adapter_trim:
                self.adapter_trimmer = AdapterTrimmer(
                    adapter_list = sample.adapter_list,
                    min_size = self.min_size,
                    min_match_len = self.min_match_len,
                    min_match_score = self.min_match_score,
                    ssw_match = self.ssw_match,
                    ssw_mismatch = self.ssw_mismatch,
                    ssw_gapO = self.ssw_gapO,
                    ssw_gapE = self.ssw_gapE)

                self.adapt_total = Value('i', 0)
                self.adapt_untrimmed = Value('i', 0)
                self.adapt_trimmed = Value('i', 0)
                self.adapt_fail = Value('i', 0)
                self.adapt_base_trimmed = Value('i', 0)

            # Init queues for input file reading and output file writing (limited to 10000 objects)
            self.inq = Queue(maxsize=10000)
            self.outq = Queue(maxsize=10000)

            # Init processes for file reading, distributed filtering and file writing
            self.pin = Process(target=self.reader, args=(sample.R1_path, sample.R2_path))
            self.ps = [Process(target=self.filter, args=(i,)) for i in range(self.n_thread)]
            self.pout = Process(target=self.writer, args=(sample.R1_outname, sample.R2_outname))

            # Start processes
            print ("\tStarting fastq trimming")
            self.pin.start()
            self.pout.start()
            for p in self.ps:
                p.start()

            # Blocks until the process is finished
            self.pin.join()
            for i in range(len(self.ps)):
                self.ps[i].join()
            self.pout.join()
            print ("\tFastq trimming done")

            if self.write_report:
                self._write_report(sample.name, len(sample.adapter_list))

        print ("Done in {}s".format(round(time()-start_time, 3)))
        return(0)

    def reader(self, R1_path, R2_path):
        """
        Initialize HTSseq FastqReader to iterate over paired fastq files. Data are send in the
        in queue for the workers. Add n_thread STOP pills at the end of the inq for each worker.
        """

        # Init FastqReader generators
        R1_gen = FastqReader(R1_path)
        R2_gen = FastqReader(R2_path)
        n = 0

        # Iterate over reads in fastq files until exhaustion
        try:
            while True:
                read1 = R1_gen.next()
                read2 = R2_gen.next()

                # Add a tuple read1 and read2 to the end of the queue
                self.inq.put( ( read1, read2 ) )
                n+=1

                # update the progress bar
                self.progress_bar(n)

        except StopIteration as E:
            print(E)

        # Add a STOP pill to the queue
        for i in range(self.n_thread):
            self.inq.put("STOP")

    def filter(self, number):
        """
        Parallelized filter that take as input a sequence couple in inqueue until a STOP pill is
        found. Sequences go through a QualityFilter and a AdapterTrimmer object and ifthe couple
        is able to pass filters then it is put at the end of outqueue. at the ebd of the process
        a STOP pill is added to the outqueue.
        """
        # Consume inq and produce and fill outq
        for read1, read2 in iter(self.inq.get, "STOP"):

            with self.total.get_lock():
                self.total.value+=1

            # Quality filtering
            if self.quality_trim:
                read1 = self.quality_trimmer(read1)
                read2 = self.quality_trimmer(read2)
                if not read1 or not read2:
                    continue

                with self.pass_qual.get_lock():
                    self.pass_qual.value+=1

            # Adapter trimming
            if self.adapter_trim:
                read1 = self.adapter_trimmer(read1)
                read2 = self.adapter_trimmer(read2)
                if not read1 or not read2:
                    continue

                with self.pass_adapt.get_lock():
                    self.pass_adapt.value+=1

            # If both filters passed = add to the output queue
            self.outq.put( (read1, read2) )

        # Add a STOP pill to the queue
        self.outq.put("STOP")
        #print ("Filter N° {} done".format(number))

        # Collecting Trimmer informations
        if self.quality_trim:
            q_dict = self.quality_trimmer.get_summary()
            with self.qual_total.get_lock():
                self.qual_total.value += q_dict["total"]
            with self.qual_untrimmed.get_lock():
                self.qual_untrimmed.value+=q_dict["untrimmed"]
            with self.qual_trimmed.get_lock():
                self.qual_trimmed.value+=q_dict["trimmed"]
            with self.qual_fail.get_lock():
                self.qual_fail.value+=q_dict["fail"]
            with self.qual_base_trimmed.get_lock():
                self.qual_base_trimmed.value+=q_dict["base_trimmed"]
            with self.qual_mean_sum.get_lock():
                self.qual_mean_sum.value+=q_dict["qual_mean_sum"]

        if self.adapter_trim:
            a_dict = self.adapter_trimmer.get_summary()
            with self.adapt_total.get_lock():
                self.adapt_total.value+=a_dict["total"]
            with self.adapt_untrimmed.get_lock():
                self.adapt_untrimmed.value+=a_dict["untrimmed"]
            with self.adapt_trimmed.get_lock():
                self.adapt_trimmed.value+=a_dict["trimmed"]
            with self.adapt_fail.get_lock():
                self.adapt_fail.value+=a_dict["fail"]
            with self.adapt_base_trimmed.get_lock():
                self.adapt_base_trimmed.value+=a_dict["base_trimmed"]
            #for n, adapter in enumerate (a_dict["adapter_found"]):
                #with self.adapter_found.get_lock():
                    #self.adapter_found[n] += adapter

    def writer(self, R1_outname, R2_outname):
        """
        Write sequence couples from outqueue in a pair of fastq files. Sequences will remains
        paired (ie at the same index in the 2 files) but they may not be in the same order
        than in the input fastq files. The process will continue until n = n_thead STOP pills were
        found in the outqueue (ie. the queue is empty)
        """
        # Open output fastq streams for writing
        try:
            out_R1 = gopen(R1_outname, "wb") if self.compress_output else open(R1_outname, "wb")
            out_R2 = gopen(R2_outname, "wb") if self.compress_output else open(R2_outname, "wb")

            current_seq = 0
            buffer_R1 = ""
            buffer_R2 = ""

            # Keep running until all thread STOP pills has been passed
            for works in range(self.n_thread):
                # Will exit the loop as soon as a STOP pill will be found
                for read1, read2 in iter(self.outq.get, "STOP"):

                    with self.total_pass.get_lock():
                        self.total_pass.value+=1

                    buffer_R1 += read1.get_fastq_str()
                    buffer_R2 += read2.get_fastq_str()

                    if self.total_pass.value%self.buffer_size == 0:
                        out_R1.write(buffer_R1)
                        out_R2.write(buffer_R2)
                        buffer_R1 = ""
                        buffer_R2 = ""

            out_R1.write(buffer_R1)
            out_R2.write(buffer_R2)
            buffer_R1 = ""
            buffer_R2 = ""

            out_R1.close()
            out_R2.close()

        except IOError as e:
            print "I/O error({}): {}".format(e.errno, e.strerror)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_values(self):
        """
        Test the validity of options in the configuration file
        """
        # Verify values from the quality section
        assert self.min_size >= 0, "Authorized values for min_size : >= 0"
        assert self.n_thread > 0, "Authorized values for n_thread : > 0"

        if self.quality_trim:
            assert self.win_size > 0, "Authorized values for win_size : > 0"
            assert self.step > 0, "Authorized values for step : > 0"
            assert 0 <= self.qual_cutdown <= 40, "Authorized values for qual_cutdown : 0 to 40"

        if self.adapter_trim:
            assert self.ssw_match
            assert self.ssw_mismatch
            #assert self.ssw_ambiguous
            assert self.ssw_gapO
            assert self.ssw_gapE

            max_score = self.ssw_match
            min_score = - max([self.ssw_mismatch, self.ssw_gapO, self.ssw_gapE])

            assert 0 < self.min_match_len <= 1, "Authorized values for min_match_len : > 0 to 1"
            assert min_score <= self.min_match_score <= max_score, "Authorized values for min_match_score : - higher penalty to ssw_match"

    def _count_fastq (self, fastq):
        """
        Basic fastq line counter
        """
        try:
            fp = gopen(fastq, "rb") if fastq[-2:].lower() == "gz" else open(fastq, "rb")
            nline = 0

            for line in fp:
                nline+=1
            fp.close()

            return nline/4

        except IOError as e:
            print "I/O error({}): {}".format(e.errno, e.strerror)

    def _write_report (self, sample_name, n_adapter):

        with open ("{}_trimming_report.txt".format(sample_name), "wb") as report:
            report.write("Sample name\t{}\n".format(sample_name))
            report.write("Generic section\n")
            report.write("Total pair\t{}\n".format(self.total.value))
            report.write("Pass quality trimming\t{}\n".format(self.pass_qual.value))
            report.write("Pass adapter trimming\t{}\n".format(self.pass_adapt.value))
            report.write("Pass total\t{}\n".format(self.total_pass.value))

            # Define Quality Trimmer Object and specific shared memory counters
            if self.quality_trim:
                report.write("\nQuality trimming section\n")
                report.write("Total read\t{}\n".format(self.qual_total.value))
                report.write("Untrimmed\t{}\n".format(self.qual_untrimmed.value))
                report.write("Trimmed\t{}\n".format(self.qual_trimmed.value))
                report.write("Fail\t{}\n".format(self.qual_fail.value))
                report.write("Base trimmed\t{}\n".format(self.qual_base_trimmed.value))
                report.write("Mean quality\t{}\n".format(self.qual_mean_sum.value/self.qual_total.value))

            # Define Adapter Trimmer Object and specific shared memory counters
            if self.adapter_trim:
                report.write("\nAdapter trimming section\n")
                report.write("Total read\t{}\n".format(self.adapt_total.value))
                report.write("Untrimmed\t{}\n".format(self.adapt_untrimmed.value))
                report.write("Trimmed\t{}\n".format(self.adapt_trimmed.value))
                report.write("Fail\t{}\n".format(self.adapt_fail.value))
                report.write("Base trimmed\t{}\n".format(self.adapt_base_trimmed.value))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    sekator = Sekator.class_init()
    sekator()
