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
    from multiprocessing import Value, Process, Queue, cpu_count
    from time import time
    import gzip
    from os import path
    import ConfigParser
    import optparse
    import sys
    import os
    import gzip

    # Third party imports
    import numpy
    from HTSeq import FastqReader

    # Local Package import
    #from pyDNA.Utilities import count_seq
    from AdapterTrimmer import AdapterTrimmer
    from QualityTrimmer import QualityTrimmer

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sekator (object):
    """
    Sekator is a multiprocessing fastq file trimmer, handling bot quality trimming and
    adaptor sequence trimming
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, conf_file = None):
        """
        Initialization function, parse options in command line and configuration file and verify
        their values. All self.variables are initialized explicitly in init
        conf_file arg is mandatory in interactive interpreter and import, else optparse will parse
        command line arguments
        """

        print("Initialize Sekator")

        # Define fundamental variables
        self.version = "Sekator 0.1"
        self.usage = "Usage: %prog -c Conf.txt"

        # Use Conf file if interactive interpreter else parse arguments with optparse
        self.conf = conf_file if conf_file else self._optparser()

        # Parse the configuration file and verify the values of variables
        try:
            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=False)
            cp.read(self.conf)

            # General section
            self.qual_scale = cp.get("general", "qual_scale")
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
                self.find_reverse = cp.getboolean("adapter", "find_reverse")
                self.min_match_len = cp.getfloat("adapter", "min_match_len")
                self.min_match_score = cp.getfloat("adapter", "min_match_score")
                self.ssw_match = cp.getint("adapter", "ssw_match")
                self.ssw_mismatch = cp.getint("adapter", "ssw_mismatch")
                self.ssw_ambiguous = cp.getint("adapter", "ssw_ambiguous")
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

        """
        # Init shared memory counters
        self.total = Value('i', 0)
        self.pass_qual = Value('i', 0)
        self.pass_trim = Value('i', 0)
        self.total_pass = Value('i', 0)
        if self.quality_trim:
            self.min_qual_found = Value('i', 100)
            self.max_qual_found = Value('i', 0)
            self.weighted_mean = Value('d', 0.0)
        if self.adapter_trim:
            self.seq_untrimmed = Value('i', 0)
            self.seq_trimmed = Value('i', 0)
            self.base_trimmed = Value('i', 0)
            self.len_pass = Value('i', 0)
            self.len_fail = Value('i', 0)
        """

    def __repr__(self):
        msg = "SEKATOR CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)


    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def __call__ (self):
        """
        Main function of the script
        """
        # Start a timer
        start_time = time()

        for sample in self.sample_list:

            if self.quality_filter:
                qt = QualityTrimmer(
                    qual_cutdown = self.qual_cutdown,
                    win_size = self.win_size,
                    step = self.step,
                    min_size = self.min_size,
                    left_trim = self.left_trim,
                    right_trim = self.right_trim)

            if self.adapter_trim:
                at = AdapterTrimmer(
                    adapter_list = sample.adapter_list,
                    find_reverse = self.find_reverse,
                    min_size = self.min_size,
                    min_match_len = self.min_match_len,
                    min_match_score = self.min_match_score,
                    ssw_match = self.ssw_match,
                    ssw_mismatch = self.ssw_mismatch,
                    ssw_ambiguous = self.ssw_ambiguous,
                    ssw_gapO = self.ssw_gapO,
                    ssw_gapE = self.ssw_gapE)

        ## Count lines in fastq file to prepare a counter of progression
        #print ("Count the number of fastq sequences")
        #self.nseq = count_seq(R1, "fastq")
        #print("fastq files contain {} sequences to align".format(self.nseq))
        #self.nseq_list = [int(self.nseq*i/100.0) for i in range(5,101,5)] # 5 percent steps

        ## Init queues for input file reading and output file writing (limited to 10000 objects)
        #self.inq = Queue(maxsize=10000)
        #self.outq = Queue(maxsize=10000)

        ## Init processes for file reading, distributed filtering and file writing
        #self.pin = Process(target=self.reader, args=())
        #self.ps = [Process(target=self.filter, args=()) for i in range(self.numprocs)]
        #self.pout = Process(target=self.writer, args=())

        ## Start processes
        #self.pin.start()
        #self.pout.start()
        #for p in self.ps:
            #p.start()

        ## Blocks until the process is finished
        #self.pin.join()
        #print ("\tReading done")
        #for i in range(len(self.ps)):
            #self.ps[i].join()
        #print ("\tFiltering done")
        #self.pout.join()
        #print ("\tWriting done\n")

        ## Stop timer and store the value
        #self.exec_time = round(time()-start_time, 3)

    #def reader(self):
        #"""
        #Initialize SeqIO.parse generators to iterate over paired fastq files. Data ara sent over
        #inqueue for the workers to do their thing and a n = numprocs STOP pills are added at the
        #end of the queue for each worker.
        #"""
        #try:
            ## Open input fastq streams for reading
            #if self.R1_in[-2:].lower() == "gz":
                #in_R1 = gzip.open(self.R1_in, "rb")
            #else:
                #in_R1 = open(self.R1_in, "rb")

            #if self.R2_in[-2:].lower() == "gz":
                #in_R2 = gzip.open(self.R2_in, "rb")
            #else:
                #in_R2 = open(self.R2_in, "rb")

        #except (IOError, TypeError, ValueError) as E:
            #print E
            #exit

        ## Init generators to iterate over files
        #genR1 = SeqIO.parse(in_R1, self.qual_scale)
        #genR2 = SeqIO.parse(in_R2, self.qual_scale)

        #i = 0
        #while True:
            ## Parse sequences in generators until one of then is empty
            #seqR1 = next(genR1, None)
            #seqR2 = next(genR2, None)
            #if not seqR1 or not seqR2:
                #break
            ## Add a tuple position, seqR1 and seqR2 to the end of the queue
            #self.inq.put( (seqR1, seqR2) )

            #i+=1
            #if i in self.nseq_list:
                #print ("\t{} sequences: {}%".format(i, int(i*100.0/self.nseq)))

        ## Close files
        #in_R1.close()
        #in_R2.close()

        ## Add a STOP pill to the queue
        #for i in range(self.numprocs):
            #self.inq.put("STOP")

    #def filter(self):
        #"""
        #Parallelized filter that take as input a sequence couple in inqueue until a STOP pill is
        #found. Sequences go through a QualityFilter and a AdapterTrimmer object and ifthe couple
        #is able to pass filters then it is put at the end of outqueue. at the ebd of the process
        #a STOP pill is added to the outqueue.
        #"""
        ## Consume inq and produce answers on outq
        #for seqR1, seqR2 in iter(self.inq.get, "STOP"):

            #with self.total.get_lock():
                #self.total.value+=1

            ## Quality filtering
            #if self.qual:
                #seqR1 = self.qual.filter(seqR1)
                #seqR2 = self.qual.filter(seqR2)
                #if not seqR1 or not seqR2:
                    #continue

            #with self.pass_qual.get_lock():
                #self.pass_qual.value+=1

            ## Adapter trimming and size filtering
            #if self.adapt:
                #seqR1 = self.adapt.trimmer(seqR1)
                #seqR2 = self.adapt.trimmer(seqR2)
                #if not seqR1 or not seqR2:
                    #continue

            #with self.pass_trim.get_lock():
                #self.pass_trim.value+=1

            ## If both filters passed = add to the output queue
            #self.outq.put( (seqR1, seqR2) )

        ## Add a STOP pill to the queue
        #self.outq.put("STOP")

        ## Fill shared memomory counters from process specific object instances.
        #if self.qual:
            #with self.weighted_mean.get_lock():
                #self.weighted_mean.value += (self.qual.get_mean_qual()*self.qual.get('total'))
            #if self.qual.get_min_qual() < self.min_qual_found.value:
                #self.min_qual_found.value = self.qual.get_min_qual()
            #if self.qual.get_max_qual() > self.max_qual_found.value:
                #self.max_qual_found.value = self.qual.get_max_qual()

        #if self.adapt:
            #with self.seq_untrimmed.get_lock():
                #self.seq_untrimmed.value += self.adapt.get('seq_untrimmed')
            #with self.seq_trimmed.get_lock():
                #self.seq_trimmed.value += self.adapt.get('seq_trimmed')
            #with self.base_trimmed.get_lock():
                #self.base_trimmed.value += self.adapt.get('base_trimmed')
            #with self.len_pass.get_lock():
                #self.len_pass.value += self.adapt.get('len_pass')
            #with self.len_fail.get_lock():
                #self.len_fail.value += self.adapt.get('len_fail')

    #def writer(self):
        #"""
        #Write sequence couples from outqueue in a pair of compressed fastq.gz files. Sequences will
        #remains paired (ie at the same index in the 2 files) but they may not be in the same order
        #than in the input fastq files. The process will continue until n = numprocs STOP pills were
        #found in the outqueue (ie. the queue is empty)
        #"""
        ## Open output fastq streams for writing
        #if self.compress_output:
            #out_R1 = gzip.open(self.R1_out, "wb")
            #out_R2 = gzip.open(self.R2_out, "wb")
        #else:
            #out_R1 = open(self.R1_out, "wb")
            #out_R2 = open(self.R2_out, "wb")

        ## Keep running until all numprocs STOP pills has been passed
        #for works in range(self.numprocs):
            ## Will exit the loop as soon as a Stop pill will be found
            #for seqR1, seqR2 in iter(self.outq.get, "STOP"):
                #out_R1.write(seqR1.format("fastq-sanger"))
                #out_R2.write(seqR2.format("fastq-sanger"))
                #with self.total_pass.get_lock():
                    #self.total_pass.value+=1

        #out_R1.close()
        #out_R2.close()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self):
        """
        Parse CLI and return a valid configuration file path
        """
        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.usage, version = self.version)
        optparser.add_option('-c', dest="conf_file", help= "Path to the configuration file")

        # Parse arguments
        options, args = optparser.parse_args()

        # Verify conf file
        if not options.conf_file:
            optparser.print_help()
            optparser.error("incorrect number of arguments")

        return options.conf_file

    def _test_values(self):
        """
        Test the validity of options in the configuration file
        """
        # Verify values from the quality section
        assert self.qual_scale in ["solexa", "solexa-old", 'phred'], "Authorized values for quality_scale : solexa, solexa-old, phred"
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    SAMPLE_NAMES = []

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def ADD_TO_SAMPLE_NAMES(self, name):
        self.SAMPLE_NAMES.append(name)

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, R1_path, R2_path, adapter_list, compress_output):

        # Create self variables
        self.name = name
        self.R1_path = R1_path
        self.R2_path = R2_path
        self.adapter_list = adapter_list

        self._test_values()

        self.R1_outname = "{}_R1_filtered.fastq{}".format(self.name, ".gz" if compress_output else "")
        self.R2_outname = "{}_R2_filtered.fastq{}".format(self.name, ".gz" if compress_output else "")
        self.ADD_TO_SAMPLE_NAMES(self.name)

    # Fundamental class methods str and repr
    def __repr__(self):
        msg = "SAMPLE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_values(self):
        assert self.name not in self.SAMPLE_NAMES, "Sample name <{}> is duplicated".format(self.name)
        assert self._is_readable_file (self.R1_path), "R1_path in Sample <{}> is not valid".format(self.name)
        assert self._is_readable_file (self.R2_path), "R2_path in Sample <{}> is not valid".format(self.name)
        for adapter in self.adapter_list:
            assert self._is_dna(adapter), "<{}> in Sample <{}> is not a valid DNA sequence".format(adapter, self.name)

    def _is_readable_file (self, fp):
        return os.access(fp, os.R_OK)

    def _is_dna (self, sequence):
        for base in sequence:
            if base not in ["A","T","C","G","N","a","t","c","g","n"]:
                return False
        return True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    sekator = Sekator()
    sekator()
