# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      Contain a class to model a fastq sequence and an iterator function to read fastq files
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from gzip import open as gopen

# Third party imports
import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqSeq (object):
    """
    Simple Representation of a fastq file. The object support slicing and addition operations
    The quality score is a numpy array to facilitate further data manipulation
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~FUNDAMENTAL AND MAGIC METHODS~~~~~~~#

    def __init__ (self, name, seq, qual, descr=""):
        """
        @param name Name of the sequence (without spaces)
        @param seq DNA sequence string
        @param qual nunpy array of the Phred Quality of bases
        @param descr Facultative description
        """
        self.name = name
        self.seq = seq
        self.qual = qual
        self.descr = descr

        assert len(seq) == len(qual), "Sequence length and quality string length are not equal."

    def __len__ (self):
        return len(self.seq)

    def __repr__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    def __str__(self):
        if len(self) > 20:
            return "{} : {}...{}  {}...{}".format(self.name, self.seq[:10], self.seq[-10:], self.qual[:10], self.qual[-10:])
        else:
            return "{} : {}  {}".format(self.name, self.seq, self.qual)

    def __getitem__( self, item ):
        return FastqSeq(name = self.name, seq = self.seq[ item ], qual = self.qual[ item ])

    @property
    def qualstr( self ):
        qualstr_phred = ""
        for i in self.qual:
            qualstr_phred += str(unichr(i+33))
        return qualstr_phred

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_fastq_str(self):
        """ Return string formated as a fastq sequence """
        return "@{}\n{}\n+\n{}\n".format(self.name, self.seq, self.qualstr)

    def __add__(self, other):
        return FastqSeq(
            name = "{}_{}".format(self.name, other.name),
            seq = self.seq+other.seq,
            qual = np.concatenate((self.qual, other.qual)),
            descr = self.descr+other.descr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def is_gz(fp):
    """ Indicate if a file is gziped """
    return fp[-2:].lower() == "gz"

def FastqReader (fastq_file):
    """ Simple fastq reader returning a generator over a fastq file """
    try:

        # Open the file depending of the compression status
        fastq = gopen(fastq_file, "rb") if is_gz(fastq_file) else open(fastq_file, "rb")
        i=0

        # Iterate on the file until the end
        while True:

            # Extract informations from the fastq file
            name = fastq.readline()
            seq = fastq.readline()
            sep = fastq.readline()
            qual = fastq.readline()

            # Condition of loop exit
            if not name or not seq or not sep or not qual:
                break

            # Try to generate a valid FastqSeq object
            try:
                yield FastqSeq(
                name = name.rstrip()[1:].split()[0],
                seq = seq.rstrip(),
                qual = np.array([ord(x)-33 for x in qual.rstrip()]))
                i+=1

            except AssertionError as E:
                print(E)
                print ("Skipping the sequence")

        raise StopIteration("\t{} sequences parsed".format(i))

    except IOError as E:
        print(E)
        exit()
