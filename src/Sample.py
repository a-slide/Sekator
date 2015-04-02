# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      Helper class for Sekator to represent Samples
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

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
    def __str__(self):
        msg = "SAMPLE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_values(self):
        assert self.name not in self.SAMPLE_NAMES, "Sample name <{}> is duplicated".format(self.name)
        for adapter in self.adapter_list:
            assert self._is_dna(adapter), "<{}> in Sample <{}> is not a valid DNA sequence".format(adapter, self.name)

    def _is_dna (self, sequence):
        for base in sequence:
            if base not in ["A","T","C","G","N","a","t","c","g","n"]:
                return False
        return True
