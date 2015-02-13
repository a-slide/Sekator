#~~~~~~~GLOBAL IMPORTS~~~~~~~#

import numpy as np
import HTSeq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class QualityFilter(object):
    """
    Simple quality filtering of fastq reads based on the overall quality of reads. If bellow the
    threshold no read will be returned
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, qual_cutdown, win_size, left_trim, right_trim, min_size):
        """
        Simple storage of object variables and counters init
        """
        # Init object variables
        self.qual_cutdown = qual_cutdown
        self.win_size = win_size
        self.left_trim = left_trim
        self.right_trim = right_trim
        self.min_size = min_size
        
        # Counters
        self.total = 0
        self.qual_pass = 0
        self.qual_fail = 0
        self.cumulative_sum = 0

    @property
    def mean_qual(self):
        return self.cumulative_sum / self.total
    
    def __repr__(self):
        msg = "Quality filer CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, seq):
        """
        Compute mean quality score and compare to the minimal quality required
        seq is and HTSeq.SequenceWithQualities object
        """
        
        self.total += 1
        self.cumulative_sum += seq.qual.mean()
        
        if self.left_trim:
            seq = self._left_trim(seq)
        
        if self.right_trim
            seq = self._right_trim(seq)
        
        # Return the record if its quality is high enough after trimming
        if len(seq) < self.min_size:
            return None
            self.qual_fail += 1
            
        else:
            self.qual_pass += 1
            return seq

    def _right_trim (self, seq):
        
        for i in np.arange(a.size, 0, -1):
            print (i)
            if a[i:i+3].mean() < 20:
                print "STOP"
                break
        
        
        

    def _left_trim (self, seq):
        
        

