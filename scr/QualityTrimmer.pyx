#~~~~~~~IMPORTS~~~~~~~#

# Third party imports
import numpy as np
from HTSeq import SequenceWithQualities as HTSfastq

#~~~~~~~CIMPORTS~~~~~~~#

# C standard library import
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cdef class QualityTrimmer(object):

#    Read quality trimmer using a sliding window to scan the sequence starting by left and/or
#    right extremities. Invalid bases and trimmed from the returned sequence
#    If all base are invalid or if the size after trimming is bellow a minimal size, None is returned

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~SELF VARIABLES DEFINITION~~~~~~~#

    cdef:
        # Init object variables
        uint32_t qual_cutdown, win_size, step, min_size, total, untrimmed, trimmed, fail
        uint64_t base_trimmed, qual_mean_sum
        bint left_trim, right_trim
    
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
    
    def __init__(self, uint32_t qual_cutdown=25, uint32_t win_size=5, uint32_t step=1,\
        uint32_t min_size=30, bint left_trim=True, bint right_trim=True):

#       Init quality trimmer
#       @param qual_cutdown Minimal quality in a given windows
#       @param win_size Size of the sliding windows
#       @param step step of the sliding window trimming 
#       @param min_size Minimal size of read to be considered as valid
#       @param left_trim Triming starting from left extremity of reads
#       @param right_trim Triming starting from right extremity of reads

        # Init object variables
        self.qual_cutdown = qual_cutdown
        self.win_size = win_size
        self.left_trim = left_trim
        self.right_trim = right_trim
        self.min_size = min_size
        self.step = step

        # Counters
        self.total = 0
        self.untrimmed = 0
        self.trimmed = 0
        self.fail = 0
        self.base_trimmed = 0
        self.qual_mean_sum = 0


#    def __repr__(self):
#        msg = "QUALITY TRIMMER CLASS\n"
#        # list all values in object dict in alphabetical order
#        keylist = [key for key in self.__dict__.keys()]
#        keylist.sort()
#        for key in keylist:
#            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
#        return (msg)

    def __str__(self):
        return ("<Instance of QualityTrimmer Class>\n")

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, seq):

#        Compute mean quality score and compare to the minimal quality required
#        @param seq HTSeq.SequenceWithQualities object

        cdef:
            int32_t seq_size, start, end, i
        
        # Update counters and init border index
        self.total += 1
        self.qual_mean_sum += seq.qual.mean()
        seq_size = len(seq)
        start = 0
        end = seq_size

        # Trimming left end
        if self.left_trim:

            # Loop from the begining of seq until the windows quality is high enough
            for i in np.arange(start=0, stop=seq_size, step=self.step):

                # Extract quality windows
                win_qual = seq.qual[i:i+self.win_size]

                # If the windows arrive at the end of the sequence return None
                if win_qual.size < self.win_size:
                    self.qual_fail += 1
                    self.base_trimmed += seq_size
                    return None

                # Mark the start and leave the loop if the quality of the windows is high enough
                if win_qual.mean() >= self.qual_cutdown:
                    start = i
                    break

        # Trimming right end
        if self.right_trim:

            # Back loop from the end of seq until the windows quality is high enough
            for i in np.arange(start=seq_size, stop=0, step=-self.step):

                # Extract quality windows
                win_qual = seq.qual[i-self.win_size:i]

                # If the windows arrive at the end of the sequence return None
                if win_qual.size < self.win_size:
                    self.qual_fail += 1
                    self.base_trimmed += seq_size
                    return None

                # Mark the end and leave the loop if the quality of the windows is high enough
                if win_qual.mean() >= self.qual_cutdown:
                    end = i
                    break
        
        # In the case were no trimming was done
        if start == 0 and end == seq_size:
            self.untrimmed += 1
            return seq
        
        # Return the trimmed read if its lenghth is sufficient
        if end-start >= self.min_size:
            self.trimmed += 1
            self.base_trimmed += (start + seq_size - end)

            ## Create a new object instead of slicing
            return HTSfastq(
                seq=seq.seq[start:end],
                name=seq.name,
                qualstr=seq.qualstr[start:end])

        else:
            self.qual_fail += 1
            self.base_trimmed += seq_size
            return None

    def get_summary (self):
            
        summary = {}
        summary["total"] = int(self.total)
        summary["untrimmed"] = int(self.untrimmed)
        summary["trimmed"] = int(self.trimmed)
        summary["fail"] = int(self.fail)
        summary["base_trimmed"] = int(self.base_trimmed)
        summary["qual_mean_sum"] = int(self.qual_mean_sum)
        
        return summary
