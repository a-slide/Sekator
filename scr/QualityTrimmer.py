#~~~~~~~IMPORTS~~~~~~~#

# Third party imports
import numpy as np
from HTSeq import SequenceWithQualities as HTSfastq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class QualityTrimmer(object):
    """
    Read quality trimmer using a sliding window to scan the sequence starting by left and/or
    right extremities. Invalid bases and trimmed from the returned sequence
    If all base are invalid or if the size after trimming is bellow a minimal size, None is returned
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, qual_cutdown=25, win_size=5, step=1, left_trim=True, right_trim=True, min_size=1):
        """
        Init quality trimmer
        @param qual_cutdown Minimal quality in a given windows
        @param win_size Size of the sliding windows
        @param left_trim Triming starting from left extremity of reads
        @param right_trim Triming starting from right extremity of reads
        @param min_size Minimal size of read to be considered as valid
        """
        # Init object variables
        self.qual_cutdown = qual_cutdown
        self.win_size = win_size
        self.left_trim = left_trim
        self.right_trim = right_trim
        self.min_size = min_size
        self.step = step

        # Counters
        self.total = 0
        self.qual_pass = 0
        self.qual_fail = 0
        self.cumulative_sum = 0
        self.base_trimmed = 0

    @property
    def mean_qual(self):
        return self.cumulative_sum / self.total

    def __repr__(self):
        msg = "QUALITY TRIMMER CLASS\n"
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
        @param seq HTSeq.SequenceWithQualities object
        """

        # Update counters and init border index
        self.total += 1
        self.cumulative_sum += seq.qual.mean()
        start = 0
        end = seq.qual.size

        # Trimming left end
        if self.left_trim:

            # Loop from the begining of seq until the windows quality is high enough
            for i in np.arange(start=0, stop=seq.qual.size, step=self.step):

                # Extract quality windows
                win_qual = seq.qual[i:i+self.win_size]

                # If the windows arrive at the end of the sequence return None
                if win_qual.size < self.win_size:
                    self.qual_fail += 1
                    self.base_trimmed += end
                    return None

                # Mark the start and leave the loop if the quality of the windows is high enough
                if win_qual.mean() >= self.qual_cutdown:
                    start = i
                    break

        # Trimming right end
        if self.right_trim:

            # Back loop from the end of seq until the windows quality is high enough
            for i in np.arange(start=seq.qual.size, stop=0, step=-self.step):

                # Extract quality windows
                win_qual = seq.qual[i-self.win_size:i]

                # If the windows arrive at the end of the sequence return None
                if win_qual.size < self.win_size:
                    self.qual_fail += 1
                    self.base_trimmed += end
                    return None

                # Mark the end and leave the loop if the quality of the windows is high enough
                if win_qual.mean() >= self.qual_cutdown:
                    end = i
                    break

        # Return the record if its quality is high enough after trimming
        if end-start >= self.min_size:
            self.qual_pass += 1
            self.base_trimmed += (start + seq.qual.size - end)

            ## Create a new object instead of slicing
            return HTSfastq(
                seq=seq.seq[start:end],
                name=seq.name,
                qualstr=seq.qualstr[start:end])

        else:
            self.qual_fail += 1
            self.base_trimmed += end
            return None
