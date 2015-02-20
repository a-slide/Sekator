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

    def __init__(self, qual_cutdown=25, win_size=5, step=1, min_size=30, left_trim=True, right_trim=True):
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
        self.step = step
        self.min_size = min_size
        self.left_trim = left_trim
        self.right_trim = right_trim

        # Counters
        self.total = 0
        self.untrimmed = 0
        self.trimmed = 0
        self.fail = 0
        self.base_trimmed = 0
        self.qual_mean_sum = 0

    @property
    def mean_qual(self):
        return self.qual_mean_sum / self.total

    def __repr__(self):
        msg = "QUALITY TRIMMER CLASS\n"
        msg += "\tQuality cutdown : {}\n".format(self.qual_cutdown)
        msg += "\tSliding windows size : {}\n".format(self.win_size)
        msg += "\tSliding windows step : {}\n".format(self.step)
        msg += "\tMinimal size of sequences : {}\n".format(self.min_size)
        msg += "\tTrim from left : {}\n".format(self.left_trim)
        msg += "\tTrim from right : {}\n".format(self.right_trim)
        msg += "\tTotal sequences : {}\n".format(self.total)
        msg += "\tUntrimmed sequences : {}\n".format(self.untrimmed)
        msg += "\tTrimmed sequences : {}\n".format(self.trimmed)
        msg += "\tFailed sequence : {}\n".format(self.fail)
        msg += "\tNumber of base trimmed : {}\n".format(self.base_trimmed)
        msg += "\tCumulative sum of quality : {}\n".format(self.qual_mean_sum)
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
        self.qual_mean_sum += seq.qual.mean()
        seq_size = len(seq)
        start = 0
        end = seq_size

        # Trimming left end
        if self.left_trim:
#            print ("Left trim")

            # Loop from the begining of seq until the windows quality is high enough
            for i in np.arange(0, seq_size-self.win_size+1, self.step):

                print ("Win : {}  Qual : {}".format(seq.qual[i:i+self.win_size], seq.qual[i:i+self.win_size].mean()))
                # Mark the start and leave the loop if the quality of the windows is high enough
                if seq.qual[i:i+self.win_size].mean() >= self.qual_cutdown:
                    start = i
                    break

            # If the windows arrive at the end of the sequence return None
            self.fail += 1
            self.base_trimmed += seq_size
            return None

        # Trimming right end
        if self.right_trim:
#            print ("Right trim")

            # Back loop from the end of seq until the windows quality is high enough
            for i in np.arange(seq_size, 0+self.win_size-1, -self.step):

                print ("Win : {}  Qual : {}".format(seq.qual[i-self.win_size:i], seq.qual[i-self.win_size:i].mean()))
                # Mark the end and leave the loop if the quality of the windows is high enough
                if seq.qual[i-self.win_size:i].mean() >= self.qual_cutdown:
                    end = i
                    break

            # If the windows arrive at the beginning of the sequence return None
            self.fail += 1
            self.base_trimmed += seq_size
            return None

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
            self.fail += 1
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
