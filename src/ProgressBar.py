# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      A simple ProgressBar for Sekator
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ProgressBar (object):
    """ Simple progress bar """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, total_seq, number_step):
        """
        To be init with the total number of sequence and the desired number of steps
        """
        self.total_seq = total_seq # To match 0 base index
        self.number_step = number_step
        self.numeric_step = int(self.total_seq/self.number_step) # Non exact steps
        self.n_step = 1

        assert total_seq >= number_step

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__ (self, n):
        """
        Call each iteration of the loop to verify is the progress bar needs to be updated
        """
        if n%self.numeric_step == 0:
            if self.n_step == self.number_step :
                print("\t[{}] 100% DONE".format("X"*self.n_step))
            else:
                print("\t[{}{}] {}%".format(
                "X"*self.n_step,
                "-"*(self.number_step - self.n_step),
                self.n_step*100/self.number_step))

            self.n_step +=1
