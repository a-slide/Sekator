# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      Contain the template of the empty configuration file for Sekator
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

def write_example_conf():

    with open ("Sekator_conf_file.txt", 'wb') as fp:
        fp.write ("""
###################################################################################################
#                                  SEKATOR CONFIGURATION FILE                                     #
###################################################################################################
# Values can by customized with users values, but the file template must remain unchanged,
# otherwise the program will not be able to load default values.
# File path should be indicated as absolute path preferably and should not contain blank spaces
# Values identified with '**' in the descriptor are not recommended to be modified

###################################################################################################
[general]
# The quality encoding of your sequence have to be Illumina 1.8+ Phred+33. The program does not
# manage the other encoding scales

# Minimal size of read after trimming (POSITIVE INTEGER)
min_size = 30

# Use all available threads for parrallel processing (BOOLEAN)
auto_thread = True

# If auto_thread is "False" specify the maximal number of thread to use (POSITIVE INTEGER)
n_thread :

# Write a txt report (BOOLEAN)
write_report : True

# Compress the fastq output (BOOLEAN)
compress_output : True

###################################################################################################
[quality]

# Perform quality trimming from the left and/or extremities (BOOLEAN)
left_trim : True
right_trim : True

# Size of the sliding window in which quality will be computed (POSITIVE INTEGER)
win_size : 6

# Step of sliding window during trimming (POSITIVE INTEGER)
step : 2

# Minimal quality in a given windows to be retained during trimming (0 <= POSITIVE INTEGER <= 40)
qual_cutdown : 28

###################################################################################################
[adapter]

# Perform a step of adapter trimming (BOOLEAN)
adapter_trim : True

# Minimal fraction of the length of the adapter matching on the read (0 < FLOAT <= 1) **
min_match_len : 0.3

# Minimal SSW score/base of the adapter matching on the read (POSITIVE FLOAT <= ssw_match) **
min_match_score : 1

# Scores for stripped Smith and Waterman sequence alignment if :
# - Gain if 2 aligned bases are identical : ssw_match
# - Penalty if 2 aligned bases are different : ssw_mismatch
# - Penalty if a gap is opened : ssw_gapO
# - Penalty if a gap is extended : ssw_gapE
# All values have to be POSITIVE INTEGER. Penalty will be converted in negative scores by the
# aligner. Values were optimized for adapter match (40-60 pb) on short reads (100, 200 pb)  **

ssw_match : 2
ssw_mismatch : 2
ssw_gapO : 3
ssw_gapE : 1

###################################################################################################
# SAMPLE DEFINITIONS

# It is possible to include as many independant sample as required by duplicating a entire sample
# section and incrementing the id number in the sample section name
# Each sample section is organize as follow :
# name = Unique identifier that will be used to prefix the read files (STRING)
# - R1_path = Valid path to the fastq(.gz) file containing the forward reads of the pair (STRING)
# - R2_path = Valid path to the fastq file (gziped or not) containing the reverse reads of the pair
#   preferably absolute path without spaces) (STRING)
# - adapter_list = list of adapter DNA sequence to be trimmed if adapter_trimming is required.
# - Separate each adapter by a blank space (LIST OF STR)

[sample1]
name : S1_pass
R1_path : ../dataset/S1_R1_pass.fastq.gz
R2_path : ../dataset/S1_R2_pass.fastq.gz
adapter_list : GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT

[sample2]
name : test
R1_path : ../dataset/testR1.fastq.gz
R2_path : ../dataset/testR2.fastq.gz
adapter_list : GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  """)
