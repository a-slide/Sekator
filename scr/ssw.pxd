from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t

cdef extern from "ssw.c":

    ctypedef struct s_align:
        uint16_t score
        int32_t ref_begin
        int32_t ref_end
        int32_t query_begin
        int32_t query_end

#        @typedef    structure of the alignment result
#        @field  score       the best alignment score
#        @field  ref_begin   0-based best alignment beginning position on reference; ref_begin = -1 when the best alignment beginning
#                            position is not available
#        @field  ref_end     0-based best alignment ending position on reference
#        @field  read_begin  0-based best alignment beginning position on read; read_begin = -1 when the best alignment beginning
#                            position is not available
#        @field  read_end    0-based best alignment ending position on read


    int8_t* score_matrix (const int8_t match, const int8_t mismatch, const int8_t ambiguous)

#    @function   Generate a 1D integer array representing a 5*5 score matrix from the provided match, mismatch and ambiguous penalties
#    @param match        Bonus in case of match (POSITIVE)
#    @param mismatch     Penalty in case of mismatch (POSITIVE)
#    @param ambiguous    Value in case of ambiguous base (Should remain NULL)
#    @return 8 bits integer array encoding a score matrix (see bellow)
#    @note   If the penalty for match, mismatch and ambiguous are 2, 2 and 0, respectively, the score matrix will be:
#            A  C  G  T  N (or other ambiguous code)
#             2 -2 -2 -2  0   A
#            -2  2 -2 -2  0   C
#            -2 -2  2 -2  0   G
#            -2 -2 -2  2  0   T
#             0  0  0  0  0   N (or other ambiguous code)
#            The corresponding array is {2,-2,-2,-2,0,-2,2,-2,-2,0,-2,-2,2,-2,0,-2,-2,-2,2,0,0,0,0,0,0}


    int8_t* DNA_seq_to_int(char* DNA_seq, const int32_t seq_len)

#    @function       Transform a DNA sequence in a int8_t* integer array with A,a=0; C,c=1; G,g=2; T,t=3; other char = 4
#    @param  seq_len Length of the DNA sequence
#    @param  DNA_seq DNA sequence in char
#    @return 8 bits integer array encoding the sequence
#    @note   Example if the query sequence is: ACGTN, the array will be : {0, 1, 2, 3, 4}

    s_align ssw_align (const int8_t* query, int32_t queryLen, const int8_t* ref, int32_t refLen, const int8_t* mat, const int8_t gapO, const int8_t gapE)

#    @function   Do Striped Smith-Waterman alignment.
#    @param query    Pointer to the query sequence; the query sequence needs to be numbers
#    @param queryLen Length of the query sequence
#    @param ref      Pointer to the reference sequence; the reference sequence needs to be numbers
#    @param refLen   Length of the reference sequence
#    @param mat      Score matrix produced by score_matrix
#    @param gapO     Penalty in case of gap opening (POSITIVE)
#    @param gapE     Penalty in case of gap extension (POSITIVE)
#    @return Alignment result structure
