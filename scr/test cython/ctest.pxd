from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t

cdef extern from "test.h":

    ctypedef struct align:
        uint16_t    score1
        uint16_t    score2
        int32_t     ref_begin1
        int32_t     ref_end1
        int32_t     read_begin1
        int32_t     read_end1
        int32_t     ref_end2
        int32_t     cigarLen

cdef extern from "test.c":

    int8_t* score_matrix (const uint8_t match, const uint8_t mismatch, const uint8_t ambiguous)
    int8_t* DNA_seq_to_int(char * DNA_seq, int seq_len)
    void print_int8_t (int8_t* mat, int dim)
    void free_int8_t (int8_t* mat)
