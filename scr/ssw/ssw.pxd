from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t

cdef extern from "ssw.c":

    ctypedef struct s_align:
        uint16_t score
        int32_t ref_begin
        int32_t ref_end
        int32_t query_begin
        int32_t query_end

    # Main functions
    int8_t* score_matrix (
        const int8_t match,
        const int8_t mismatch,
        const int8_t ambiguous)

    int8_t* DNA_seq_to_int(
        char* DNA_seq,
        int seq_len)

    s_align ssw_align (
        const int8_t* query,
        int32_t queryLen,
        const int8_t* ref,
        int32_t refLen,
        const int8_t* mat,
        const int8_t gapO,
        const int8_t gapE)

    #Memory management functions
    void print_int8_t (int8_t* mat, int dim)
