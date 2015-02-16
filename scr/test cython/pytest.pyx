from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t

from ctest cimport align, score_matrix, DNA_seq_to_int, print_int8_t, free_int8_t

cpdef align align_struct (int b,int c,int d,int e,int f,int g,int h,int i):

    cdef align a

    a.score1 = b
    a.score2 = c
    a.ref_begin1 = d
    a.ref_end1 = e
    a.read_begin1 = f
    a.read_end1 = g
    a.ref_end2 = h
    a.cigarLen = i

    return a

cpdef int map (int match, int mismatch, int ambiguous, char* DNA):

    cdef int8_t *mat
    cdef int8_t *seq

    mat = score_matrix (match, mismatch, ambiguous)
    print_int8_t (mat, 25)

    seq = DNA_seq_to_int(DNA, len(DNA))
    print_int8_t (seq, len(DNA))

    free_int8_t (mat)
    free_int8_t (seq)

    return 1


# Example of C class

cdef class Query:

    # Declare self variables
    cdef:
        int32_t id, seq_len
        char *seq
        int8_t *score_mat
        int8_t *seq_int
        int8_t, int8_t mismatch, ambiguous
        s_profile *profile

    # Fast init method
    def __cinit__(self, int32_t id, char* seq, int match, int mismatch, int ambiguous):

        self.id = id
        self.seq = seq
        self.match = <int8_t>match
        self.mismatch = <int8_t>mismatch
        self.ambiguous = <int8_t>ambiguous
        self.seq_len = <int32_t>len(seq)
        self.seq_int = DNA_seq_to_int(self.seq, self.seq_len)
        #self.profile = ssw_init (self.seq_int, self.seq_len, self.score_mat, 5, 2)

    def q_print (self):
        print ("id {}\tlen{}".format(self.id, self.seq_len))
        print (self.seq)
        print_int8_t(self.seq_int, self.seq_len)
        print("match {}\tmismatch {}\tambiguous {}".format(self.match, self.mismatch, self.ambiguous))
        print_int8_t(self.score_mat, 25)


    def __dealloc__(self):
        print ("Dealloc matrix seq id : {}".format(self.id))
        score_mat_destroy (self.score_mat)
        seq_mat_destroy (self.seq_int)
