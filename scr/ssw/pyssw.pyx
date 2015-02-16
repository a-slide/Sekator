# Advanced integer type import
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, free

# Local ssw import
from ssw cimport s_align, score_matrix, DNA_seq_to_int, ssw_align, print_int8_t

# Define a C-struct to store adapter information

ctypedef struct s_query:
    int32_t size
    char* seq
    int8_t* seq_int
    int32_t found
    int32_t min_len
    int32_t min_score

cdef class AdapterTrimmer:

    # Define self object variables here
    cdef:
        int min_size, n_query
        int8_t ssw_match, ssw_mismatch, ssw_ambiguous, ssw_gapO, ssw_gapE
        int8_t* score_mat
        s_query* query_list

    def __init__(self,
        list adapter_list, # List of adapter sequence to match on all read to be analysed
        bint find_reverse = True, # If true will also search for the reverse complementary sequence of the adapter
        int min_size = 0, # Minimal size of read after trimming to be considered valid
        float min_match_len = 0.8, # Minimal fraction of adapter len that needs to be aligned on the target
        float min_match_score = 1.4, # Minimal score per base for the alignment of adapter and read
        int8_t ssw_match = 2, # Bonus in case of match (POSITIVE)
        int8_t ssw_mismatch = -2, # Malus in case of mismatch (NEGATIVE)
        int8_t ssw_ambiguous = 0, # Value in case of ambiguous base (POSITIVE, NEGATIVE or NULL)
        int8_t ssw_gapO = -3, # Malus in case of gap opening (NEGATIVE)
        int8_t ssw_gapE = -1): # Malus in case of gap extension (NEGATIVE)

        # Store self value for future usage
        self.min_size = min_size
        self.ssw_match = ssw_match
        self.ssw_mismatch = ssw_mismatch
        self.ssw_ambiguous = ssw_ambiguous
        self.ssw_gapO = ssw_gapO
        self.ssw_gapE = ssw_gapE

        # Init a score matrix
        print("Create Score Matrix")
        score_mat = score_matrix (ssw_match, ssw_mismatch, ssw_ambiguous)
        print_int8_t(score_mat, 25)

        # init a list of adapters
        self.n_query = len(adapter_list)
        self.query_list = <s_query *>malloc(self.n_query * sizeof(s_query))

        for n, seq in enumerate(adapter_list):
            self.query_list[n] = self.build_query (seq, min_match_len, min_match_score)


    cdef s_query build_query (self, char* seq, float min_match_len, float min_match_score):

        # Create a query struct and fill the fields
        cdef s_query q

        q.size = len(seq)
        q.seq = seq
        q.seq_int = DNA_seq_to_int(seq, q.size)
        q.found = 0
        q.min_len = <int32_t>(min_match_len*q.size) #################### Cast a verif
        q.min_score = <int32_t>(min_match_score*q.size)################# Cast a verif

        print(q)

        return q


    def align(self, ref_seq):
        cdef:
            int32_t ref_size
            int8_t* ref_int
            s_align res

        # Prepare the reference sequence
        ref_size = len(ref_seq)
        ref_int = DNA_seq_to_int(ref_seq, ref_size)

        print("Ref Seq")
        print_int8_t(ref_int, ref_size)

        # Iterate over the
        for i in range(self.n_query):
            res = ssw_align(
                query = self.query_list[i].seq_int,
                queryLen = self.query_list[i].size,
                ref = ref_int,
                refLen = ref_size,
                mat = self.score_mat,
                gapO = self.ssw_gapO,
                gapE = self.ssw_gapE)

            print res

        free(ref_int)
        return 1

    def __dealloc__(self):
        print ("Dealloc memory")
        # Free memory
        for i in range(self.n_query):
            free(self.query_list[i].seq_int)
            free(self.query_list[i].seq)

        free(self.query_list)
        free(self.score_mat)
