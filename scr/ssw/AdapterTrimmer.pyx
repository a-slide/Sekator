# Advanced integer type import
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, free

# Local ssw import
from ssw cimport s_align, score_matrix, DNA_seq_to_int, ssw_align, print_int8_t

# Define a C-struct to store adapter information

ctypedef struct s_query:
    int32_t id
    int32_t size
    int8_t* seq_int
    int32_t found
    int32_t min_len
    int32_t min_score
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cdef class AdapterTrimmer:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~SELF VARIABLES DEFINITION~~~~~~~#
    cdef:
        int min_size, n_query
        int8_t ssw_match, ssw_mismatch, ssw_ambiguous, ssw_gapO, ssw_gapE
        int8_t* score_mat
        s_query* ql
        
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
    
    def __init__(self,
        list adapter_list, # List of adapter sequence to match on all read to be analysed
        bint find_reverse = True, # If true will also search for the reverse complementary sequence of the adapter
        int min_size = 0, # Minimal size of read after trimming to be considered valid
        float min_match_len = 0.8, # Minimal fraction of adapter len that needs to be aligned on the target
        float min_match_score = 1.4, # Minimal score per base for the alignment of adapter and read
        int8_t ssw_match = 2, # Bonus in case of match (POSITIVE)
        int8_t ssw_mismatch = 2, # Malus in case of mismatch (POSITIVE)
        int8_t ssw_ambiguous = 0, # Value in case of ambiguous base (Should remain NULL)
        int8_t ssw_gapO = 3, # Malus in case of gap opening (POSITIVE)
        int8_t ssw_gapE = 1): # Malus in case of gap extension (POSITIVE)

        # Store self value for future usage
        self.min_size = min_size
        self.ssw_match = ssw_match
        self.ssw_mismatch = ssw_mismatch
        self.ssw_ambiguous = ssw_ambiguous
        self.ssw_gapO = ssw_gapO
        self.ssw_gapE = ssw_gapE

        # Init a score matrix
        self.score_mat = score_matrix (ssw_match, ssw_mismatch, ssw_ambiguous)
        
        # Add rev complement of adapters if needed
        if find_reverse:
            adapter_list = adapter_list + [self.rev_comp(adapter) for adapter in adapter_list]
        
        # Init a list of adapters
        self.n_query = len(adapter_list)
        self.ql = <s_query *>malloc(self.n_query * sizeof(s_query))        
        for n, seq in enumerate(adapter_list):
            self.ql[n] = self.build_query (n, seq, min_match_len, min_match_score)
        
    def __repr__(self):
        msg = "ADAPTER TRIMMER CLASS\n"
        msg += "Minimal size : {}\n".format(self.min_size)
        msg += "Number of adater (+rc) : {}\n".format(self.n_query)
        msg += "List of adapter\n"
        for i in range(self.n_query):
            msg += "\tID:{}  Size:{}  Found:{}  Min len:{}  Min score:{}\n".format(
                self.ql[i].id,
                self.ql[i].size,
                self.ql[i].found,
                self.ql[i].min_len,
                self.ql[i].min_score)
            msg += "\t"+str([self.ql[i].seq_int[j] for j in range(self.ql[i].size)])+"\n"
        msg += "SSW parameters\n"
        msg += "Match:{}  Mismatch:{}  Ambiguous:{}  Gap Open:{}  Gap extend:{}\n".format(
            self.ssw_match,
            self.ssw_mismatch,
            self.ssw_ambiguous,
            self.ssw_gapO,
            self.ssw_gapE)
        msg += "Score matrix\n"
        buf = ""
        for i in range(0, 25, 5):
            for j in range(0,5):
                buf += "{:<4}".format(self.score_mat[i+j])
            msg+=buf+"\n"
            buf=""
        return msg
    
    def __str__(self):
        return ("<Instance of AdapterTrimmer Class>\n")
    
    def __dealloc__(self):
        print ("Dealloc memory")
        # Free memory
        
        for i in range(self.n_query):
            free(self.ql[i].seq_int)

        free(self.ql)
        free(self.score_mat)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, ref_seq):
        cdef:
            int32_t ref_size
            int8_t* ref_int
            s_align res
            bint found = False
            s_query query

        # Prepare the reference sequence
        ref_size = len(ref_seq)
        ref_int = DNA_seq_to_int(ref_seq, ref_size)
        
        # Iterate over the adapto query sequence to align against the reference read
        for i in range(self.n_query):
            query = self.ql[i]
            res = ssw_align(
                query = self.ql[i].seq_int,
                queryLen = self.ql[i].size,
                ref = ref_int,
                refLen = ref_size,
                mat = self.score_mat,
                gapO = self.ssw_gapO,
                gapE = self.ssw_gapE)
            print res
            # Update bit field and counters if the score is high enough
            if res.score >= self.ql[i].min_score and res.ref_end-res.ref_begin >= self.ql[i].min_len:
                self.ql[i].found += 1
                found = True
                #update bit field (start, end)
        
        # Dealoc int8_t* allocated for ref_int
        free(ref_int)
        
        # If matches were found trim the sequence and verify its length after trimming
#        if found:
#            start, end = extract_slice (bit_field)
#            if end - start >= self.min_size:
#                return ref_seq[start, end]
#            else:
#                return None
        
#        else:
#            return ref_seq
        
        # Free memory


    #~~~~~~~PRIVATE METHODS~~~~~~~#
    
    cdef s_query build_query (self, int32_t n, char* seq, float min_match_len, float min_match_score):

        # Create a query struct and fill the fields
        cdef s_query q
        
        q.id = n 
        q.size = len(seq)
        q.seq_int = DNA_seq_to_int(seq, q.size)
        q.found = 0
        q.min_len = <int32_t>(min_match_len*q.size) # compute min len and cast in int32_t
        q.min_score = <int32_t>(min_match_score*q.size) # compute min score and cast in int32_t

        return q

    cdef char* rev_comp(self, char* dna):
        """
        One liner returning the reverse complement of a DNA sequence
        """
        return "".join([{'A':'T','T':'A','G':'C','C':'G'}.get(base, "N") for base in dna[::-1]])
