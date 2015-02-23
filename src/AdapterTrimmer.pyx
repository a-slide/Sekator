#~~~~~~~CIMPORTS~~~~~~~#

# C standard library import
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, free

# Local package import
from ssw cimport s_align, score_matrix, DNA_seq_to_int, ssw_align

#~~~~~~~STRUCTURES~~~~~~~#

ctypedef struct s_query:
    int32_t id
    int32_t size
    int8_t* seq_int
    int32_t count
    int32_t min_len
    int32_t min_score

#    @typedef struct to store adapters information
#    @field  id          Number of identification
#    @field  size        Size of the reference in base
#    @field  seq_int     Base sequence converted in numeric values (A,a=0; C,c=1; G,g=2; T,t=3; other char = 4)
#    @field  count       Number of time the adapter is found
#    @field  min_len     Minimal length of the adapter to match on the reference
#    @field  min_score   Minimal score of the adapter match on the reference


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cdef class AdapterTrimmer:
    """
    Load a list of adapter sequences and align them with reads through a fast C implemention
    of Smith and Waterman alignment algorithm
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~SELF VARIABLES DEFINITION~~~~~~~#

    cdef:
        uint32_t min_size, n_query, total, untrimmed, trimmed, fail
        uint64_t base_trimmed
        int8_t ssw_match, ssw_mismatch, ssw_ambiguous, ssw_gapO, ssw_gapE
        int8_t* score_mat
        s_query* ql


    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, list adapter_list, int32_t min_size=30,\
        float min_match_len=0.3, float min_match_score=1, int8_t ssw_match=2,\
        int8_t ssw_mismatch=2, int8_t ssw_ambiguous=0, int8_t ssw_gapO=3, int8_t ssw_gapE=1):

#        Initialize AdapterTrimmer from a list of adapter sequence and compute the score matrix
#        based on the provided ssw scores.
#        @param adapter_list    List of adapter sequence to match on all read to be analysed
#        @param min_size        Minimal size of read after trimming to be considered valid
#        @param min_match_len   Minimal fraction of adapter len that needs to be aligned on the target
#        @param min_match_score Minimal score per base for the alignment of adapter and read
#        @param ssw_match       Gain in case of match (POSITIVE)
#        @param ssw_mismatch    Penalty in case of mismatch (POSITIVE)
#        @param ssw_ambiguous   Value in case of ambiguous base (Should remain NULL)
#        @param ssw_gapO        Penalty in case of gap opening (POSITIVE)
#        @param ssw_gapE        Penalty in case of gap extension (POSITIVE)
#        @note Default values determined for 100pb reads with randomly generated 60 pb adaptors

        # Store self value for future usage
        self.min_size = min_size
        self.ssw_match = ssw_match
        self.ssw_mismatch = ssw_mismatch
        self.ssw_ambiguous = ssw_ambiguous
        self.ssw_gapO = ssw_gapO
        self.ssw_gapE = ssw_gapE

        # Init Counters
        self.total = 0
        self.untrimmed = 0
        self.trimmed = 0
        self.fail = 0
        self.base_trimmed = 0

        # Init a score matrix
        self.score_mat = score_matrix (ssw_match, ssw_mismatch, ssw_ambiguous)

        # Init a list of adapters
        self.n_query = len(adapter_list)
        self.ql = <s_query *>malloc(self.n_query * sizeof(s_query))
        for n, seq in enumerate(adapter_list):
            self.ql[n] = self.build_query (n, seq, min_match_len, min_match_score)

    def __repr__(self):
        msg = "ADAPTER TRIMMER CLASS\n"
        msg += "Minimal size:{} Total:{} Untrimmed:{} Trimmed:{} Fail:{} Base Trimmed:{}\n".format(
            self.min_size, self.total, self.untrimmed, self.trimmed, self.fail, self.base_trimmed)
        msg += "Number of adater (+rc) : {}\n".format(self.n_query)
        msg += "List of adapter\n"
        for i in range(self.n_query):
            msg += "\tID:{} Size:{} Found:{} Min len:{} Min score:{}\n".format(
                self.ql[i].id, self.ql[i].size, self.ql[i].count, self.ql[i].min_len, self.ql[i].min_score)
            msg += "\tInteger sequence : {}\n".format("".join([str(self.ql[i].seq_int[j]) for j in range(self.ql[i].size)]))
        msg += "SSW parameters\n"
        msg += "Match:{} Mismatch:{} Ambiguous:{} Gap Open:{} Gap extend:{}\n".format(
            self.ssw_match, self.ssw_mismatch, self.ssw_ambiguous, self.ssw_gapO, self.ssw_gapE)
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

    def __call__(self, object seq):

#        Find imperfect adapter matches with ssw algorithm and extract the larger interval of
#        reference sequence that do not overlap interval match
#        @param seq HTSeq.SequenceWithQualities object

        cdef:
            int32_t seq_size, i, start_max=0, end_max=0, inter_max=0, start=0, inter=0
            int8_t* seq_int
            int8_t* bool_mat
            s_align res
            int8_t found = 0
            s_query query

        self.total += 1

        # Prepare the reference sequence by converting the HTSfastq sequence in int8_t*
        seq_size = len(seq)
        seq_int = DNA_seq_to_int(seq.seq, seq_size)

        # Init a zero padded matrix of short 8 bits int
        bool_mat = <int8_t *>malloc(seq_size * sizeof(int8_t))
        for i in range(seq_size):
            bool_mat[i] = 0

        # Iterate over the adapter query sequence to align against the reference read
        for i in range(self.n_query):
            res = ssw_align(
                query = self.ql[i].seq_int,
                queryLen = self.ql[i].size,
                ref = seq_int,
                refLen = seq_size,
                mat = self.score_mat,
                gapO = self.ssw_gapO,
                gapE = self.ssw_gapE)

            # Update bool mat and counters if the score is high enough
            if res.score >= self.ql[i].min_score and res.ref_end-res.ref_begin >= self.ql[i].min_len:
                self.ql[i].count += 1
                found = 1
                for i in range (res.ref_begin, res.ref_end+1):
                    bool_mat[i] = 1

        # Dealoc int8_t* allocated for seq_int
        free(seq_int) ##

        # Print bool mat with number to verify the match
#        print ("".join(["{:<3}".format(i) for i in range (seq_size)]))
#        bool_str =""
#        for i in range(seq_size):
#            bool_str += "{}".format("XXX" if bool_mat[i] else "---")
#        print (bool_str)

        # Return an unchanged sequence if no significant match were found
        if not found:
            self.untrimmed += 1
            free(bool_mat) ##
            return seq

        # Find the longer interval without primers
        for i in range(seq_size):
            if bool_mat[i]:
                start = i+1
                inter = 0
            else:
                inter += 1
                if inter > inter_max:
                    inter_max = inter
                    start_max = start
                    end_max = i
        #print ("start {}  end {}  inter {}".format(start_max, end_max, inter_max))

        free(bool_mat) ##
        self.base_trimmed += seq_size-inter_max

        # Return None if the size of the longer interval
        if inter_max < self.min_size:
            self.fail += 1
            return None

        # Finally in the last case
        self.trimmed += 1
        return seq[start_max:end_max+1]

    def get_summary (self):

        cdef:
            dict summary
            int32_t i

        summary = {}
        summary["total"] = int(self.total)
        summary["untrimmed"] = int(self.untrimmed)
        summary["trimmed"] = int(self.trimmed)
        summary["fail"] = int(self.fail)
        summary["base_trimmed"] = int(self.base_trimmed)
        summary["adapter_found"] = []

        for i in range(self.n_query):
            summary["adapter_found"].append(int(self.ql[i].count))

        return summary


    #~~~~~~~PRIVATE METHODS~~~~~~~#

    cdef s_query build_query (self, int32_t n, char* seq, float min_match_len, float min_match_score):
#       Create a query struct and fill the fields
        cdef s_query q

        q.id = n
        q.size = len(seq)
        q.seq_int = DNA_seq_to_int(seq, q.size)
        q.count = 0
        q.min_len = <int32_t>(min_match_len*q.size) # compute min len and cast in int32_t
        q.min_score = <int32_t>(min_match_score*q.size) # compute min score and cast in int32_t

        return q
