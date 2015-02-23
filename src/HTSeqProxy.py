#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class SequenceWithQualitiesProxy (object):
    """
    Recopy HTS_seq in simple object able to go through processes
    Original HTSeq objects are enable to travel through Queues
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~FUNDAMENTAL AND MAGIC METHODS~~~~~~~#

    @classmethod
    def init_from_HTSeq (self, HTSeq):
        return SequenceWithQualitiesProxy(HTSeq.name, HTSeq.seq, HTSeq.qual)

    def __init__ (self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.size = len(self.seq)

    def __len__ (self):
        return len(self.seq)

    def __str__(self):
        if len(self) > 20:
            return "{} : {}...{}  {}...{}".format(self.name, self.seq[:10], self.seq[-10:], self.qual[:10], self.qual[-10:])
        else:
            return "{} : {}  {}".format(self.name, self.seq, self.qual)

    def __getitem__( self, item ):
        return SequenceWithQualitiesProxy(name = self.name, seq = self.seq[ item ], qual = self.qual[ item ])

    @property
    def qualstr( self ):
        qualstr_phred = ""
        for i in self.qual:
            qualstr_phred += str(unichr(i+33))
        return qualstr_phred

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_fastq_str( self):
        return "@{}\n{}\n+\n{}\n".format(self.name, self.seq, self.qualstr)
