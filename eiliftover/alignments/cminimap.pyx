import re
from libc.stdlib cimport malloc, free

cdef class Location:

    cdef public str chrom
    cdef public int start
    cdef public int end
    cdef public str transcript

    def __cinit__(self, str chrom, int start, int end, str transcript):
        self.chrom, self.transcript = chrom, transcript
        self.start, self.end = start, end


cdef class MiniAlignment:

    """Class to parse very quickly a single MiniMap alignment"""

    cdef readonly str query, target
    cdef readonly int qlen, qstart, qend
    cdef readonly int tlen, tstart, tend
    cdef readonly str strand
    cdef int residue_matches
    cdef readonly int aln_length, mapq
    cdef str line
    cdef readonly dict attrs
    cdef readonly int score
    cdef readonly Location location
    cdef readonly double as

    __trn_pattern = re.compile("^(.*)_([^:]*):([0-9]*)-([0-9]*)$")
    __reg_pattern = re.compile("^([^:]*):([0-9]*)-([0-9]*)$")

    def __cinit__(self, str line):

        cdef list __fields
        cdef int i
        cdef str tag, typ, val

        self.line = line
        __fields = line.rstrip().split("\t")
        self.query = __fields[0]
        self.location = self.define_location(self.query)
        self.qlen, self.qstart, self.qend = [int(_) for _ in __fields[1:4]]
        self.strand = __fields[4]
        self.target = __fields[5]
        self.tlen, self.tstart, self.tend = [int(_) for _ in __fields[6:9]]
        self.residue_matches = int(__fields[9])
        self.aln_length = int(__fields[10])
        self.mapq = int(__fields[11])
        self.attrs = dict()
        self.as = float("-inf")

        for i in range(12, len(__fields)):
            tag = __fields[i]
            tag, typ, val = tag.split(":")
            self.attrs[tag] = dict()
            self.attrs[tag]["type"] = typ
            self.attrs[tag]["value"] = val
            if tag == "AS":
                self.as = float(val)


    cdef Location define_location(self, str query):
        cdef str chrom
        cdef str transcript
        cdef int start, end
        cdef tuple groups

        s = re.search(self.__trn_pattern, self.query)

        if s:
            groups = s.groups()
            chrom = groups[1][:]
            start, end = int(groups[2][:]), int(groups[3][:])
            transcript = groups[0][:]

            return Location(chrom, start, end, transcript)
        else:
            r = re.search(self.__reg_pattern, self.query)
            if r:
                groups = r.groups()
                chrom = groups[0][:]
                start = int(groups[1][:])
                end = int(groups[2][:])
                transcript = None
                return Location(chrom, start, end, transcript)
            else:
                return None

    def __hash__(self):
        return hash(self.line)

    def __lt__(self, MiniAlignment other):

        # if not isinstance(other, type(self)):
        #    raise TypeError("Cannot compare objects of different classes")

        if self.query != other.query:
            return self.query < other.query
        elif self.qstart != other.qstart:
            return self.qstart < other.qstart
        elif self.qend != other.qend:
            return self.qend < other.qend
        elif self.strand != other.strand:
            return self.strand == "-"
        else:
            return False

    cpdef bint is_complete(self, int flank):

        if self.qstart <= flank and self.qend >= self.qlen - flank:
            return True
        else:
            return False


cdef class Stitcher:

    cdef public set minimaps
    cdef list sorted_minis
    cdef int flank
    cdef str __query
    # cdef str query

    def __cinit__(self, int flank):

        self.minimaps = set()
        self.sorted_minis = []
        self.__query = None
        self.flank = flank

    cpdef query_getter(self):

        return self.__query

    cpdef query_setter(self, MiniAlignment minimap):
        if self.__query is None:
            self.__query = minimap.query
        elif self.__query != minimap.query:
            raise ValueError("Trying to add an incompatible alignment: {} vs {}".format(
                self.__query, minimap.query
            ))
        else:
            pass

    query = property(query_getter, query_setter)

    def add(self, MiniAlignment other):

        self.query = other  # This will raise an error as appropriate
        self.minimaps.add(other)  # As this is a set, it will avoid re-adding lines for nothing

    def sort(self):

        self.sorted_minis = sorted(list(self.minimaps))

    def stitch(self):

        """This method will try to find the best path for each alignment.
        The basic idea is to use an overlap graph."""

        # ordered_maps = sorted(self.minimaps)
        raise NotImplementedError()

    def sort_by_aln(self):

        return sorted(list(self.minimaps), key=lambda mini: mini.as, reverse=True)

    @property
    def is_best_single_good_alignment(self):

        """This boolean property will return True if there is at least one alignment which is """

        if self.minimaps:
            best = self.sort_by_aln()[0]
            return best.is_complete(self.flank)
        else:
            return False

    @property
    def best_complete_alignment(self):

        cdef int pos
        cdef list sorted
        cdef MiniAlignment minimap
        sorted = self.sort_by_aln()

        for pos in range(len(self.minimaps)):
            minimap = sorted[pos]
            if minimap.is_complete(self.flank):
                return minimap

        return None
