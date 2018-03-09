from collections import defaultdict, namedtuple
import re


class MiniAlignment:

    __loctuple = namedtuple(["transcript", "chrom", "start", "end"])

    def __init__(self, line):
        self.__line = line.rstrip().split("\t")
        self.__query = self.__line[0]
        self.__qlen, self.__qstart, self.__qend = [int(_) for _ in self.__line[1:4]]
        self.__strand = self.__line[4]
        self.__target = self.__line[5]
        self.__tlen, self.__tstart, self.__tend = [int(_) for _ in self.__line[6:9]]
        self.__residue_matches = int(self.__line[9])
        self.__al_block_len = int(self.__line[10])
        self.__mapq = int(self.__line[11])
        self.__attrs = defaultdict(dict)
        self.__location = -1

        for tag in self.__line[12:]:
            tag, typ, val = tag.split(":")
            self.__attrs[tag]["type"] = typ
            self.__attrs[tag]["value"] = val

    query = property(lambda self: self.__query)
    qlen = property(lambda self: self.__qlen)
    qstart = property(lambda self: self.__qstart)
    qend = property(lambda self: self.__qend)

    target = property(lambda self: self.__target)
    tstart = property(lambda self: self.__tstart)
    tend = property(lambda self: self.__tend)
    tlen = property(lambda self: self.__tlen)

    strand = property(lambda self: self.__strand)
    matches = property(lambda self: self.__residue_matches)
    aln_length = property(lambda self: self.__al_block_len)
    mapq = property(lambda self: self.__mapq)
    attrs = property(lambda self: self.__attrs.copy())

    score = property(lambda self: float(self.attrs.get("AS", dict()).get("value", -1)))
    @property
    def primary(self):
        return self.attrs.get("tp", dict()).get("value", "S") == "P"

    def __lt__(self, other):

        if not isinstance(other, type(self)):
            raise TypeError("Cannot compare objects of different classes")

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

    def __hash__(self):
        return hash(self.__line)

    # TraesCS1A01G000100LC_chr1A:6085-47819
    __trn_pattern = re.compile("^(.*)_([^:]*):([0-9]*)-([0-9]*)$")
    __reg_pattern = re.compile("^([^:]*):([0-9]*)-([0-9]*)$")

    @property
    def location(self):
        if self.__location == -1:

            s = re.search(self.__trn_pattern, self.query)
            if s:
                self.__location = self.__loctuple(s.groups()[0], s.groups()[1], int(s.groups()[2]), int(s.groups()[3]))
            else:
                r = re.search(self.__reg_pattern, self.query)
                if r:
                    self.__location = self.__loctuple(s.groups()[0], int(s.groups()[1]), int(s.groups()[2]))
                else:
                    self.__location = None
        return self.__location


class Stitcher:

    def __init__(self):

        self.minimaps = set()
        self.__query = None

    @property
    def query(self):

        return self.__query

    @query.setter
    def query(self, minimap: MiniAlignment):
        if self.__query is None:
            self.__query = minimap.query
        elif self.__query != minimap.query:
            raise ValueError("Trying to add an incompatible alignment: {} vs {}".format(
                self.__query, minimap.query
            ))
        else:
            pass

    def __add__(self, other):

        if not isinstance(other, MiniAlignment):
            raise ValueError("Stritchers only accept MiniAlignments as input!")
        self.query = other.query  # This will raise an error as appropriate
        self.minimaps.add(other)  # As this is a set, it will avoid re-adding lines for nothing

    def stitch(self):

        """This method will try to find the best path for each alignment.
        The basic idea is to use an overlap graph."""

        ordered_maps = sorted(self.minimaps)

        






