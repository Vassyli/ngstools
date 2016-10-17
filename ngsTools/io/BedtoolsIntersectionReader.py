from .BaseAlignedRead import BaseAlignedRead
from .BaseReader import BaseReader

class BedtoolsIntersectionReader(BaseReader):
    """A reader for tabulated intersection files from bedtools
    """
    def __iter__(self):
        with open(self._filename) as fh:
            for line in fh:
                yield BedtoolsIntersectionItem.from_line(line)


class BedtoolsIntersectionItem(BaseAlignedRead):
    @classmethod
    def from_line(cls, line):
        return cls(*line.split("\t"))

    def __init__(self, chromosome, start, end, qname, strand, type, fstart, fstop, fid):
        self._chromosome = chromosome
        self._start = int(start)
        self._stop = int(end)
        self._qname = qname
        self._strand = strand
        self._type = type
        self._fstart = int(fstart) - 1 # fstart and fstop comes from gff, which is 1-based..
        self._fstop = int(fstop) - 1
        self._fid = fid

        super().__init__("N" * (self._stop - self._start), self._start, self._chromosome,
                         True if self._strand == "-" else False)

    @property
    def chromosome(self):
        return self._chromosome

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    def __len__(self):
        return self.stop - self.start

    @property
    def queryName(self):
        return self._qname

    @property
    def feature_id(self):
        return self._fid

