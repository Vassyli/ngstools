"""
Provides a reader for .sam stored Sequencing AlignMents.

.sam files use 1-based coordinates.
"""
import os
from .BaseAlignedRead import BaseAlignedRead
from .BaseReader import BaseReader


class SamReader(BaseReader):
    """Represents the information stored in a SAM file.

    Examples
    --------

    >>> sam = SamReader.open("file.sam")
    >>> for read in sam:
    >>>     print("A read in the sam file ", read.sequence)
    >>> for read in sam.aligned:
    """
    def __iter__(self):
        """ Yields all reads, whether aligned or not. """
        with open(self._filename, "r") as fh:
            for line in fh:
                if line.startswith("@"):
                    # @ are additional information
                    # ToDo: Add this information later!
                    continue

                read = SamAlignedRead.from_line(line)
                yield read

    @property
    def aligned(self):
        """ Yields only aligned reads. """
        for read in self:
            if not read.isAligned:
                continue
            else:
                yield read


class SamAlignedRead(BaseAlignedRead):
    """An aligned read from a SAM file.

    Parameters
    ----------
    qname : str
        query template name, the read identifier.
    flag : int
        the flag column.
    rname : str
        reference name, typically the chromosome the read was assigned to
    pos : int
        position on the reference
    cigar : str
        the CIGAR string
    rnext : str
        reference name of the next read
    pnext : int
        position on the ference of the next read
    tlen : int
        observed Template LENgth
    seq : str
        sequence of the read
    qual : str
        quality description of the read per base (typically in SANGER fastq format), phred+33
    """
    FLAG_MULTIPLE_SEGMENTS = 0x1
    FLAG_PROPERLY_ALIGNED = 0x2
    FLAG_UNMAPPED = 0x4
    FLAG_NEXT_UNMAPPED = 0x8
    FLAG_SEQ_REVERSE_COMPLEMENTED = 0x10
    FLAG_NEXT_SEQ_REVERSE_COMPLEMENTED = 0x20
    FLAG_FIRST = 0x40
    FLAG_LAST = 0x80
    FLAG_SECONDARY_ALIGNMENT = 0x100
    FLAG_NOT_PASSING = 0x200
    FLAG_DUPLICATE = 0x400
    FLAG_SUPPLEMENTARY_ALIGNMENT = 0x800

    @classmethod
    def from_line(cls, line):
        """ Creates a SamAlignedRead object by parsing a SAM line. """
        return cls(*line.split("\t"))

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *kwarg):
        self._qname = None if qname == "*" else qname
        self._flag = int(flag)
        self._rname = None if rname == "*" else rname
        self._pos = int(pos) - 1
        self._mapq = int(mapq)
        self._cigar = cigar
        self._rnext = None if rnext == "*" else (self.rname if rnext == "=" else rnext)
        self._pnext = None if pnext == "0" else int(pnext)
        self._tlen = tlen
        self._qual = qual

        is_rc = True if (self._flag & self.FLAG_SEQ_REVERSE_COMPLEMENTED) > 0 else False

        super().__init__(seq, self._pos, self._rname, is_rc)

    @property
    def queryName(self):
        return self._qname

    @property
    def referenceName(self):
        return self._rname

    @property
    def length(self):
        return len(self._seq)

    @property
    def isAligned(self):
        return False if (self._flag & self.FLAG_UNMAPPED) > 0 else True