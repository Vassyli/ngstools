from ..utils import get_reverse_complement


class BaseAlignedRead:
    """Represents a read sequence mapped to a genom.

    Parameters
    ----------
    seq : str
        The sequence of the read
    start : int
        The position of the first base in seq on the alignment map. 0-based index.
    chromosome : str
        The chromosome identifier name of the alignment
    is_rc : bool
        True if the actual read is the reverse complement of seq, not seq itself.
    """
    def __init__(self, seq, start, chromosome, is_rc = False):
        self._start = start
        self._seq = seq.upper()
        self._stop = start + len(seq)
        self._chromosome = chromosome
        self._is_rc = is_rc

    @property
    def start(self):
        """Position of the first nucleotide on the genom map."""
        return self._start

    @property
    def stop(self):
        """Position of the last nucleotide on the genom map + 1

        Since stop - start gives the length, the position stop position is actually the position
        after the last base. This makes it compatible with python counting as well, since [0:10]
        does only include the numbers from 0-9, excluding the last."""
        return self._stop

    @property
    def step(self):
        """ -1 if reverse Complemented should be used instead. """
        return -1 if self.isReverseComplemented is True else None

    @property
    def isReverseComplemented(self):
        """ Is True if the actual sequence is the reverse complement."""
        return self._is_rc

    @property
    def chromosome(self):
        """ The chromosome identifier name."""
        return self._chromosome

    @property
    def sequence(self):
        """ The actual sequence of the read.

        This property represents the actual sequence of the read, not necessarily the sequence
        in the original file or read. """
        return get_reverse_complement(self._seq) if self.isReverseComplemented is True else self._seq

    def super(self, left = 0, right = 0):
        """ Creates a read instance which is longer than the current one.

        Creates a read which the current read is a subset of. Bases left and right of the current read sequence
        are filled with N.

        Parameters
        ----------
        left : int
            Number of bases which get added on the left end.
        right: int
            Number of bases which get added on the right end.

        Returns
        -------
        BaseAlignedRead
            A new BaseAlignedRead with the extended range.
        """
        if left < 0 or right < 0:
            raise TypeError("right and left need to be >= 0")

        new_start = self.start - left

        if new_start < 0:
            raise Exception("Cannot read at starts < 0!")

        new_sequence = "N" * left + self._seq + "N" * right

        if self.isReverseComplemented:
            return __class__(new_sequence, new_start, self._chromosome, True)
        else:
            return __class__(new_sequence, new_start, self._chromosome, False)

    def supra(self, left = 0, right = 0):
        """ Creates a read instance which is shorter than the current one.

        Creates a read which the current read is a superset of. Bases left and right of the current read sequence
        are removed.

        Parameters
        ----------
        left : int
            Number of bases which get removed on the left end.
        right: int
            Number of bases which get removed on the right end.

        Returns
        -------
        BaseAlignedRead
            A new BaseAlignedRead with the shortened range.
        """
        if left < 0 or right < 0:
            raise TypeError("right and left need to be >= 0")

        new_start = self.start + left
        new_stop = self.start - right

        if new_stop <= new_start:
            raise TypeError("The new right end is at the position or smaller than the new small end.")

        new_sequence = self._seq[left:-right]

        if self.isReverseComplemented:
            return __class__(new_sequence, new_start, self._chromosome, True)
        else:
            return __class__(new_sequence, new_start, self._chromosome, False)
