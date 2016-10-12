"""Utility package for working with sequences

"""

_nucleotide_complement_map = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "S": "S",
    "W": "W",
    "R": "Y",
    "Y": "R",
    "K": "M",
    "M": "K",
    "D": "H",
    "H": "D",
    "V": "B",
    "B": "V"
}


def get_reverse_complement(sequence):
    """Returns the reverse complement of a given sequence.

    This function takes a sequence and translates the iupac nucleotide
    1-base-code to it's reverse complement. It not only regnognizes A, C, G and T,
    but also superset of these (like N, S (for G|C), etc).

    Parameters
    ----------
    sequence : string
        The sequence which get complemented and reversed.

    Returns
    -------
    string
        The complemented and reversed sequence.

    Examples
    --------

    >>> get_reverse_complement("ATCGATCG")
    'CGATCGAT'
    >>> get_reverse_complement("AANTDHW")
    'WDHANTT'
    """
    revCompl = ""

    for nucleotide in sequence[::-1]:
        revCompl += _nucleotide_complement_map[nucleotide]

    return revCompl
