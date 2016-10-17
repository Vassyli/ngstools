"""Subpackage for handling import and export of sequencing data.

"""

import os

from .GenomFeatureReader import GenomFeatureReader
from .GenomReader import GenomReader
from .SamReader import SamReader, SamAlignedRead
from .BaseAlignedRead import BaseAlignedRead
from .BedtoolsIntersectionReader import BedtoolsIntersectionReader, BedtoolsIntersectionItem

_extension_to_reader = {
    "sam": SamReader.open,
    "fasta": GenomReader.open,
    "fa": GenomReader.open,
    "fna": GenomReader.open,
    "gff": GenomFeatureReader.open
}


def read(filename):
    """Tries to read a given filename.

    This function guesses the filetype based in the filename and returns
    the appropriate object or raises an exception if the filetype is unknown.

    Parameters
    ----------
    filename : str
        The filename of the file to open.

    Returns
    -------
    Reader
        A reader representing the file and offering methods to work with the file.

    Supported file formats
    ----------------------
    A list of files supported by this function.

        .sam      Sequence Alignment/Map     ngsTools.io.SamReader
        .fa
        .fasta    FASTA file format          ngsTools.io.GenomReader
    """
    basename = os.path.basename(filename)
    extension = basename.split(".")[-1]

    if extension not in _extension_to_reader:
        raise TypeError("ngsTools.io does not support this file extension «" + extension + "»")

    return _extension_to_reader[extension](filename)
