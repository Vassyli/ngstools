"""
Provides a reader genomic .fasta files
"""
import os
from .BaseReader import BaseReader
from .BaseAlignedRead import BaseAlignedRead
from .SamReader import SamAlignedRead
from ..utils import get_reverse_complement


class GenomReader(BaseReader):
    def __init__(self, filename):
        super().__init__(filename)
        self.prepare()

    def prepare(self):
        self._chromosomes = {}
        self._listOfChromosomes = []

        lastChromosome = "(empty)"
        oldPos = -1

        with open(self._filename, "r") as fh:
            while True:
                pos = fh.tell()
                line = fh.readline().strip()
                if pos == oldPos:
                    break
                else:
                    oldPos = pos

                if line.startswith(">"):
                    lastChromosome = line[1:].split(" ")[0]
                    self._chromosomes[lastChromosome] = []
                    continue

                if len(line) > 0:
                    self._chromosomes[lastChromosome].append((pos, len(line)))

        self._listOfChromosomes = self._chromosomes.keys()

    def __contains__(self, item):
        if type(item) == SamAlignedRead:
            return True if item.isAligned and item.referenceName in self._listOfChromosomes else False
        else:
            return True if item in self._listOfChromosomes else False

    def __getitem__(self, index):
        if type(index) == tuple:
            if len(index) < 2:
                raise TypeError("Invalid argument: Index must not be a one-sized tuple.")

            if type(index[0]) != str:
                raise TypeError("First index must be a string which references a chromosome")

            chromosome = self._chromosomes[index[0]]
            listOfSlices = index[1:]
        else:
            if type(index) == str:
                # only a chromosome string, do something else.
                if index in self:
                    length = 0
                    for row in self._chromosomes[index]:
                        length += row[1]
                    return self[index,0:length]
                else:
                    raise TypeError("Chromosome is not in this genom.")
            elif not isinstance(index, BaseAlignedRead):
                raise TypeError("A ShortRead can only be used as the sole index")

            chromosome = self._chromosomes[index.chromosome]
            listOfSlices = [index]

        ret = []
        fh = open(self._filename)

        for slicePiece in listOfSlices:
            if type(slicePiece) == int:
                slicePiece = slice(slicePiece, slicePiece + 1)

            length = slicePiece.stop - slicePiece.start
            lengthLeft = length
            chromosome_pos = - 1
            sequence = ""
            startFound = False
            stopFound = False

            for i in range(0, len(chromosome)):
                seek, linelength = chromosome[i]
                # add linelength to chromosome position
                if chromosome_pos < 0:
                    chromosome_pos = 0
                else:
                    chromosome_pos += chromosome[i-1][1]

                # Search for start position first
                if startFound == False:
                    # skip if position is further away than current line
                    if slicePiece.start >= (chromosome_pos + linelength):
                        continue

                    # now we should be on the correct line
                    pos = slicePiece.start - chromosome_pos
                    fh.seek(seek + pos)
                    startFound = True

                    if slicePiece.stop > (chromosome_pos + linelength):
                        sequence += fh.readline().strip()
                        lengthLeft -= len(sequence)
                    else:
                        sequence += fh.readline().strip()[0:length]
                        lengthLeft -= len(sequence)
                        stopFound = True
                        break

                if startFound == True and stopFound == False and lengthLeft > 0:
                    if lengthLeft - linelength <= 0:
                        # stop found
                        line = fh.readline().strip()[0:lengthLeft]
                        lengthLeft -= len(line)
                        sequence += line
                        stopFound = True
                        break
                    else:
                        lengthLeft -= linelength
                        sequence += fh.readline().strip()

            if startFound and stopFound:
                if slicePiece.step == -1:
                    ret.append(get_reverse_complement(sequence).upper())
                else:
                    ret.append(sequence.upper())
            elif startFound and not stopFound:
                raise Exception("Out of Boundaries exception: ", slicePiece)

        fh.close()

        if len(ret) == 1:
            return ret[0]
        else:
            return ret