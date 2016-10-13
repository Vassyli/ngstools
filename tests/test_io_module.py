import unittest

from ngsTools import io


class TestSamReader(unittest.TestCase):
    def test_read_iteration(self):
        reader = io.read("tests/test_data/test.sam")
        counter = 0
        for read in reader:
            counter += 1
            assert type(read) == io.SamAlignedRead

        assert counter == 7

    def test_read_aligned_iteration(self):
        reader = io.read("tests/test_data/test.sam")
        counter = 0
        for read in reader.aligned:
            counter += 1
            assert type(read) == io.SamAlignedRead

        assert counter == 2

    def test_sam_read_parsing(self):
        # unaligned read
        read = io.SamAlignedRead.from_line("D00535:175:C9NNTANXX:6:2213:3882:83320	4	*	0	0	*	*	0	0	CTGCAGTTCCTGTGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCT	BBC0=11=1E111;=1EGD11=11=<CGGG<E/F01:>F@111111<FG1:	XM:i:0")
        assert type(read) == io.SamAlignedRead
        assert read.isAligned == False
        assert read.queryName == "D00535:175:C9NNTANXX:6:2213:3882:83320"
        assert read.referenceName is None

        # aligned read
        read = io.SamAlignedRead.from_line("D00535:175:C9NNTANXX:6:2213:3776:83396	0	accn|JRYM01000050	63	255	51M	*	0	0	AAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAG	BCCCBGGGGGGGGGGGGFDGDGGG/1F>G1FGDGGGGGGGGGBFGGG/><E	XA:i:0	MD:Z:51	NM:i:0")
        assert type(read) == io.SamAlignedRead
        assert read.isAligned == True
        assert read.queryName == "D00535:175:C9NNTANXX:6:2213:3776:83396"
        assert read.referenceName == "accn|JRYM01000050"

class TesGenomReader(unittest.TestCase):
    def test_read_genom_file(self):
        reader = io.read("tests/test_data/genom.fasta")

        assert "chrI" in reader
        assert type(reader) == io.GenomReader

    def test_genom_slices(self):
        reader = io.read("tests/test_data/genom.fasta")

        assert reader["chrI",0:1] == "A"
        assert reader["chrI",5:6] == "G"
        assert reader["chrI",9:10] == "A"
        assert reader["chrI",10:11] == "T"
        assert reader["chrI",0] == "A"
        assert reader["chrI",5] == "G"
        assert reader["chrI",9] == "A"
        assert reader["chrI",10] == "T"
        assert reader["chrI",0:10] == "ATCGTGCGTA"
        assert reader["chrI",10:20] == "TGCGATGTAC"
        assert reader["chrI",20:30] == "TGCGAGGCAT"
        assert reader["chrI",30:35] == "GTAGT"
        assert reader["chrI",0:20] == "ATCGTGCGTATGCGATGTAC"
        assert reader["chrI",30:40] == "GTAGT"
        assert reader["chrI",0:35] == "ATCGTGCGTATGCGATGTACTGCGAGGCATGTAGT"
        assert "".join(reader["chrI",0,1,2,3,4,5,6,7,8,9]) == "ATCGTGCGTA"
        assert reader["chrI",18:27] == "ACTGCGAGG"
        assert reader["chrII",24:27] == "GTA"

    def test_genom_get_whole_chromosome(self):
        reader = io.read("tests/test_data/genom.fasta")

        assert reader["chrI"] == "ATCGTGCGTATGCGATGTACTGCGAGGCATGTAGT"

    def test_sam_reads_in_genom(self):
        goodRead = io.SamAlignedRead.from_line("alpha\t0\tchrI\t1\t255\t*\t*\t0\t0\tATCGT\tGGGGG")
        badRead = io.SamAlignedRead.from_line("alpha\t4\t*\t1\t255\t*\t*\t0\t0\tATCGT\tGGGGG")
        genom = io.read("tests/test_data/genom.fasta")

        assert goodRead in genom
        assert badRead not in genom

    def test_sam_reads_as_slices(self):
        genom = io.read("tests/test_data/genom.fasta")

        normalRead = io.SamAlignedRead.from_line("alpha\t0\tchrI\t1\t255\t*\t*\t0\t0\tATCGT\tGGGGG")
        assert normalRead.isReverseComplemented is False
        assert normalRead.sequence == "ATCGT"
        assert genom[normalRead] == "ATCGT"

        invertRead = io.SamAlignedRead.from_line("alpha\t16\tchrI\t11\t255\t*\t*\t0\t0\tTGCGA\tGGGGG")
        assert invertRead.isReverseComplemented is True
        assert invertRead.sequence == "TCGCA"
        assert genom[invertRead] == "TCGCA"

    def test_sam_reads_as_extended_slices(self):
        genom = io.read("tests/test_data/genom.fasta")

        normalRead = io.SamAlignedRead.from_line("alpha\t0\tchrI\t11\t255\t*\t*\t0\t0\tTGCGA\tGGGGG")
        assert genom[normalRead.chromosome, normalRead.super(1, 1)] == "ATGCGAT"

        invertRead = io.SamAlignedRead.from_line("alpha\t16\tchrI\t11\t255\t*\t*\t0\t0\tTGCGA\tGGGGG")
        assert genom[invertRead.chromosome, invertRead.super(1, 1)] == "ATCGCAT"

class TestGenomFeatureReader(unittest.TestCase):
    def test_opening(self):
        features = io.read("tests/test_data/features.gff")

        count = 0
        for feature in features:
            count += 1

        assert count == 4
        assert "chrI" in features
        assert "chrII" in features
        assert "GEN001A" in features
        assert "GEN002A" in features

    def test_gff_to_genom_slice(self):
        features = io.read("tests/test_data/features.gff")
        genom = io.read("tests/test_data/genom.fasta")

        assert genom[features["chrI"]] == "ATCGTGCGTATGCGATGTACTGCGAGGCATGTAGT"
        assert genom[features["chrII"]] == "GTTACTGATCAGCTAGTGAGAATCGTAGCTAGCT"
        assert genom[features["GEN001A"]] == "ATGCGA"
        assert genom[features["GEN002A"]] == "TAGCTG"

class TestBedtoolsIntersectionReader(unittest.TestCase):
    def test_opening(self):
        intersection = io.BedtoolsIntersectionReader.open("tests/test_data/intersect.tab")

        count = 0
        for read in intersection:
            count += 1

        assert count == 2

    def test_intersection_to_genom_slice(self):
        intersection = io.BedtoolsIntersectionReader.open("tests/test_data/intersect.tab")
        genom = io.read("tests/test_data/genom.fasta")

        items = []
        for item in intersection:
            items.append(item)

        assert genom[items[0]] == "TGCGT"
        assert genom[items[1]] == "CTGATC"

if __name__ == '__main__':
    unittest.main()