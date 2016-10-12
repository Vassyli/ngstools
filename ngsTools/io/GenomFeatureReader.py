
from .BaseAlignedRead import BaseAlignedRead
from .BaseReader import BaseReader


class GenomFeatureReader(BaseReader):
    """A reader for genom feature files.

    """
    def __iter__(self):
        with open(self._filename) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                yield GenomFeatureItem.from_line(line)

    def __contains__(self, item):
        self.prepare()

        if item in self._feature_dict:
            return True
        else:
            return False

    def __getitem__(self, item):
        self.prepare()

        if item in self:
            return self._feature_dict[item]
        else:
            raise IndexError("{0} feature not found in .gff file".format(item))

    def prepare(self):
        if "_prepared" in self.__dict__:
            return None

        self._feature_dict = {}

        for feature in self:
            self._feature_dict[feature.get_attribute("ID")] = feature

        self._prepared = True


class GenomFeatureItem(BaseAlignedRead):
    @classmethod
    def from_line(cls, line):
        return cls(*line.split("\t"))

    @staticmethod
    def parse_attributes(attributes):
        ret = {}
        for item in attributes.split(";"):
            keyval = item.split("=")
            ret[keyval[0].upper()] = keyval[1]
        return ret

    def __init__(self, chromosome, source, feature, start, end, score, strand, frame, attributes):
        self._chromosome = chromosome
        self._source = source
        self._feature = feature
        self._start = int(start) - 1
        self._stop = int(end)
        self._score = score
        self._strand = strand
        self._frame = int(frame) if frame != "." else None
        self._attributes = __class__.parse_attributes(attributes)

        super().__init__("N" * (self._stop - self._start), self._start, self._chromosome, True if self._strand == "-" else False)

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
    def step(self):
        return -1 if self._strand == "-" else None

    def get_attribute(self, attribute):
        if attribute in self._attributes:
            return self._attributes[attribute]
        else:
            return None

