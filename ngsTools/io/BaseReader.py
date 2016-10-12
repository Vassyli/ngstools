import os


class BaseReader:
    @classmethod
    def open(cls, filename):
        """ Returns a representation of the information in the given file.

        Parameters
        ----------
        filename : str
            A file to represent

        Returns
        -------
        self.__class__
            A reader object representing the information contained in the file.
        """
        if not os.path.exists(filename):
            raise FileNotFoundError("File «{0}» not found".format(filename))

        reader = cls(filename)
        return reader

    def __init__(self, filename):
        self._filename = filename