"""
FASTQ_Tools.py v1.0.0
    April 29, 2018
    Dennis A. Simpson
    Cleaned up file and set the version to 1.0.0 for paper.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import pathlib
import gzip
import ntpath
import magic

__author__ = 'Dennis A. Simpson'
__version__ = "1.1.0"


class Writer:
    """
    Write new FASTQ file.
    """
    __slots__ = ['log', 'file']

    def __init__(self, log, out_file_string):
        """

        :param log:
        :param out_file_string:
        """
        self.file = open(out_file_string, "w")
        self.log = log

    def lethal_write(self, read):
        outstring = ""

        try:
            assert len(read.seq) == len(read.qual)
        except AssertionError:
            self.log.error("Sequence and quality scores of different lengths! Read Name {0}; Seq Length {1}; Qual "
                           "Length {2}".format(read.name, len(read.seq), len(read.qual)))
            raise SystemExit(1)
        outstring += "@{0}\n{1}\n{2}\n{3}\n".format(read.name, read.seq, read.index, read.qual)

        self.file.write(outstring)
        return True

    def write(self, read_list):
        """
        Write a block of text to new FASTQ file
        :param read_list:
        :return:
        """
        outstring = ""
        for read in read_list:
            try:
                assert len(read.seq) == len(read.qual)
            except AssertionError:
                self.log.error("Sequence and quality scores of different lengths! Read Name {0}; Seq Length {1}; Qual "
                               "Length {2}".format(read.name, len(read.seq), len(read.qual)))
                raise SystemExit(1)
            outstring += "@{0}\n{1}\n{2}\n{3}\n".format(read.name, read.seq, read.index, read.qual)

        self.file.write(outstring)
        read_list.clear()

        return True

    def close(self):
        """
        Close the FASTQ file
        :return:
        """
        self.file.close()
        return True


class FASTQ_Reader:
    __slots__ = ['input_file', 'log', 'name', 'seq', 'index', 'qual', 'read_block', 'file_name', 'fq_file']

    def __init__(self, input_file, log=None):
        """
        Splits the FASTQ read list from the FASTQ Iterator into the lines to be manipulated.  Also does a check to make
        sure the sequence length = quality string length.
        :param input_file:
        :return:
        """

        self.name = None
        self.seq = None
        self.index = None
        self.qual = None
        self.input_file = input_file
        self.log = log
        self.read_block = []
        self.file_name = ntpath.basename(input_file)
        self.fq_file = self.__fastq_file()

    def __fastq_file(self):
        """
        Handles opening the FASTQ file
        :return:
        """
        if len(self.input_file) < 3:
            self.log.warning("FASTQ file parameter missing from options file. Correct error and try again.")
            raise SystemExit(1)

        elif not pathlib.Path(self.input_file).is_file():
            self.log.warning("FASTQ file {0} not found.  Correct error and run again.".format(self.input_file))
            raise SystemExit(1)

        try:
            mime_type = magic.from_file(self.input_file, mime=True).decode()
        except AttributeError:
            mime_type = magic.from_file(self.input_file, mime=True)

        if "text" in mime_type:
            fq_file = open(self.input_file, 'rU')
        elif "gzip" in mime_type:
            fq_file = gzip.open(self.input_file, 'rt', encoding='utf-8')
        else:
            self.log.warning("Unsupported file-type for {0}.  Only TEXT or GZIP Allowed.".format(self.input_file))
            raise SystemExit(1)
        return fq_file

    def line_reader(self):
        """
        Part of the generator to read the FASTQ files
        """
        for line in self.fq_file:
            while True:
                yield line

    def seq_read(self):
        """
        generator to get sequence reads from FAST file into the appropriate blocks of 4 lines.
        """
        read_block = []
        count = 0
        eof = False
        try:
            # for i in range(4):
            while count < 4:
                read_block.append(next(FASTQ_Reader.line_reader(self)))
                count += 1
        except StopIteration:
            eof = True

        if len(read_block) == 4 and not eof:

            self.name = read_block[0].strip("\n").strip("@")
            self.seq = read_block[1].strip("\n").strip()
            self.index = read_block[2].strip("\n").strip()
            self.qual = read_block[3].strip("\n").strip()

            if len(self.seq) != len(self.qual):
                self.log.error("Sequence and quality scores of different lengths! \n{0:s}\n{1:s}\n{2:s}\n{3:s}"
                               .format(self.name, self.seq, self.index, self.qual))
                raise ValueError("Sequence and quality scores of different lengths! \n{0:s}\n{1:s}\n{2:s}\n{3:s}"
                                 .format(self.name, self.seq, self.index, self.qual))
            yield self

        # I am using this as my EOF.  Not so sure the code ever reaches this.
        self.name = None


def read_trim(fastq_read, trim5=None, trim3=None):
    """
    Trim any additional sequences
    :param fastq_read:
    :param trim5:
    :param trim3:
    """
    if trim5 and trim3:
        fastq_read.seq = fastq_read.seq[trim5:-trim3]
        fastq_read.qual = fastq_read.qual[trim5:-trim3]
    elif trim5:
        fastq_read.seq = fastq_read.seq[trim5:]
        fastq_read.qual = fastq_read.qual[trim5:]
    elif trim3:
        fastq_read.seq = fastq_read.seq[:-trim3]
        fastq_read.qual = fastq_read.qual[:-trim3]
