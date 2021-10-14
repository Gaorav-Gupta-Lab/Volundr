#!python
#cython: language_level=3

import gzip
import ntpath
import magic

cdef class Reader:
    """
    Main class that creates FASTQ reads
    """
    cdef public input_file, eof, name, seq, index, qual, read_block, file_name, fq_file, batch_size

    def __init__(self, input_file, batch_size):
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
        self.eof = False
        self.input_file = input_file
        self.batch_size = batch_size
        self.file_name = ntpath.basename(input_file)
        self.fq_file = self.__fastq_file()

    def __fastq_file(self):
        """
        Create FASTQ file object.
        :return:
        """

        try:
            mime_type = magic.from_file(self.input_file, mime=True).decode()
        except AttributeError:
            mime_type = magic.from_file(self.input_file, mime=True)

        if "text" in mime_type:
            fq_file = open(self.input_file, 'rU')
        elif "gzip" in mime_type:
            fq_file = gzip.open(self.input_file, 'rt', encoding='utf-8')

        return fq_file

    def line_reader(self):
        """
        Part of generator for FASTQ reads
        """
        for line in self.fq_file:
            while True:
                yield line

    def seq_read(self):
        """
        Generator reads FASTQ file creating read block.
        """
        read_block = []
        count = 0
        eof = False
        batch_block = []
        batch_counter = 0
        try:
            while count < 4:
                read_block.append(next(Reader.line_reader(self)))
                count += 1
        except StopIteration:
            self.eof = True

        if len(read_block) == 4 and not eof:
            count = 0
            self.name = read_block[0].strip("\n").strip("@")
            self.seq = read_block[1].strip("\n").strip()
            self.index = read_block[2].strip("\n").strip()
            self.qual = read_block[3].strip("\n").strip()

            if len(self.seq) != len(self.qual):
                self.log.error("Sequence and quality scores of different lengths! \n{0:s}\n{1:s}\n{2:s}\n{3:s}"
                               .format(self.name, self.seq, self.index, self.qual))
                raise ValueError("Sequence and quality scores of different lengths! \n{0:s}\n{1:s}\n{2:s}\n{3:s}"
                                 .format(self.name, self.seq, self.index, self.qual))

    def grouper(self):
        read_batch = []
        if self.eof:
            return read_batch
        for i in range(self.batch_size):
            self.seq_read()

            read_batch.append((self.name, self.seq, self.index, self.qual))
        yield read_batch

