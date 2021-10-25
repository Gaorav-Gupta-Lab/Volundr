"""
Some helpers to make it easier to manipulate FASTQ files.
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27517
@copyright: 2020
"""

import pathlib
import gzip
import ntpath
import itertools
from operator import add
import collections
import time
import natsort
import pathos
import Levenshtein
import magic
from Valkyries import Tool_Box, FileWriter

__author__ = 'Dennis A. Simpson'
__version__ = "0.17.7"


class FastqSplitter:
    """
    This class will chop up a FASTQ file or a pair of FASTQ files for use in the aligner or quality analysis.
    The class name here no longer reflects the primary function of this class.  The file split function is not needed
    when using the parallel version of the aligners.
    """

    def __init__(self, args, log, fastq1, fastq2, index1=None, index2=None, paired_end=False):
        self.log = log
        self.args = args
        self.fastq1_file = fastq1
        self.fastq2_file = fastq2
        self.index1_file = index1
        self.index2_file = index2
        self.paired_end = paired_end
        self.read_count = None

    def new_file_size(self, line_count):
        """

        :param line_count:
        :return:
        """
        # Adjust the number of lines to be written in each file.
        self.log.info("Begin computing file size of new FASTQ files.")
        self.read_count = line_count/4
        file_size = self.read_count/int(self.args.Split)

        if not file_size % 1 == 0:
            file_size = 1 + int(file_size)

        self.log.info("FASTQ file size computation done.  {0} files containing {1} reads will be written, last file "
                      "may differ.".format(self.args.Split, file_size))

        return file_size

    def file_line_counter(self):
        """
        This function counts the number of lines in FASTQ1.  This information is needed to accurately split the FASTQ
        file(s) for parallel alignment.
        :return:
        """

        infile = FASTQ_Reader(self.args.FASTQ1, self.log)
        start_time = time.time()
        self.log.info("Begin counting lines in {}.".format(self.args.FASTQ1))

        bufgen = \
            itertools.takewhile(lambda x: x, (infile.fq_file.read(1024 * 1024) for _ in itertools.repeat(None)))

        line_count = sum(buf.count('\n') for buf in bufgen)

        self.log.info("Found {0} reads in {1}. Count took {2} seconds"
                      .format(line_count/4, infile.file_name, int(time.time() - start_time)))

        infile.fq_file.close()
        return line_count

    def temp_file_writer(self, limit):
        """
        Write the temporary FASTQ files.  Also create list of temporary BAM file names for use later.
        :return:
        """

        self.log.info("Begin writing temporary FASTQ files.")
        i = 0
        temp_file1 = None
        temp_file2 = None
        fastq_file_list = []
        bam_file_list = []
        read_count = 0
        limit_counter = 0

        while read_count <= self.read_count:
            try:
                # This generator is returning actual reads not lines.
                fastq1_read = next(self.fastq1_file.seq_read())
                fastq2_read = next(self.fastq2_file.seq_read())
                if self.index1_file is not None:
                    fastq3_read = next(self.index1_file.seq_read())
            except StopIteration:
                read_count += 1
                continue

            read_count += 1

            try:
                fastq1_n_frac = fastq1_read.seq.count("N")/len(fastq1_read.seq)
                fastq2_n_frac = fastq2_read.seq.count("N")/len(fastq2_read.seq)
            except ZeroDivisionError:
                continue

            # Apply Filters
            if (len(fastq1_read.seq) < int(self.args.Minimum_Length) or
                    len(fastq2_read.seq) < int(self.args.Minimum_Length) or
                    fastq1_n_frac >= float(self.args.N_Limit) or
                    fastq2_n_frac >= float(self.args.N_Limit)):
                continue

            if limit_counter % limit == 0:
                if temp_file1:
                    temp_file1.close()
                    limit_counter = 0
                if temp_file2:
                    temp_file2.close()

                file1 = "{0}{1}_R1_tmp_{2}.fastq.gz".format(self.args.WorkingFolder, self.args.Job_Name, i)
                file2 = "{0}{1}_R2_tmp_{2}.fastq.gz".format(self.args.WorkingFolder, self.args.Job_Name, i)
                bam_file_list.append("{0}{1}_R1_tmp_{2}.bam".format(self.args.WorkingFolder, self.args.Job_Name, i))
                fastq_file_list.append((file1, file2))
                temp_file1 = Writer(self.log, file1)
                temp_file2 = Writer(self.log, file2)

                self.log.info("Writing {0} and {1}".format(file1, file2))
                i += 1

            limit_counter += 1

            # BAM files are missing the barcodes because of a space in some of the header files.
            # fastq1_read.name = fastq1_read.name.replace(" ", ":")
            # fastq2_read.name = fastq2_read.name.replace(" ", ":")

            # Add the UMT's to the header.
            if self.args.HaloPLEX:
                umi = fastq3_read.seq
                header1 = "{0}|{1}:{2}".format(fastq1_read.name.split(":")[-1], umi, fastq1_read.name)
                header2 = "{0}|{1}:{2}".format(fastq2_read.name.split(":")[-1], umi, fastq2_read.name)

            elif self.args.ThruPLEX:
                # header1 = "{0}|{1}".format(fastq1_read.name.split(":")[-1], fastq1_read.name)
                umt1 = fastq1_read.seq[:6]
                umt2 = fastq2_read.seq[:6]
                header1 = "{0}|{1}:{2}".format(umt1, umt2, fastq1_read.name)
                header2 = "{0}|{1}:{2}".format(umt1, umt2, fastq2_read.name)
            else:
                Tool_Box.debug_messenger("Only HaloPLEX or ThruPLEX currently enabled.")
                self.log.error("Only HaloPLEX or ThruPLEX currently enabled.")
                raise SystemExit(1)

            # Trim adapter sequences from 5' end if needed.
            if int(self.args.trim) > 0:
                fastq1_read.seq = fastq1_read.seq[int(self.args.trim):]
                fastq1_read.qual = fastq1_read.qual[int(self.args.trim):]
                fastq2_read.seq = fastq2_read.seq[int(self.args.trim):]
                fastq2_read.qual = fastq2_read.qual[int(self.args.trim):]

            fastq1_read.name = header1
            fastq2_read.name = header2

            temp_file1.write(self.fastq1_file)
            temp_file2.write(self.fastq2_file)

        if temp_file1:
            temp_file1.close()
        if temp_file2:
            temp_file2.close()

        self.log.info("All temporary FASTQ files written")

        return fastq_file_list, bam_file_list

    def file_writer(self):
        """
        Process FASTQ file(s) and write new version(s) suitable for aligners.  Return file name.
        :return:
        """

        self.log.info("Begin writing temporary FASTQ files.")
        current_read_count = 0
        file1 = "{0}{1}_R1_processed.fastq.gz".format(self.args.WorkingFolder, self.args.Job_Name)
        file2 = "{0}{1}_R2_processed.fastq.gz".format(self.args.WorkingFolder, self.args.Job_Name)
        temp_file1 = Writer(self.log, file1)
        temp_file2 = Writer(self.log, file2)
        self.log.info("Writing {0} and {1}".format(file1, file2))
        fastq1_list = []
        fastq2_list = []
        eof = False
        Read = collections.namedtuple('Read', 'name, seq, index, qual')

        # This generator returns objects not lines.
        while not eof:
            try:
                fastq1_read = next(self.fastq1_file.seq_read())
                fastq2_read = next(self.fastq2_file.seq_read())
                if self.index1_file is not None:
                    index1_read = next(self.index1_file.seq_read())
                if self.index2_file is not None:
                    index2_read = next(self.index2_file.seq_read())
            except StopIteration:
                eof = True
                continue

            current_read_count += 1

            # Apply Filters
            trim_5 = int(self.args.Trim5)
            trim_3 = int(self.args.Trim3)
            min_length = int(self.args.Minimum_Length) + trim_5 + trim_3

            # Filter reads based on length and number of N's.
            if (len(fastq1_read.seq) < min_length or len(fastq2_read.seq) < min_length
                    or fastq1_read.seq.count("N") / len(fastq1_read.seq) >= float(self.args.N_Limit)
                    or fastq2_read.seq.count("N")/len(fastq2_read.seq) >= float(self.args.N_Limit)):
                continue

            # Add the UMT's to the header.
            if self.args.HaloPLEX:
                header1 = "{0}|{0}:{1}".format(index1_read.seq, fastq1_read.name)
                header2 = "{0}|{0}:{1}".format(index1_read.seq, fastq2_read.name)

                # Fixme: This needs to be exposed to the user.
                # Short HaloPLEX reads have issues.  Found that reads <= 100 all show a 3' -1 or -2 error
                if len(fastq1_read.seq) <= 100:
                    read_trim(fastq1_read, trim5=0, trim3=3)
                if len(fastq2_read.seq) <= 100:
                    read_trim(fastq2_read, trim5=0, trim3=3)

            elif self.args.ThruPLEX:
                # header1 = "{0}|{1}".format(fastq1_read.name.split(":")[-1], fastq1_read.name)
                umt1 = fastq1_read.seq[:6]
                umt2 = fastq2_read.seq[:6]
                header1 = "{0}|{1}:{2}".format(umt1, umt2, fastq1_read.name)
                header2 = "{0}|{1}:{2}".format(umt1, umt2, fastq2_read.name)
                read_trim(fastq1_read, trim5=len(umt1), trim3=0)
                read_trim(fastq2_read, trim5=len(umt2), trim3=0)

            elif self.args.FASTQ_PreProcess:
                # The indices are after the last ":" in the header.
                header1 = "{}:{}+{}".format(fastq1_read.name, index1_read.seq, index2_read.seq)
                header2 = "{}:{}+{}".format(fastq2_read.name, index1_read.seq, index2_read.seq)

            else:
                self.log.error("Only HaloPLEX or ThruPLEX currently enabled.")
                raise SystemExit(1)

            # Trim sequences from ends if needed.
            if trim_5 > 0 or trim_3 > 0:
                read_trim(fastq1_read, trim_5, trim_3)
                read_trim(fastq2_read, trim_5, trim_3)

            fastq1_read.name = header1
            fastq2_read.name = header2

            fastq1_list.append(Read(fastq1_read.name, fastq1_read.seq, fastq1_read.index, fastq1_read.qual))
            fastq2_list.append(Read(fastq2_read.name, fastq2_read.seq, fastq2_read.index, fastq2_read.qual))

            # empirically determined for UNC Longleaf cluster.  May need to expose this to user.  Writes blocks of data
            # to disk speeding up entire process.
            if current_read_count % 1000000 == 0:
                temp_file1.write(fastq1_list)
                temp_file2.write(fastq2_list)
                fastq1_list.clear()
                fastq2_list.clear()

        # Cleans up any writes still needed and closes files
        if fastq1_list:
            temp_file1.write(fastq1_list)
            fastq1_list.clear()
        if fastq2_list:
            temp_file2.write(fastq2_list)
            fastq2_list.clear()
        if temp_file1:
            temp_file1.close()
        if temp_file2:
            temp_file2.close()

        self.log.info("Modified FASTQ file(s) written")
        if self.args.FASTQ_PreProcess:
            Tool_Box.compress_files(file1, self.log)
            Tool_Box.compress_files(file2, self.log)
        return file1, file2


class FastqProcessing:
    def __init__(self, args, log, paired_end=False):
        self.log = log
        self.args = args
        self.fq1_batch = None
        self.fq2_batch = None
        self.fastq3_file = None
        self.paired_end = paired_end
        self.master_index_dict = {}
        self.sample_manifest_list = None
        self.sample_manifest_dictionary = collections.defaultdict(list)
        self.demultiplex_file_names = collections.defaultdict(list)
        self.summary_data = collections.defaultdict(list)

    def dataframe_build(self):
        """
        Build the dataframes containing the worker ID, indices, and file names.
        """
        self.log.info("Begin building dataframes")

        with open(self.args.Master_Index_File) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                self.master_index_dict[l_list[0]] = [l_list[1], l_list[2]]

        self.sample_manifest_list = Tool_Box.FileParser.indices(self.log, self.args.SampleManifest)

        outfile_list_dict = collections.defaultdict(list)

        '''
        self.demultiplex_file_names["Unknown"] = \
            [Writer(self.log, "{}{}_Unknown_R1.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "Unknown"),
             Writer(self.log, "{}{}_Unknown_R2.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "Unknown")]

        self.demultiplex_file_names["GhostIndex"] = \
            [Writer(self.log, "{}{}_GhostIndex_R1.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex"),
             Writer(self.log, "{}{}_GhostIndex_R2.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex")]
        '''
        # Initialize our summary dataframe
        self.summary_data["Unknown"] = ["", "", 0]
        self.summary_data["GhostIndex"] = ["", "", 0]
        self.summary_data['total_reads'] = ["", "", 0]

        for i in range(int(self.args.Spawn)):
            outfile_list_dict["Unknown"].append(
                ["{}{}_Unknown_R1.fastq.gz".format(self.args.WorkingFolder, i),
                 "{}{}_Unknown_R2.fastq.gz".format(self.args.WorkingFolder, i)])

            '''
            I noticed with Nextera Libraries we occasionally found an index in a read that was valid but never used
            in the library prep.  I don't know, yet, how prevalent these Ghost Indices are so I am capturing them for
            later analysis.
            '''
            outfile_list_dict["GhostIndex"].append(
                ["{}{}_GhostIndex_R1.fastq.gz".format(self.args.WorkingFolder, i),
                 "{}{}_GhostIndex_R2.fastq.gz".format(self.args.WorkingFolder, i)])

        for sample in self.sample_manifest_list:
            sample_index = sample[0]
            self.sample_manifest_dictionary[sample_index] = [sample[1], sample[2], sample[3]]

            # Check for errors in the sample manifest
            if sample_index not in self.master_index_dict:
                self.log.error("Sample index {} for sample {}, replicate {} is not in Master Index File.".
                               format(sample_index, sample[1], sample[2]))
                raise SystemExit(1)

            for i in range(int(self.args.Spawn)):
                sample_key = "{}|{}".format(i, sample_index)

                # Check for more errors in the sample manifest
                if sample_key in self.demultiplex_file_names:
                    self.log.error("The index {} is duplicated.  Correct the error in {} and try again."
                                   .format(sample_index, self.args.SampleManifest))
                    raise SystemExit(1)

                # Get a list of the temporary FASTQ file names generated.
                outfile_list_dict[sample_index].append(
                    ["{}{}_{}_R1.fastq.gz".format(self.args.WorkingFolder, i, sample_index)])

                # Initialize FASTQ file objects
                # tmp1 = gzip.open("{}{}_{}_R1.fastq.gz".format(self.args.WorkingFolder, i, sample_index), "wb")
                # tmp2 = gzip.open("{}{}_{}_R2.fastq.gz".format(self.args.WorkingFolder, i, sample_index), "wb")
                # tmp1.close()
                # tmp2.close()
                '''
                self.demultiplex_file_names[sample_key] = \
                    [Writer(self.log,
                            "{}{}_{}_R1.fq.gz".format(self.args.WorkingFolder, i, sample_index),
                            sample_index),
                     Writer(self.log,
                            "{}{}_{}_R2.fq.gz".format(self.args.WorkingFolder, i, sample_index),
                            sample_index)]
            '''
            # Add sample tracking framework
            self.summary_data[sample_index] = [sample[1], sample[2], 0]

        self.log.info("Dataframes built")

        return

    def file_writer(self, worker_id, fq1_batch):
        """
        Called by multiprocessor in MainEntry.  Process FASTQ file(s) and write temp version(s).
        :param tmp_file_list:
        :param worker_id:
        :param fq1_batch:
        :return:
        """

        self.output_process(worker_id, FileWriter.file_writer(self, fq1_batch))

    def output_process(self, worker_id, temp_data_dict):
        """
        Write data to the appropriate temporary compressed file.
        :param file_list:
        :param worker_id:
        :param temp_data_dict:
        """

        for sample_index in temp_data_dict:

            tmp_count_file = open("{}{}_{}_mismatch.tmp".format(self.args.WorkingFolder, worker_id, sample_index), "a")
            tmp_count_file.write("{}\t{}\t{}\n"
                                 .format(temp_data_dict[sample_index]["QC"][0], temp_data_dict[sample_index]["QC"][1],
                                         temp_data_dict[sample_index]["QC"][2]))
            tmp_count_file.close()

            outstring = ""
            for block in temp_data_dict[sample_index]["R1"]:
                outstring += "{}".format(block)
            outfile = open("{}{}_{}.tmp".format(self.args.WorkingFolder, worker_id, sample_index), "a")
            outfile.write(outstring)
            outfile.close()

        # Not clearing this results in duplications in the output.  Strange since the dictionary is rebuilt everytime
        # within the job.
        temp_data_dict.clear()


class FastqQuality:
    """
    This class will look at various aspects of the library.  As of Aug. 22 2017 it will analyze ThruPLEX data only.
    """

    # ToDo: Add ability to look at average read quality and positional read qualities grouped by index.

    def __init__(self, args, log, paired_end):
        self.index_list = Tool_Box.FileParser.indices(log, args.Index_File)
        self.index_list.append(("Unknown", "Unknown"))
        self.log = log
        self.args = args
        self.paired_end = paired_end
        self.file1_anchor_seq = "TCAGTAGCTCA"
        self.file2_anchor_seq = "TCAGTAGCTCA"
        self.anchor_dict = None
        self.umt_counts_dict = None

    def module_director(self, splitter_data):
        """
        Parallel job coordination and data processing call.
        """

        self.log.info("Spawning {0} parallel jobs for quality analysis of {1} temporary FASTQ files."
                      .format(self.args.Spawn, len(splitter_data.fastq_file_list)))

        dict_list = []
        data_bundle = int(self.args.prog_check), self.index_list, self.file1_anchor_seq, self.file2_anchor_seq

        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        dict_list += p.starmap(self.quality_check, zip(itertools.repeat(data_bundle), splitter_data.fastq_file_list))

        # Data captured from the multiprocessing pool in this manner is messy.  This sorts it.
        anchor_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        umt_counts_dict = collections.defaultdict(lambda: collections.defaultdict(int))

        for data_dicts in dict_list:
            for dd in data_dicts:
                for index_key in dd:
                    for key2 in dd[index_key]:
                        if key2 in ("R1", "R2"):
                            # Total the anchor seq lengths.
                            r2_list = anchor_dict[index_key]["R2"]
                            r1_list = anchor_dict[index_key]["R1"]
                            anchor_dict[index_key]["R1"] = list(map(add, data_dicts[0][index_key]["R1"], r1_list))
                            anchor_dict[index_key]["R2"] = list(map(add, data_dicts[0][index_key]["R2"], r2_list))

                            # if index_key in anchor_dict:
                            #     r2_list = anchor_dict[index_key]["R2"]
                            #     r1_list = anchor_dict[index_key]["R1"]
                            #     anchor_dict[index_key]["R1"] = list(map(add, data_dicts[0][index_key]["R1"], r1_list))
                            #     anchor_dict[index_key]["R2"] = list(map(add, data_dicts[0][index_key]["R2"], r2_list))
                            # else:
                            #     anchor_dict[index_key]["R1"] = data_dicts[0][index_key]["R1"]
                            #     anchor_dict[index_key]["R2"] = data_dicts[0][index_key]["R2"]
                        else:
                            # Total the UMT count data
                            umt_counts_dict[index_key][key2] += data_dicts[1][index_key][key2]

        self.anchor_dict = anchor_dict
        self.umt_counts_dict = umt_counts_dict
        self.log.info("All parallel quality analysis jobs complete. Begin compiling data.")

        self.data_processing()

    def data_processing(self):
        """
        Gather the data from all the dictionaries and format it for the output files.
        """

        # Build the header string and write it to the output file.
        anchor_outstring = "Sample_Name\tIndex\tTotal_Reads"
        if self.paired_end:
            anchor_outstring = "Sample_Name\tIndex\tTotal_Reads\tPair"

        for i in range(len(self.file1_anchor_seq) + 1):
            anchor_outstring += "\t{0}_Mismatch".format(i)

        # Process anchor sequence data and write file.
        for index in self.index_list:
            for key in natsort.natsorted(self.anchor_dict[index[0]]):
                tmp_string = "\t".join(str(x) for x in self.anchor_dict[index[0]][key])
                total_reads = sum(self.anchor_dict[index[0]][key]) * 0.5

                anchor_outstring += "\n{0}\t{1}\t{2}\t{3}\t{4}".format(index[1], index[0], total_reads, key, tmp_string)

        # Initialize the output file.
        outfile = "{0}{1}_FASTQ_QualityAssessment.txt".format(self.args.WorkingFolder, self.args.Job_Name)
        quality_outfile = open(outfile, "w")
        quality_outfile.write(anchor_outstring)
        quality_outfile.close()
        self.log.info("{0}_FASTQ_QualityAssessment.txt".format(self.args.Job_Name))

        self.anchor_dict.clear()

    @staticmethod
    def quality_check(data_bundle, fastq_files):
        """
        Called by the multiprocessor pool.  Examines the indices and determines the mismatches and N counts.

        :param data_bundle:
        :param fastq_files:
        :return:
        """

        prog_check = data_bundle[0]
        index_list = data_bundle[1]
        file1_anchor_seq = data_bundle[2]
        file2_anchor_seq = data_bundle[3]
        fastq1 = FASTQ_Reader(fastq_files[0])
        fastq2 = FASTQ_Reader(fastq_files[1])

        umt_dict = collections.defaultdict(lambda: collections.defaultdict(int))
        anchor_dict = Tool_Box.VivifiedDictionary()
        read_count = 0

        try:
            while True:
                fastq1_read = next(fastq1.seq_read())
                fastq2_read = next(fastq2.seq_read())
                read_count += 1

                if read_count % int(prog_check) == 0:
                    print("      -->Processed {0} reads in file {1} and {2}."
                          .format(read_count, fastq_files[0], fastq_files[1]))

                # Get read index and UMT.
                umt = "{0}{1}".format(fastq1_read.name.split("|")[0], fastq2_read.name.split("|")[1].split(":")[0])
                read_index = fastq1_read.name.split(":")[-1]

                # Quantify anchor lengths.
                unknown_anchor1 = fastq1_read.seq[7: 18]
                unknown_anchor2 = fastq2_read.seq[7: 18]
                match1 = Levenshtein.distance(file1_anchor_seq, unknown_anchor1)
                match2 = Levenshtein.distance(file2_anchor_seq, unknown_anchor2)

                for index in index_list:
                    index_match = Levenshtein.distance(read_index, index[0][: 6])

                    # Add anchor and UMT data to dictionaries.
                    if index[0] in anchor_dict and index_match < 2:
                        anchor_dict[index[0]]["R1"][match1] += 1
                        anchor_dict[index[0]]["R2"][match2] += 1
                        umt_dict[index[0]][umt] += 1
                        # if umt in umt_dict[index[0]]:
                        #     umt_dict[index[0]][umt] += 1
                        # else:
                        #     umt_dict[index[0]][umt] = 1

                    elif index_match < 2:
                        anchor_dict[index[0]]["R1"] = [0] * len(file1_anchor_seq)
                        anchor_dict[index[0]]["R2"] = [0] * len(file2_anchor_seq)
                        anchor_dict[index[0]]["R1"][match1] += 1
                        anchor_dict[index[0]]["R2"][match2] += 1
                        umt_dict[index[0]][umt] = 1
        except StopIteration:
            return anchor_dict, umt_dict


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
        # self.file = gzip.open(out_file_string, "rb")
        self.log = log

    def lethal_write(self, read):
        """
        This is potentially dead code.
        :param read:
        :return:
        """
        # ToDo: Potentially dead code block.

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
        Writes our new FASTQ file.
        :param read_list:
        :return:
        """
        outstring = ""
        for read in read_list:
            try:
                assert len(read[1]) == len(read[2])
            except AssertionError:
                self.log.error("Sequence and quality scores of different lengths! Read Name {0}; Seq Length {1}; Qual "
                               "Length {2}".format(read[0], len(read[1]), len(read[2])))
                raise SystemExit(1)
            outstring += "@{}\n{}\n+\n{}\n".format(read[0], read[1], read[2])

        self.file.write(outstring)
        read_list.clear()

        return True

    def close(self):
        """
        Closes FASTQ file
        :return:
        """
        self.file.close()
        return True


class FASTQ_Reader:
    """
    Main class that creates FASTQ reads using a generator
    """
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
        Create FASTQ file object.
        :return:
        """
        if len(self.input_file) < 3:
            self.log.warning("FASTQ file parameter missing from options file. Correct error and try again.")
            raise SystemExit(1)

        if not pathlib.Path(self.input_file).is_file():
            self.log.warning("FASTQ file {} not found.  Correct error and run again.".format(self.input_file))
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
            self.log.warning("Unsupported file-type for {}.  Only TEXT or GZIP Allowed.".format(self.input_file))
            raise SystemExit(1)
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
    Provide additional trimming to reads beyond the adaptor trim.
    :param fastq_read:
    :param trim5:
    :param trim3:
    """
    # fastq_read.seq = fastq_read.seq[trim5:-trim3]
    # fastq_read.qual = fastq_read.qual[trim5:-trim3]
    #
    if trim5 and trim3:
        fastq_read.qual = fastq_read.qual[trim5:-trim3]
    elif trim5:
        fastq_read.seq = fastq_read.seq[trim5:]
        fastq_read.qual = fastq_read.qual[trim5:]
    elif trim3:
        fastq_read.seq = fastq_read.seq[:-trim3]
        fastq_read.qual = fastq_read.qual[:-trim3]
