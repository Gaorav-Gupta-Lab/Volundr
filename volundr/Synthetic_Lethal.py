"""
Synthetic_Lethal.py 0.10.0
    April 28, 2019
    Dennis A. Simpson
    Modified so that variable length sgRNA target sequences can be used.
Synthetic_Lethal.py 0.8.0
    May 11, 2018
    Dennis A. Simpson
    Major change to index file structure.  Index file uses Index ID's instead of sequence now.  Requires a master index
    file to lookup sequence.  As part of this the internal logic was changed and the count file names were changed to
    use the Index ID instead of the index sequence.  Because Python was not releasing the memory in a timely fashion, an
    explicit call to the garbage collector was added to clean things up.
Synthetic_Lethal.py 0.7.0
    October 30, 2017
    Dennis A. Simpson
    Fixed a looping bug that was causing target search to run for days.  Refactored code so that the statistics is all
    done in memory now.  Cleaned up many lines of dead code.  This has resulted in a great speed improvement.  Added
    some quality checks to the targets.  Removed the requirement that only one target be allowed per mismatch.  Now
    returns first target with lowest mismatch score found.
Synthetic_Lethal.py 0.6.0
    September 10, 2017
    Dennis A. Simpson
    Target search is done based on the location of an anchor sequence.  If this sequence is not found then no search is
    attempted.  Search no longer uses padded sequences.  A label error has been corrected in the output files.  Target
    search stops when a mismatch of 0 or 1 is found.  Greater mismatches require searching to the end of the target
    file.  Code has been cleaned up.  Smith-Waterman and Muscle options have been removed.
Synthetic_Lethal.py 0.5.0
    July 25, 2017
    Dennis A. Simpson
    Major changes.  Added statistical analysis routines, folded permutation function into class.
Synthetic_Lethal.py 0.2.0
    Dec. 1, 2016
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).  Added version dependency checks for
    modules.
Synthetic_Lethal.py v2.0
    Aug. 21, 2016
    Dennis A. Simpson
    Added the ability to search for targets by Smith-Waterman and Muscle alignments.  Fixed output summary that was not
    aligning the total reads and targeted reads correctly.
Synthetic_Lethal.py v1.6
    Aug. 19, 2016
    Dennis A. Simpson
    Fixed a bug that caused the no index analysis to loop the same number of times as indexes present.  Also corrected
    a bug that over counted the number of reads without targets.
Synthetic_Lethal.py v1.5
    Aug. 7, 2016
    Dennis A. Simpson
    Refactored code.  Made this formally a class now.
Synthetic_Lethal.py v1.2
    April 2, 2016
    Dennis A. Simpson
    Added a control to determine if sequence is reverse complimented prior to target search.  Added a Try/Except near
    line 413 to capture and bypass a rare error seen in the CRISPR searches.
Synthetic_Lethal.py v1.1
    January 20, 2016
    Dennis A. Simpson
    Minor tweak to add a column in the param output file that contains the total targeted reads for a given index.
Synthetic_Lethal.py v1.0
    December 29, 2015
    Dennis A. Simpson
    Program now requires that sequencing indices and shRNA targets all be unique.  Unknown sequence target can range in
    length from the length of the shRNA target to any size defined by the user.  Expected start location and desired
    padding values are available to the user now.  These changes have slowed the execution of the code substantially.
    CPU time for 20,000 reads now approximately 131 seconds using a mismatch of 2 and shRNA padding of 2.  Reads with
    no shRNA target are tabulated.
Synthetic_Lethal.py v0.8
    December 22, 2015
    Dennis A. Simpson
    Still experimenting with the best distance calculator.  Now using python-Levenshtein to do the comparison.  This 
    module is a full 2 seconds faster per 2000 reads than editdistance.  This should shave at least 15 minutes from the 
    cpu time on a full run.  Additional major changes.  Demultiplexed FASTQ can be given as input.  This will analyzed,
    demultiplexed, the demultiplexed files will then be searched for shRNA targets.  Overall these changes have resulted
    in more cpu time.  20,000 reads now take 85 seconds vs 67 seconds.  The trade off is much shorter real world time.
Synthetic_Lethal.py v0.6
    December 18, 2015.
    Dennis A. Simpson
    Changed the match_maker function from using re.finall() to editdistance.eval().  The latter is an implementation of
    Levenshtein distance using C++ and Cython.  Dramatic improvement; 48 seconds for 200 reads down to 1.3 seconds.
Synthetic_Lethal.py v0.5
    December 18, 2015
    Dennis A. Simpson
    This magical little module will score shRNA targets in indexed cell populations.  While it does run, it is slow.  
    I think the code may benefit from changing the mismatch calling from regex to blat.  Will require some fairly
    substantial additions to make it happen.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2018
"""
import gc
import ntpath
import Valkyries.FASTQ_Tools as FASTQ_Tools
import Valkyries.Tool_Box as Tool_Box
import Valkyries.Sequence_Magic as Sequence_Magic
from time import clock
import datetime
import collections
import pathos
import itertools
import os
import natsort
from scipy.stats import gmean
from scipy.stats import ks_2samp
import statistics
import numpy
import math

__author__ = 'Dennis A. Simpson'
__version__ = '0.10.1'
__package__ = 'Volundr'


class SyntheticLethal:
    def __init__(self, log, args):
        if getattr(args, "FASTQ1", False):
            self.fastq = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)
        self.date_format = "%a %b %d %H:%M:%S %Y"
        self.run_start = datetime.datetime.today().strftime(self.date_format)
        self.gene_data_dict = None

        self.control_td_norm_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        self.gtc_norm_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        self.sample_mapping_dict = collections.defaultdict(list)
        self.target_dict = collections.defaultdict(str)
        self.targets = collections.defaultdict(list)
        self.args = args
        self.sample_td_norm_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        self.index_list = Tool_Box.FileParser.indices(log, self.args.Index_File)
        self.fastq_read_counts = [0, 0, 0, 0]  # Tracks number of reads in input file.
        self.log = log
        self.index_dict, self.fastq_file_dict, self.fastq_out_list, self.master_index_dict = \
            SyntheticLethal.dictionary_builds(self)

    def combo_test(self):
        """
        This will combine the replicates done on different days into a single file.  Not terribly user friendly yet.
        :return:
        """

        sample_name = "PolQ"
        input_file = open("{0}{1}_{2}_combined.txt".format(self.args.Working_Folder, self.args.Job_Name, sample_name))
        input_file_dict = collections.defaultdict(float)
        input_file_gene_dict = collections.defaultdict(list)
        control_list = []
        first_line = True

        for l in input_file:
            if first_line:
                first_line = False
                continue

            l_list = [x for x in l.strip("\n").split("\t")]
            gene_name = l_list[0].split("_")[0]
            input_file_dict[l_list[0]] = float(l_list[1])

            if gene_name == self.args.Species:
                control_list.append(float(l_list[1]))

            if gene_name in input_file_gene_dict:
                input_file_gene_dict[gene_name].append(float(l_list[1]))
            else:
                input_file_gene_dict[gene_name] = [float(l_list[1])]

        out_file = open("{0}{1}_KS_Log2_{2}_Combined_Genes.txt"
                        .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
        out_string = "Gene\tLog2\tp-value"

        for gene in natsort.natsorted(input_file_gene_dict):
            v = ks_2samp(control_list, input_file_gene_dict[gene])
            gene_value = statistics.mean(input_file_gene_dict[gene])

            out_string += "\n{0}\t{1}\t{2}".format(gene, gene_value, round(v[1], 6))
        out_file.write(out_string)
        out_file.close()

        # Define some parameters and start the permutations.  Permutations are done based on the control target values.
        working_dict = collections.defaultdict(float)
        permutation_group_size = 10000
        count = 0

        while count < int(permutation_group_size):
            permuted_array = []
            group_key = "m{0}".format(count)
            permuted_array.append(numpy.random.choice(control_list, 10))
            count += 1

            # Process each array of permuted data.
            for permuted_group in permuted_array:
                working_dict[group_key] = sum(permuted_group)/float(len(permuted_group))
        permuted_outstring = "Group\tLog2"
        for group in natsort.natsorted(working_dict):
            permuted_outstring += "\n{0}\t{1}".format(group, working_dict[group])
        permuted_outfile = open("{0}{1}_Permuted_Test.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        permuted_outfile.write(permuted_outstring)
        permuted_outfile.close()

    def fastq_analysis(self):
        """
        This will send the FASTQ file off to be demultiplexed and quantified.  When that is done this will spawn
        parallel jobs to search for target sequences in the demultiplexed FASTQ files.  Each parallel job processes a
        single FASTQ file.
        :return:
        """

        self.fastq_processor()

        self.log.info("Spawning \033[96m{0}\033[m parallel job(s) to search \033[96m{1}\033[m FASTQ files for targets"
                      .format(self.args.Spawn, len(self.fastq_out_list)))

        multiprocessor_tmp_data_list = []
        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        multiprocessor_tmp_data_list.append(
            p.starmap(self.target_search,
                      zip(self.fastq_out_list,
                          itertools.repeat((self.args, self.targets, self.log, self.index_dict)))))

        self.log.info(" ***All Parallel Jobs Complete.***")

        # Process summary data
        self.__summary_output(multiprocessor_tmp_data_list)

        if eval(self.args.Delete_Demultiplexed_FASTQ):
            self.log.debug("Deleting modified FASTQ files from system.")
            Tool_Box.delete(self.fastq_out_list)

        elif eval(self.args.compress):
            self.log.info("Begin compressing FASTQ files with gzip.")
            p = pathos.multiprocessing.Pool(int(self.args.Spawn))
            p.starmap(Tool_Box.compress_files, zip(self.fastq_out_list))

        return

    def statistics(self):

        self.log.info("Begin TCnorm calculations")
        self.tc_norm()
        self.log.debug("TCnorm calculations complete.")

        self.log.info("Begin TDnorm calculations and Log2 transformation")
        self.td_norm()
        self.log.info("TDnorm calculations and Log2 transformation complete.")

        self.control_permutation()

        self.log.info("Begin collapsing gene target groups into individual genes")
        self.__gene_group()
        self.log.info("Gene data manipulations complete.")

        #  Determine significance and output some nice data.
        self.kolmogorov_smirnov()

    def control_permutation(self):
        """
        This bad boy does a permutation analysis using the control targets.  Outputs several data files.
        :param self:
        :return:
        """

        self.log.info("\033[93mBegin Control Permutation Analysis.\033[m")
        big_data_out_string = "Perm_Group\tControl_Target"
        gmean_out_string = "Perm_Group"
        log2_out_string = "Perm_Group"
        raw_data_header = "Control_Target"
        raw_data_string = ""
        selection_space = []
        selection_data_dict = collections.defaultdict(list)
        library_count = 0
        working_library_key_list = []

        # Get the TDnorm data into a dictionary.
        for sample_name in natsort.natsorted(self.control_td_norm_dict):
            if sample_name == "Plasmid" or sample_name == "Unknown":
                continue

            self.log.info("{0} permutation".format(sample_name))
            working_library_key_list.append(sample_name)
            big_data_out_string += "\t{0}".format(sample_name)
            raw_data_header += "\t{0}".format(sample_name)
            library_count += 1

            if sample_name != self.args.Control_Sample:
                gmean_out_string += "\t{0}".format(sample_name)
                log2_out_string += "\t{0}_log2".format(sample_name)

            for control_target_key in self.control_td_norm_dict[sample_name]:
                selection_space.append(control_target_key)
                selection_data_dict[control_target_key].append(self.control_td_norm_dict[sample_name][control_target_key])

        for control_target_key in selection_data_dict:
            control_target_name = self.target_dict[control_target_key]
            raw_data_string += "\n{0}".format(control_target_name)
            for v in selection_data_dict[control_target_key]:
                raw_data_string += "\t{0}".format(math.log2(float(v)))

        # Write Log2 control target data to a file by itself.
        raw_data_out_string = "{0}{1}".format(raw_data_header, raw_data_string)
        raw_data_outfile = open("{0}{1}_Log2_Control_Targets.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        raw_data_outfile.write(raw_data_out_string)
        raw_data_outfile.close()

        # Define some parameters and start the permutations.  Permutations are done based on the control target labels.
        working_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        permutation_group_size = 100000
        count = 0

        while count < int(permutation_group_size):
            permuted_array = []
            group_key = "m{0}".format(count)
            permuted_array.append(numpy.random.choice(selection_space, 10))
            count += 1

            # Process each array of permuted data.
            for permuted_group in permuted_array:
                for control_target_key in permuted_group:

                    big_data_out_string += "\n{0}\t{1}\t{2}" \
                        .format(group_key, control_target_key, "\t".join(map(str, selection_data_dict[control_target_key])))

                    for sample_name in working_library_key_list:
                        try:
                            td_norm_control_ratio = \
                                self.control_td_norm_dict[sample_name][control_target_key]/self.sample_td_norm_dict[self.args.Control_Sample][control_target_key]
                        except ZeroDivisionError:
                            self.log.error("Cannot complete Permutation Analysis.  {} missing from indices."
                                           .format(self.args.Control_Sample))
                            return

                        working_dict[group_key][sample_name].append(td_norm_control_ratio)

        # big_data_outfile = open("{0}{1}_Permuted_BigData.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        # big_data_outfile.write(big_data_out_string)
        # big_data_outfile.close()

        for group_key in natsort.natsorted(working_dict):
            gmean_out_string += "\n{0}".format(group_key)
            log2_out_string += "\n{0}".format(group_key)

            for sample_name in natsort.natsorted(working_dict[group_key]):
                if sample_name != self.args.Control_Sample:
                    gmean_data = gmean(working_dict[group_key][sample_name])
                    gmean_out_string += "\t{0}".format(gmean_data)
                    log2_out_string += "\t{0}".format(round(math.log2(gmean_data), 4))

        # gmean_outfile = open("{0}{1}_Permuted_GMeans.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        # gmean_outfile.write(gmean_out_string)
        # gmean_outfile.close()

        log2_outfile = open("{0}{1}_Permuted_Log2_GMeans.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        log2_outfile.write(log2_out_string)
        log2_outfile.close()

        self.log.info("Permutation Analysis Complete.")

    def kolmogorov_smirnov(self):
        """
        Do a Kolmogorov_Smirnov test on the target sets for each library excluding the no index.  Done on the difference
        with control library.  Writes a file for each library continuing the Log2 delta value and the p-value for each
        gene.

        :return:
        """

        self.log.info("\033[93mBegin Kolmogorov-Smirnov analysis.\033[m")

        for sample_name in self.sample_td_norm_dict:
            if sample_name == "Unknown" or sample_name == "Plasmid" or sample_name == self.args.Control_Sample:
                continue

            # Read TD_norm data for working library into dictionary.
            working_library_dict = collections.defaultdict(list)
            control_list = []

            for target_key in self.sample_td_norm_dict[sample_name]:
                target_name = self.target_dict[target_key]
                gene_name = target_name.split("_")[0]
                log2_td_norm = math.log2(self.sample_td_norm_dict[sample_name][target_key])

                try:
                    log2_control_value = math.log2(self.sample_td_norm_dict[self.args.Control_Sample][target_key])
                except ValueError:
                    self.log.error("{} control sample not present.  Kolmogorov-Smirnov testing aborted."
                                   .format(self.args.Control_Sample))
                    return

                if gene_name == self.args.Species:
                    control_list.append(log2_td_norm-log2_control_value)
                working_library_dict[gene_name].append(log2_td_norm - log2_control_value)

            out_file = open("{0}{1}_{2}_KS_Log2_Delta_Genes.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
            out_string = "Gene\tLog2\tp-value"

            for gene in natsort.natsorted(working_library_dict):
                v = ks_2samp(control_list, working_library_dict[gene])
                gene_value = \
                    round(self.gene_data_dict[gene][sample_name]-self.gene_data_dict[gene][self.args.Control_Sample], 6)
                out_string += "\n{0}\t{1}\t{2}".format(gene, gene_value, round(v[1], 6))
            out_file.write(out_string)
            out_file.close()

        self.log.info("Kolmogorov-Smirnov analysis complete.")

    def __gene_group(self):
        """
        Collapse set of TDnorm target values for each gene into a single Log2 value.  Uses Geometric Mean to collapse
        data set.  Writes a single file containing all the data and creates a dictionary of the date for use later.

        :return:
        """

        out_string = "Gene"
        log_out_string = "Gene"
        gene_data_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        delta_gene_data_dict = collections.defaultdict(lambda: collections.defaultdict(float))

        for sample_name in natsort.natsorted(self.sample_td_norm_dict):

            if sample_name == "Unknown" or sample_name == "Plasmid":
                continue

            tmp_delta_dict = collections.defaultdict(list)
            tmp_dict = collections.defaultdict(list)

            for target_key in self.sample_td_norm_dict[sample_name]:
                target_name = self.target_dict[target_key]
                gene_name = target_name.split("_")[0]

                td_norm = self.sample_td_norm_dict[sample_name][target_key]
                control_td_norm = self.sample_td_norm_dict[self.args.Control_Sample][target_key]
                process_delta = True
                try:
                    tmp_delta_dict[gene_name].append((td_norm / control_td_norm) + 1e-9)
                except ZeroDivisionError:
                    process_delta = False

                tmp_dict[gene_name].append(td_norm + 1e-9)

            for gene in natsort.natsorted(tmp_dict):
                gene_value = math.log2(gmean(tmp_dict[gene]) - 1e-9)
                gene_data_dict[gene][sample_name] = gene_value
                if process_delta:
                    delta_gene_value = math.log2(gmean(tmp_delta_dict[gene]) - 1e-9)
                    delta_gene_data_dict[gene][sample_name] = delta_gene_value

            out_string += "\t{0}".format(sample_name)
            if process_delta or sample_name != self.args.Control_Sample:
                log_out_string += "\t{0}".format(sample_name)

        for gene in natsort.natsorted(gene_data_dict):
            out_string += "\n{0}".format(gene)
            log_out_string += "\n{0}".format(gene)
            for sample_name in natsort.natsorted(gene_data_dict[gene]):
                log_out_string += "\t{0}".format(gene_data_dict[gene][sample_name])

                if process_delta:
                    out_string += "\t{0}".format(delta_gene_data_dict[gene][sample_name])

        if process_delta:
            out_file = \
                open("{0}{1}_Log2_Delta_{2}_Genes.txt".format(self.args.Working_Folder, self.args.Job_Name,
                                                              self.args.Control_Sample), "w")
            out_file.write(out_string)
            out_file.close()

        log_out_file = open("{0}{1}_Log2_Genes.txt".format(self.args.Working_Folder, self.args.Job_Name), "w")
        log_out_file.write(log_out_string)
        log_out_file.close()

        if not process_delta:
            self.log.error("Cannot determine control sample {0} deltas.  {0} sample missing"
                           .format(self.args.Control_Sample))

        self.gene_data_dict = gene_data_dict

    def individual_td_norm_and_gene(self):
        """
        This is currently dead code.  I want to create heatmaps with the individual data eventually.  This may or may
        not hold that code.
        :return:
        """
        # ToDo: Remove this block if we are not going to use it!!!!
        Tool_Box.debug_messenger("This is Intended to be DEAD CODE.")

        control_data_file = open("{0}{1}_Plasmid_gTC_norm.txt".format(self.args.Working_Folder, self.args.Job_Name))
        plasmid_dict = {}
        first_line = True

        for control_file_line in control_data_file:
            control_file_line_list = [x for x in control_file_line.strip("\n").split("\t")]
            if first_line:
                first_line = False
            else:
                plasmid_dict[control_file_line_list[0]] = float(control_file_line_list[2])

        for sample_name in self.gtc_norm_dict:

            if sample_name == "Unknown":
                continue

            tc_norm_data_file = open("{0}{1}_{2}_TC_norm.txt"
                                     .format(self.args.Working_Folder, self.args.Job_Name, sample_name))

            first_line = True
            for data_file_line in tc_norm_data_file:

                if first_line:
                    first_line = False
                    continue

                data_file_line_list = [x for x in data_file_line.strip("\n").split("\t")]
                target_name = data_file_line_list[0]

                self.individual_td_norm_dict[sample_name][target_name] = \
                    float(data_file_line_list[2])/plasmid_dict[data_file_line_list[0]]

        tmp_dict = collections.defaultdict(list)
        # Fixme: I don't think this block does anything.
        for library_index in natsort.natsorted(self.individual_td_norm_dict):
            for target_name in natsort.natsorted(self.individual_td_norm_dict[library_index]):
                gene_name = target_name.split("_")[0]

                if gene_name in tmp_dict:
                    if self.individual_td_norm_dict[library_index][target_name] == 0:
                        pass
                    if library_index in tmp_dict[gene_name]:
                        tmp_dict[gene_name].append(self.individual_td_norm_dict[library_index][target_name])

    def td_norm(self):
        """
        Processes the data in the gtc_norm_dict to produce the td_norm data.  Writes the data to a file.
        :return:
        """

        for sample_name in self.gtc_norm_dict:

            if sample_name == "Unknown":
                continue

            self.log.debug("TDnorm for {0}". format(sample_name))

            out_string = "Gene\tTarget\tgTC_norm for {0} mismatches\tTD_norm for {0} mismatches\tLog2_TD_norm"\
                .format(self.args.Target_Mismatch)

            for target_key in collections.OrderedDict(natsort.natsorted(self.target_dict.items(), key=lambda t: t[1])):

                sample_gtc_norm = self.gtc_norm_dict[sample_name][target_key]
                plasmid_gtc_norm = self.gtc_norm_dict["Plasmid"][target_key]
                target_name = self.target_dict[target_key]

                out_string += "\n{0}\t{1}\t{2}".format(target_name, target_key, sample_gtc_norm)

                try:
                    td_norm = sample_gtc_norm/plasmid_gtc_norm
                except ZeroDivisionError:
                    td_norm = 1

                gene_name = target_name.split("_")[0]

                if gene_name == self.args.Species:
                    self.control_td_norm_dict[sample_name][target_key] = td_norm

                out_string += "\t{0}\t{1}".format(td_norm, math.log2(td_norm))
                self.sample_td_norm_dict[sample_name][target_key] = td_norm

            out_file = open("{0}{1}_{2}_TD_norm.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
            out_file.write(out_string)
            out_file.close()

    def tc_norm(self):
        # grouped_data_dict = collections.defaultdict(int)
        library_index_target_counts = collections.defaultdict(int)
        sample_tc_data = collections.defaultdict(float)
        library_tc_norm_values = collections.defaultdict(list)
        bad_targets_dict = collections.defaultdict(int)
        bad_targets_list = []

        # Check plasmid library for drop outs.
        for library_index in self.sample_mapping_dict["Plasmid"]:
            try:
                tmp_data_file = open("{0}{1}_{2}_target_counts.txt"
                                     .format(self.args.Working_Folder, self.args.Job_Name, library_index))
            except FileNotFoundError:
                self.log.error("{0}_{1}_target_counts.txt not found".format(self.args.Job_Name, library_index))
                continue
            first_line = True

            # Go through each target in Plasmid target counts file.
            for line in tmp_data_file:
                if first_line:
                    first_line = False
                    continue

                line_list = [x for x in line.strip("\n").split("\t")]
                target_key = line_list[1]

                # Count the number of reads for each target in each file for each mismatch.
                for i in range(int(self.args.Target_Mismatch) + 1):
                    bad_targets_dict[target_key] += int(line_list[i + 2])
            tmp_data_file.close()

        # Check for missing Plasmid Targets.
        for target_key in bad_targets_dict:
            if target_key not in self.target_dict:
                self.log.warning("{} not found in {}.  Confirm correct Target File is being used."
                                 .format(target_key, self.args.Target_File))
            if bad_targets_dict[target_key] == 0:
                bad_targets_list.append(target_key)
                self.log.warning("{}  {} has no counts in Plasmid".format(self.target_dict[target_key], target_key))

        # Process each library.
        for sample_name in self.sample_mapping_dict:
            if sample_name == "Unknown":
                continue

            for library_index in self.sample_mapping_dict[sample_name]:
                self.log.debug("processing {} sample {}".format(sample_name, library_index))

                try:
                    tmp_data_file = open("{0}{1}_{2}_target_counts.txt"
                                         .format(self.args.Working_Folder, self.args.Job_Name, library_index))
                except FileNotFoundError:
                    self.log.error("{0}_{1}_target_counts.txt not found".format(self.args.Job_Name, library_index))
                    continue

                first_line = True
                target_data_list = []

                # Go through each target in target counts file.
                for line in tmp_data_file:
                    if first_line:
                        first_line = False
                        continue

                    line_list = [x for x in line.strip("\n").split("\t")]

                    target_key = line_list[1]
                    target_data_list.append(line_list)

                    # Count the number of reads for each target in each file for each mismatch.
                    for i in range(int(self.args.Target_Mismatch)+1):
                        # grouped_data_dict[target_key] += int(line_list[i+2])
                        library_index_target_counts[target_key] += int(line_list[i+2])

                    # Adding 1 to prevent errors from 0 counts.
                    library_index_target_counts[target_key] += 1
                tmp_data_file.close()

                library_index_total_count = sum(library_index_target_counts.values())
                # Normalize for reads in individual libraries
                for target_key in library_index_target_counts:
                    library_tc_norm_values[target_key].append(library_index_target_counts[target_key]/library_index_total_count)
                library_index_target_counts.clear()

            # Determine the gTC_norm for each sample
            for target_key in natsort.natsorted(library_tc_norm_values):
                # Skip bad guides
                bad_target = False
                if len(bad_targets_list) > 0:
                    for k in bad_targets_list:
                        if k == target_key:
                            bad_target = True
                    if bad_target:
                        continue

                gtc_norm = statistics.mean(library_tc_norm_values[target_key])

                if gtc_norm == 0:
                    self.log.warning("gtc_norm=0|{} {} {}".format(sample_name, target_key, sample_tc_data[target_key]))
                self.gtc_norm_dict[sample_name][target_key] = gtc_norm
            # grouped_data_dict.clear()
            library_tc_norm_values.clear()

    def __summary_output(self, multiprocessor_tmp_data_list):
        """
        Process temporary data file into final output.
        :return:
        """

        self.log.debug("Begin Writing Data Output File")
        args = self.args

        index_header = "Index Name\tSample Name\tSample Replica\tTotal\tShort Reads\tFull Reads"

        for i in range(int(args.Index_Mismatch)+1):
            index_header += "\t{0}_mismatches".format(i)
        index_header += "\tTargeted\tNot Targeted\tFraction Targeted\n"

        run_param_out = open("{0}{1}_summary.txt".format(args.Working_Folder, args.Job_Name), "w")
        run_stop = datetime.datetime.today().strftime(self.date_format)
        run_param_out.write("Running:\t{0} Synthetic_Lethal v{1}\nStart_Time:\t{2}\nStop_Time\t{3}\nFASTQ_File:\t{4}\n"
                            "INDEX_File:\t{5}\nTarget_File:\t{6}\nIndex_Mismatches\t{7}\nTarget_Mismatches\t{8}\n"
                            "Target_Padding\t{9}\nExpected_Position\t{10}\nMin_Read_length\t{11}\nTarget_Start\t{12}\n"
                            "Target_Length\t{13}\nTotal_Reads:\t{14}\nIndexed_Reads:\t{15}\nUnknown_Count:\t{16}\n"
                            "Unknown_Short_Count:\t{17}\nUnknown_Full_Length_Count:\t{18}\n\n{19}"
                            .format(__package__, __version__, self.run_start, run_stop, args.FASTQ1, args.Index_File,
                                    args.Target_File, args.Index_Mismatch, args.Target_Mismatch, args.Target_Padding,
                                    args.Expected_Position, args.Min_Length, args.Target_Start, args.Target_Length,
                                    self.fastq_read_counts[0], self.fastq_read_counts[0] - self.fastq_read_counts[1],
                                    self.fastq_read_counts[1], self.fastq_read_counts[2], self.fastq_read_counts[3],
                                    index_header))

        # Data captured from the multiprocessor calls is a list of lists.  This is to remove the outer level without
        # the need to modify my existing code below.
        multiprocessor_data_list = []

        for b in multiprocessor_tmp_data_list:
            for c in b:
                multiprocessor_data_list.append(c)

        param_string = ""
        for data in multiprocessor_data_list:

            # If we are not analyzing the unknowns there will be a NoneType entry in the list.
            try:
                index_name = data[1]
            except TypeError:
                continue

            index_seq = self.master_index_dict[index_name]
            # Index Name, Sample Name, Sample Replica
            param_string += "{}\t{}\t{}"\
                .format(index_name, self.index_dict[index_seq][2], self.index_dict[index_seq][3])

            # Mismatched counts
            for i in range(len(self.index_dict[index_seq][0])):
                param_string += "\t{0}".format(self.index_dict[index_seq][0][i])

            # Targeted, Not Targeted, Fraction Targeted.
            param_string += "\t{}\t{}\t{}\n"\
                .format(data[0], int(self.index_dict[index_seq][0][1]) - int(data[0]),
                        int(data[0])/int(self.index_dict[index_seq][0][1]))

        run_param_out.write(param_string)
        run_param_out.close()

        self.log.debug("Data Summary File Written.")

        return

    def dictionary_builds(self):
        """
        Build dictionaries, initialize output files and capture list of output file names.
        :return:
        """

        # Read master index file into a dictionary.
        fastq_out_list = []
        master_index_dict = {}
        fastq_file_dict = collections.defaultdict(object)
        index_dict = {"Unknown": [[0] * (int(self.args.Index_Mismatch) + 4), "Unknown", "Unknown", "Unknown"]}

        with open(self.args.Master_Index_File) as f:
            for l in f:
                if "#" in l or len(l) == 0:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                master_index_dict[l_list[0]] = l_list[1]

        for sample in self.index_list:
            index_sequence = master_index_dict[sample[0]]
            if index_sequence in index_dict:
                self.log.error("The index {0} is duplicated.  Correct the error in {1} and try again."
                               .format(sample[0], self.args.Index_File))
                raise SystemExit(1)

            sample.insert(0, [0]*(int(self.args.Index_Mismatch)+4))
            # for each sample name append a list of all index ID's
            self.sample_mapping_dict[sample[2]].append(sample[1])

            index_dict[index_sequence] = sample
            if eval(self.args.Target_Search):
                fastq_file_dict[index_sequence] = \
                    FASTQ_Tools.Writer(self.log, "{0}{1}_{2}.fastq"
                                           .format(self.args.Working_Folder, self.args.Job_Name, index_sequence))

                fastq_out_list.append("{0}{1}_{2}.fastq"
                                          .format(self.args.Working_Folder, self.args.Job_Name, index_sequence))

        # This is for no index found.
        self.index_list.append(("Unknown", "Unknown", "Unknown"))

        if eval(self.args.Analyze_Unknowns):
            master_index_dict["Unknown"] = "Unknown"
            fastq_file_dict["Unknown"] = \
                FASTQ_Tools.Writer(self.log, "{0}{1}_Unknown.fastq".format(self.args.Working_Folder, self.args.Job_Name))

            fastq_out_list.append("{0}{1}_Unknown.fastq".format(self.args.Working_Folder, self.args.Job_Name))

        # Fill target list and dictionary.  Do initial quality check on target file for duplicates.
        target_list = []

        for target in Tool_Box.FileParser.indices(self.log, self.args.Target_File):
            try:
                target_key = target[2][int(self.args.Target_Start):][:int(self.args.Target_Length)]
            except ValueError:
                target_key = target[2]

            target_name = target[1]
            target_list.append(target_key)
            self.targets[len(target_key)].append(target_key)

            if target_key in self.target_dict:
                self.log.error("The target sequence in {0} is duplicated.  Correct the error in {1} and try again."
                               .format(target, self.args.Target_File))
                raise SystemExit(1)
            elif target_name in self.target_dict[target_key]:
                self.log.error("Target name in {0} is duplicated.  Correct the error in {1} and try again."
                               .format(target, self.args.Target_File))
                raise SystemExit(1)

            self.target_dict[target_key] = target_name

        similarity_count = 0
        off_target_count = 0
        # Do a fine scale quality analysis of targets looking for similar sequences and off target guides.
        for target_key in self.target_dict:
            for target in target_list:
                mismatch_index = Sequence_Magic.match_maker(target, target_key)
                if 0 < mismatch_index <= 3:
                    target_gene_name = self.target_dict[target].split("_")[0]
                    query_gene_name = self.target_dict[target_key].split("_")[0]

                    if target_gene_name != query_gene_name:
                        off_target_count += 1
                        self.log.debug("!!!POTENTIAL OFF TARGET!!! {0} differs from {1} by {2}"
                                       .format(self.target_dict[target_key], self.target_dict[target], mismatch_index))
                    else:
                        similarity_count += 1
                        self.log.debug("{0} {1} differs from {2} {3} by only {4}"
                                       .format(self.target_dict[target_key], target_key, self.target_dict[target],
                                               target, mismatch_index))
        if similarity_count > 0:
            self.log.info("{0} targets similar to each other".format(similarity_count))
            self.log.info("{0} potential off target guides".format(off_target_count))

        self.log.info("Dictionaries Built")
        return index_dict, fastq_file_dict, fastq_out_list, master_index_dict

    def fastq_processor(self):
        """
        Demultiplex FASTQ file.  Writes new files every 1.5 million reads processed.
        :return:
        """

        self.log.info("\033[96mDemultiplexing Input FASTQ File.\033[m|Writing Data to \033[93m{0}\033[m Files."
                      .format(len(self.fastq_out_list)))
        args = self.args
        t0 = clock()
        eof = False
        output_dict = collections.defaultdict(list)
        index_count = len(self.index_dict)

        while not eof:
            try:
                fastq_read = next(self.fastq.seq_read())
            except StopIteration:
                # Clean up any reads not written to files yet.
                for index_seq in output_dict:
                    self.fastq_file_dict[index_seq].write(output_dict[index_seq])
                output_dict.clear()
                eof = True
                continue

            self.fastq_read_counts[0] += 1
            index_found = False  # Flag for marking a read as having a known index.

            if self.fastq_read_counts[0] % 1500000 == 0:
                for index_seq in output_dict:
                    self.fastq_file_dict[index_seq].write(output_dict[index_seq])
                output_dict.clear()
                gc.collect()
                self.log.info("Processed 1,500,000 reads in {} seconds for a total of {:,} reads in file {}"
                              .format((clock() - t0), self.fastq_read_counts[0], ntpath.basename(self.fastq.file_name)))

                t0 = clock()

                # Tool_Box.debug_messenger("Limiting reads here to 1.5 million")
                # self.log.warning("Limiting reads here to 1.5 million")
                # eof = True
            index_loop_count = 0
            while not index_found:
                for index_seq in self.index_dict:
                    index_loop_count += 1

                    # No need to search for the "Unknown" key.
                    unknown = True
                    mismatch_index = int(args.Index_Mismatch)+1
                    if index_seq != "Unknown":
                        unknown = False
                        index_key_length = len(index_seq)
                        mismatch_index = Sequence_Magic.match_maker(index_seq, fastq_read.seq[:index_key_length])

                    if mismatch_index <= int(args.Index_Mismatch) and not unknown:
                        index_found = True
                        self.index_dict[index_seq][0][0] += 1  # Total reads with index.

                        # Check read length and only output if the reads are long enough.
                        if len(fastq_read.seq) > int(args.Min_Length):
                            FASTQ_Tools.read_trim(fastq_read, trim5=index_key_length)
                            fastq_read.name = "{0}|{1}".format(index_seq, fastq_read.name)
                            Read = collections.namedtuple('Read', 'name, seq, index, qual')
                            new_read = Read(fastq_read.name, fastq_read.seq, fastq_read.index, fastq_read.qual)
                            output_dict[index_seq].append(new_read)
                            self.index_dict[index_seq][0][1] += 1  # Full reads with index.
                            self.index_dict[index_seq][0][3 + mismatch_index] += 1  # Mismatch count.
                            break
                        else:
                            self.index_dict[index_seq][0][2] += 1  # Short reads with index.
                        # Break the loop if nothing matches.
                    elif index_loop_count == index_count:
                        index_found = True
                        self.fastq_read_counts[1] += 1  # Total Unknown index reads in FASTQ file.
                        self.index_dict["Unknown"][0][0] += 1  # Total reads without index.

                        if len(fastq_read.seq) > int(args.Min_Length):
                            self.fastq_read_counts[3] += 1  # Full reads with no index.
                            self.index_dict["Unknown"][0][1] += 1  # Full reads with no index.
                            if eval(self.args.Analyze_Unknowns):
                                head_seq = fastq_read.seq[:index_key_length]
                                fastq_read.name = "Unknown|{0}|{1}".format(head_seq, fastq_read.name)
                                FASTQ_Tools.read_trim(fastq_read, trim5=index_key_length)
                                Read = collections.namedtuple('Read', 'name, seq, index, qual')
                                new_read = Read(fastq_read.name, fastq_read.seq, fastq_read.index, fastq_read.qual)
                                output_dict["Unknown"].append(new_read)
                        else:
                            self.fastq_read_counts[2] += 1  # Short reads with no index.
                            self.index_dict["Unknown"][0][2] += 1  # Short reads with no index.

        # Close all the open FASTQ files before proceeding.
        for index_seq in self.fastq_file_dict:
            self.fastq_file_dict[index_seq].close()

        self.log.info("Demultiplexing Complete. Processed a total of {:,} reads in {}"
                      .format(self.fastq_read_counts[0], ntpath.basename(self.fastq.file_name)))

    @staticmethod
    def target_search(fq_file, argvs):
        """
        I intend this to be called using multiprocessor.  Searches a demultiplexed FASTQ file for target sequences.
        :param fq_file:
        :param argvs:
        :return:
        """

        args, targets_dict, log, index_dict = argvs
        log.info("Begin Target Search in {}".format(ntpath.basename(fq_file)))

        # If the FASTQ file is missing we need to get out of here
        if not os.path.isfile(fq_file):
            log.warning("\033[1;31m{0} file not found for target search.\033[m" .format(fq_file))
            return

        # If the FASTQ file is empty remove it and get out of here.
        elif os.stat(fq_file).st_size < 50:
            log.warning("{0} is empty; File removed." .format(fq_file))
            os.remove(str(fq_file))
            return

        t0 = clock()

        # Retrieve the index sequence and index name for the file we are processing.
        index_key = ntpath.basename(fq_file).split(".")[0].split("_")[-1]
        index_name = index_dict[index_key][1]
        multiprocessor_tmp_data_list = [0, index_name]
        index_key_length = len(index_key)

        if not eval(args.Analyze_Unknowns) and index_key == "Unknown":
            log.info("\033[1;31mNotice:\033[m  Skipping {0} at user request.".format(ntpath.basename(fq_file)))
            return multiprocessor_tmp_data_list

        target_data_outstring = "Target\tTarget_Key"
        for i in range(int(args.Target_Mismatch)+1):
            target_data_outstring += "\t{0}_mismatches".format(i)

        # Iterate the FASTQ file; extract the target region; reverse-compliment it; check it against the target list
        # for matches; tabulate results.

        total_target_pos_list = []
        target_found_pos_list = []
        no_target_pos_list = []
        fastq_read_count = 0
        target_count = 0
        no_anchor_count = 0
        anchor_count = 0
        eof = False
        fastq = FASTQ_Tools.FASTQ_Reader(fq_file, log)
        target_count_dict = collections.defaultdict(lambda: collections.defaultdict(int))
        target_file = Tool_Box.FileParser.indices(log, args.Target_File)

        while not eof:
            try:
                fastq_read = next(fastq.seq_read())
            except StopIteration:
                eof = True
                continue

            fastq_read_count += 1
            if len(fastq_read.seq) <= int(args.Min_Length):
                continue

            anchor_found = False
            start_pos = int(args.AnchorStart) - len(index_key)

            while not anchor_found:
                mismatch_index = Sequence_Magic.match_maker(
                    args.AnchorSeq, fastq_read.seq[start_pos:][:len(args.AnchorSeq)])

                # If we do not find the anchor sequence exit the loop and go to the next read.
                if start_pos > int(args.AnchorStop):
                    unknown_seq_start = int(args.Expected_Position) - int(args.Target_Padding)
                    # unknown_seq_start = int(args.Expected_Position) - int(args.Target_Padding) - len(index_key)
                    no_anchor_count += 1
                    break

                elif mismatch_index <= int(args.AnchorMismatch):
                    anchor_found = True
                    unknown_seq_start = start_pos + len(args.AnchorSeq)
                    anchor_count += 1
                    total_target_pos_list.append(unknown_seq_start)

                start_pos += 1

            target_seq, mismatch_index = \
                SyntheticLethal.__position_search(targets_dict, fastq_read, int(args.Target_Mismatch),
                                                  unknown_seq_start, args)

            if target_seq:
                target_count_dict[target_seq][mismatch_index] += 1
                target_found_pos_list.append(unknown_seq_start)
                target_count += 1
            else:
                no_target_pos_list.append(unknown_seq_start)

            if fastq_read_count % 250000 == 0:
                log.info("Searched 250,000 reads in {} seconds for a total of {:,} reads in file {}"
                          .format(int((clock() - t0)), fastq_read_count, ntpath.basename(fq_file)))
                t0 = clock()

        # Process frequency data and write output file.
        log.debug("Processing data for {}".format(ntpath.basename(fq_file)))
        freq_pos_outstring = "Total_Anchors\tPosition\tCount\tFrequency"
        total_target_pos_counter = collections.Counter(total_target_pos_list)
        for k in total_target_pos_counter.items():
            freq_pos_outstring += "\n{0}\t{1}\t{2}".format(k[0]+index_key_length, k[1], round((k[1]/anchor_count), 4))

        freq_pos_outstring += "\nTargets_Found\tPosition\tCount\tFrequency"
        target_found_pos_counter = collections.Counter(target_found_pos_list)
        for k in target_found_pos_counter.items():
            freq_pos_outstring += "\n{0}\t{1}\t{2}".format(k[0]+index_key_length, k[1], round((k[1]/anchor_count), 4))

        freq_pos_outstring += "\nNo_Targets_Found\tPosition\tCount\tFrequency"
        no_target_pos_counter = collections.Counter(no_target_pos_list)
        for k in no_target_pos_counter.items():
            freq_pos_outstring += "\n{0}\t{1}\t{2}".format(k[0]+index_key_length, k[1], round((k[1]/anchor_count), 4))

        target_position_freq_outfile = open("{0}{1}_{2}_Target_Position_Freq.txt"
                                            .format(args.Working_Folder, args.Job_Name, index_name), "w")
        target_position_freq_outfile.write(freq_pos_outstring)
        target_position_freq_outfile.close()

        # Format target count data for output file and write data to file.
        for line in target_file:
            target_name = line[1]
            if args.Target_Length == 'Variable':
                target_key = line[2]
            else:
                target_key = line[2][int(args.Target_Start):][:int(args.Target_Length)]
            target_data_outstring += "\n{0}\t{1}".format(target_name, target_key)
            for i in range(int(args.Target_Mismatch)+1):
                target_data_outstring += "\t{0}".format(target_count_dict[target_key][i])

        target_data_file_name = "{0}{1}_{2}_target_counts.txt".format(args.Working_Folder, args.Job_Name, index_name)
        target_data_out = open(target_data_file_name, "w")
        target_data_out.write(target_data_outstring)
        target_data_out.close()

        log.info(target_data_file_name)
        multiprocessor_tmp_data_list[0] = target_count

        return multiprocessor_tmp_data_list

    @staticmethod
    def __position_search(targets_dict, fastq_read, target_mismatch, unknown_seq_start, args):
        """
        Does a Levenshtein search for items in target list.
        :return:
        :param targets_dict:
        :param fastq_read:
        :param target_mismatch:
        """

        targets_found_dict = collections.defaultdict(list)

        # Go through targets based on size.
        for target_length in targets_dict:

            if eval(args.RevComp):
                unknown_seq = Sequence_Magic.rcomp(fastq_read.seq[unknown_seq_start:][:target_length])
            else:
                unknown_seq = fastq_read.seq[unknown_seq_start:][:target_length]

            for target in targets_dict[target_length]:
                mismatch_index = Sequence_Magic.match_maker(target, unknown_seq)

                if mismatch_index <= 1:
                    return target, mismatch_index

                elif mismatch_index <= target_mismatch:
                    targets_found_dict[mismatch_index].append(target)

        if targets_found_dict:
            for i in range(2, target_mismatch+1):
                if i in targets_found_dict:
                    return targets_found_dict[i][0], i

        return False, False
