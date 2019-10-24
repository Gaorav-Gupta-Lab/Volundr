"""
Synthetic_Lethal.py 2.0.0
    August 30, 2019
    Added multiple sample p-value correction.  Added percentile output.  Added output file for masked sgRNA sequences.
    Added a library control and sample control option.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import gc
import ntpath
from time import clock
import datetime
import collections
import itertools
import os
import statistics
import statsmodels.stats.multitest as stats
import math
import numpy
from scipy.stats import gmean
from scipy.stats import ks_2samp
import natsort
import pathos
import Valkyries.FASTQ_Tools as FASTQ_Tools
import Valkyries.Tool_Box as Tool_Box
import Valkyries.Sequence_Magic as Sequence_Magic

__author__ = 'Dennis A. Simpson'
__version__ = '2.0.0'
__package__ = 'Völundr'


class SyntheticLethal:
    """
    Main class that coordinates the target searching and analysis.
    """
    def __init__(self, log, args):
        """

        :param log:
        :param args:
        """
        if getattr(args, "FASTQ1", False):
            self.fastq = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)
        self.date_format = "%a %b %d %H:%M:%S %Y"
        self.run_start = datetime.datetime.today().strftime(self.date_format)
        self.gene_data_dict = None

        self.control_td_norm_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        self.gtc_norm_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        self.tc_norm_dict = collections.defaultdict(lambda: collections.defaultdict(list))
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

    def fastq_analysis(self):
        """
        This will send the FASTQ file off to be demultiplexed and quantified.  When that is done this will spawn
        parallel jobs to search for target sequences in the demultiplexed FASTQ files.  Each parallel job processes a
        single FASTQ file.
        """

        self.fastq_processor()
        if not self.args.Delete_Demultiplexed_FASTQ and self.args.Compress:
            self.log.info("Begin compressing FASTQ files with gzip.")
            p = pathos.multiprocessing.Pool(int(self.args.Spawn))
            p.starmap(Tool_Box.compress_files, zip(self.fastq_out_list))

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

        if self.args.Delete_Demultiplexed_FASTQ:
            self.log.debug("Deleting modified FASTQ files from system.")
            Tool_Box.delete(self.fastq_out_list)

        elif self.args.Compress:
            self.log.info("Begin compressing FASTQ files with gzip.")
            p = pathos.multiprocessing.Pool(int(self.args.Spawn))
            p.starmap(Tool_Box.compress_files, zip(self.fastq_out_list))

    def statistics(self):
        """
        Runs through the methods to do the analysis in the correct order.
        """
        self.log.info("\033[93mBegin TCnorm calculations\033[m")
        self.tc_norm()
        self.log.info("TCnorm calculations complete.")

        self.log.info("\033[93mBegin TDnorm calculations and Log2 transformation\033[m")
        self.td_norm()
        self.log.info("TDnorm calculations and Log2 transformation complete.")

        self.control_permutation()

        self.log.info("\033[93mBegin collapsing gene target groups into individual genes\033[m")
        self.gene_group()
        self.log.info("Gene data manipulations complete.")

        #  Determine significance and output some nice data.
        # self.kolmogorov_smirnov2()

        self.kolmogorov_smirnov3()

    def control_permutation(self):
        """
        This bad boy does a permutation analysis using the control targets.  Outputs several data files.
        :param self:
        :return:
        """

        self.log.info("\033[93mBegin Control Permutation Analysis.\033[m")

        log2_out_string = \
            "Sample Permutations\nFile Generated:\t{}\nControl Sample:\t{}\nLibrary Control:\t{}\n" \
            "Upper Percentile:\t{}\nLower Percentile:\t{}\n\n"\
            .format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.args.Control_Sample,
                    self.args.Library_Control, self.args.UpperPercentile, self.args.LowerPercentile)

        raw_data_header = "Control_Target"
        raw_data_string = ""
        selection_space = []
        selection_data_dict = collections.defaultdict(list)
        library_count = 0
        working_library_key_list = []

        # Get the TDnorm data into a dictionary.
        for sample_name in natsort.natsorted(self.control_td_norm_dict):
            if sample_name in ["Unknown", "Plasmid", self.args.Library_Control]:
                continue

            self.log.info("{0} permutation".format(sample_name))
            working_library_key_list.append(sample_name)
            raw_data_header += "\t{0}".format(sample_name)
            library_count += 1

            for control_target_key in self.control_td_norm_dict[sample_name]:
                selection_space.append(control_target_key)
                selection_data_dict[control_target_key].append(self.control_td_norm_dict[sample_name]
                                                               [control_target_key])

        for control_target_name in selection_data_dict:
            raw_data_string += "\n{0}".format(control_target_name)
            for v in selection_data_dict[control_target_name]:
                raw_data_string += "\t{0}".format(math.log2(float(v)))

        # Write Log2 control target data to a file by itself.
        raw_data_out_string = "{0}{1}".format(raw_data_header, raw_data_string)
        raw_data_outfile = \
            open("{0}{1}_TDnorm_Log2_Control_Targets.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        raw_data_outfile.write(raw_data_out_string)
        raw_data_outfile.close()

        # Define some parameters and start the permutations.  Permutations are done based on the control target labels.
        working_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        percentile_dict = collections.defaultdict(list)
        permutation_group_size = int(self.args.PermutationCount)
        count = 0

        while count < int(permutation_group_size):
            permuted_array = []
            group_key = "Iteration_{0}".format(count)
            permuted_array.append(numpy.random.choice(selection_space, 10))
            count += 1

            # Process each array of permuted data.
            for permuted_group in permuted_array:
                for control_target_key in permuted_group:

                    for sample_name in working_library_key_list:
                        try:
                            td_norm_control_ratio = \
                                self.control_td_norm_dict[sample_name][control_target_key]/self.sample_td_norm_dict[self.args.Control_Sample][control_target_key]
                        except ZeroDivisionError:
                            self.log.error("Cannot complete Permutation Analysis.  {} missing from indices."
                                           .format(self.args.Control_Sample))
                            return

                        working_dict[group_key][sample_name].append(td_norm_control_ratio)

        log2_perm_data_string = ""
        for group_key in natsort.natsorted(working_dict):
            log2_perm_data_string += "\n{0}".format(group_key)

            for sample_name in natsort.natsorted(working_dict[group_key]):

                gmean_data = gmean(working_dict[group_key][sample_name])
                log2_perm_data_string += "\t{}".format(round(math.log2(gmean_data), 4))
                percentile_dict[sample_name].append(math.log2(gmean_data))

        sample_name_list = []
        upper_limit_list = []
        lower_limit_list = []
        for sample_name in natsort.natsorted(percentile_dict):
            lower_limit = \
                str(numpy.percentile(numpy.array(percentile_dict[sample_name]), float(self.args.LowerPercentile), interpolation='linear'))
            upper_limit = \
                str(numpy.percentile(numpy.array(percentile_dict[sample_name]), float(self.args.UpperPercentile), interpolation='linear'))
            sample_name_list.append(sample_name)
            lower_limit_list.append(lower_limit)
            upper_limit_list.append(upper_limit)

        log2_out_string += "Sample:\t{}\n".format("\t".join(sample_name_list))
        log2_out_string += "Upper Limit:\t{}\n".format("\t".join(upper_limit_list))
        log2_out_string += "Lower Limit:\t{}\n\n{}".format("\t".join(lower_limit_list), log2_perm_data_string)

        log2_outfile = open("{0}{1}_Permuted_Log2_GMeans.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        log2_outfile.write(log2_out_string)
        log2_outfile.close()

        self.log.info("Permutation Analysis Complete.")

    def kolmogorov_smirnov2(self):
        """
        Do a Kolmogorov_Smirnov test on the target sets for each library excluding the no index.  Done on the difference
        with control library.  Writes a file for each library continuing the Log2 delta value and the p-value for each
        gene.

        :return:
        """

        self.log.info("\033[93mBegin Kolmogorov-Smirnov2 analysis.\033[m")

        for sample_name in self.tc_norm_dict:
            if sample_name == "Unknown":
                continue

            # Read TC_norm data for working library into dictionary.
            tc_control_list = []

            for tc_control_value in self.tc_norm_dict[sample_name][self.args.Species]:
                tc_control_list.append(tc_control_value-self.gtc_norm_dict[sample_name][self.args.Species])

            out_file = open("{0}{1}_{2}_KS2_Log2_Delta_Genes.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
            out_string = "Gene\tLog2\tp-value"

            for gene in natsort.natsorted(self.tc_norm_dict[sample_name]):
                ks2_list = []
                for tc_norm_value in self.tc_norm_dict[sample_name][gene]:
                    ks2_list.append(tc_norm_value-self.gtc_norm_dict[sample_name][gene])

                try:
                    v = ks_2samp(tc_control_list, ks2_list)
                except RuntimeWarning:
                    v = [1, 1]

                gene_value = \
                    round(self.gene_data_dict[gene][sample_name]-self.gene_data_dict[gene][self.args.Control_Sample], 6)
                out_string += "\n{0}\t{1}\t{2}".format(gene, gene_value, round(v[1], 6))

            out_file.write(out_string)
            out_file.close()

        self.log.info("Kolmogorov-Smirnov2 analysis complete.")

    def kolmogorov_smirnov3(self):
        """
        Do a Kolmogorov_Smirnov test on the target sets for each library excluding the no index.  Done on the difference
        with control library.  Writes a file for each library continuing the Log2 delta value and the p-value for each
        gene.

        :return:
        """

        self.log.info("\033[93mBegin Kolmogorov-Smirnov3 analysis.\033[m")

        for sample_name in self.tc_norm_dict:
            if sample_name in ["Unknown", "Plasmid"]:
                continue

            working_library_dict = collections.defaultdict(list)
            run_date = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
            out_string = "{}\nSample:\t{}\nControl Sample:\t{}\nLibrary Control:\t{}\nAlpha:\t{}\n\nGene\tLog2\t" \
                         "Original p-value\tCorrected p-value\tReject Null Hypothesis"\
                .format(run_date, sample_name, self.args.Control_Sample, self.args.Library_Control, self.args.Alpha)

            p_value_list = []
            out_string_list = []
            null_set = []
            for target_name in self.sample_td_norm_dict[sample_name]:
                gene_name = target_name.split("_")[0]
                sample_lg2 = math.log2(self.sample_td_norm_dict[sample_name][target_name])

                try:
                    ctrl_lg2 = math.log2(self.sample_td_norm_dict[self.args.Control_Sample][target_name])
                except ValueError:
                    self.log.error("{} control sample not present.  Kolmogorov-Smirnov testing aborted."
                                   .format(self.args.Control_Sample))
                    return

                delta_value = sample_lg2-ctrl_lg2

                if gene_name == self.args.Species:
                    null_set.append(delta_value)

                working_library_dict[gene_name].append(delta_value)

            for gene in working_library_dict:
                try:
                    v = ks_2samp(null_set, working_library_dict[gene])
                except RuntimeWarning:
                    v = [1, 1]

                p_value_list.append(v[1])
                gene_value = \
                    round(self.gene_data_dict[gene][sample_name]-self.gene_data_dict[gene][self.args.Control_Sample], 6)

                out_string_list.append(["\n{0}\t{1}\t{2}".format(gene, gene_value, round(v[1], 6))])

            fdr_data = stats.fdrcorrection_twostage(p_value_list, alpha=float(self.args.Alpha), method="bky")
            for v1, corrected_p, null_rejection in zip(out_string_list, fdr_data[1], fdr_data[0]):
                out_string += "{}\t{}\t{}".format(v1[0], round(corrected_p, 6), null_rejection)

            out_file = open("{0}{1}_{2}_KS3_Log2_Delta_Genes.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
            out_file.write(out_string)
            out_file.close()

        self.log.info("Kolmogorov-Smirnov3 analysis complete.")

    def kolmogorov_smirnov(self):
        """
        Do a Kolmogorov_Smirnov test on the target sets for each library excluding the no index.  Done on the difference
        with control library.  Writes a file for each library continuing the Log2 delta value and the p-value for each
        gene.

        :return:
        """

        self.log.info("\033[93mBegin Kolmogorov-Smirnov analysis.\033[m")

        for sample_name in self.sample_td_norm_dict:
            if sample_name in ["Unknown", "Plasmid", self.args.Control_Sample]:
                continue

            # Read TD_norm data for working library into dictionary.
            working_library_dict = collections.defaultdict(list)
            control_list = []

            for target_name in self.sample_td_norm_dict[sample_name]:
                gene_name = target_name.split("_")[0]
                log2_td_norm = math.log2(self.sample_td_norm_dict[sample_name][target_name])

                try:
                    log2_control_value = math.log2(self.sample_td_norm_dict[self.args.Control_Sample][target_name])
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

    def gene_group(self):
        """
        Collapse set of TDnorm target values for each gene into a single Log2 value.  Uses Geometric Mean to collapse
        data set.  Writes a single file containing all the data and creates a dictionary of the date for use later.

        :return:
        """

        delta_out_string = "Gene"
        log_out_string = "Gene"
        gene_data_dict = collections.defaultdict(lambda: collections.defaultdict(float))
        delta_gene_data_dict = collections.defaultdict(lambda: collections.defaultdict(float))

        for sample_name in natsort.natsorted(self.sample_td_norm_dict):
            if sample_name in ["Unknown", "Plasmid"]:
                continue

            tmp_delta_dict = collections.defaultdict(list)
            tmp_dict = collections.defaultdict(list)

            for target_name in self.sample_td_norm_dict[sample_name]:
                gene_name = target_name.split("_")[0]
                td_norm = self.sample_td_norm_dict[sample_name][target_name]

                tmp_delta_dict[gene_name]\
                    .append((td_norm + 1.0e-10) / (self.sample_td_norm_dict[self.args.Control_Sample][target_name] + 1.0e-10))
                tmp_dict[gene_name].append(td_norm + 1.0e-10)

            for gene in natsort.natsorted(tmp_dict):
                gene_value = math.log2(gmean(tmp_dict[gene]) - 1.0e-10)
                gene_data_dict[gene][sample_name] = gene_value

                delta_gene_value = math.log2(gmean(tmp_delta_dict[gene]))
                delta_gene_data_dict[gene][sample_name] = delta_gene_value

            delta_out_string += "\t{0}".format(sample_name)
            if sample_name != self.args.Control_Sample:
                log_out_string += "\t{0}".format(sample_name)

        for gene in natsort.natsorted(gene_data_dict):
            delta_out_string += "\n{0}".format(gene)
            log_out_string += "\n{0}".format(gene)

            for sample_name in natsort.natsorted(gene_data_dict[gene]):
                log_out_string += "\t{0}".format(gene_data_dict[gene][sample_name])

                delta_out_string += "\t{0}".format(delta_gene_data_dict[gene][sample_name])

        delta_out_file = \
            open("{0}{1}_Log2_Delta_{2}_Genes.txt".format(self.args.Working_Folder, self.args.Job_Name,
                                                          self.args.Control_Sample), "w")
        delta_out_file.write(delta_out_string)
        delta_out_file.close()

        log_out_file = open("{0}{1}_Log2_Genes.txt".format(self.args.Working_Folder, self.args.Job_Name), "w")
        log_out_file.write(log_out_string)
        log_out_file.close()

        self.gene_data_dict = gene_data_dict

    def td_norm(self):
        """
        Processes the data in the gtc_norm_dict to produce the td_norm data.  Writes the data to a file.
        :return:
        """

        for sample_name in self.gtc_norm_dict:
            if sample_name == "Unknown":
                continue

            self.log.info("TDnorm for {0}". format(sample_name))

            out_string = "Gene\tTarget\tgTC_norm for {0} mismatches\tTD_norm for {0} mismatches\tLog2_TD_norm"\
                .format(self.args.Target_Mismatch)

            for target_name in natsort.natsorted(self.gtc_norm_dict[sample_name]):
                target_key = self.target_dict[target_name]
                sample_gtc_norm = self.gtc_norm_dict[sample_name][target_name]
                library_control_gtc_norm = self.gtc_norm_dict[self.args.Library_Control][target_name]

                out_string += "\n{0}\t{1}\t{2}".format(target_name, target_key, sample_gtc_norm)

                try:
                    td_norm = sample_gtc_norm/library_control_gtc_norm
                except ZeroDivisionError:
                    td_norm = 1

                gene_name = target_name.split("_")[0]

                if gene_name == self.args.Species:
                    self.control_td_norm_dict[sample_name][target_name] = td_norm

                out_string += "\t{0}\t{1}".format(td_norm, math.log2(td_norm))
                self.sample_td_norm_dict[sample_name][target_name] = td_norm

            out_file = open("{0}{1}_{2}_TD_norm.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, sample_name), "w")
            out_file.write(out_string)
            out_file.close()

    def tc_norm(self):
        """
        This does the calculations to normalize the raw CRISPR sgRNA counts to the total counts for the library.
        """
        library_index_target_counts = collections.defaultdict(int)
        sample_tc_data = collections.defaultdict(float)
        library_tc_norm_values = collections.defaultdict(list)
        bad_targets_dict = collections.defaultdict(list)
        tmp_bad_targets_dict = collections.defaultdict(int)
        bad_targets_list = []
        percentile_list = []
        expected_diversity_list = []
        # Check library control for drop outs.
        for library_index in self.sample_mapping_dict[self.args.Library_Control]:
            try:
                tmp_data_file = open("{0}{1}_{2}_target_counts.txt"
                                     .format(self.args.Working_Folder, self.args.Job_Name, library_index))
            except FileNotFoundError:
                self.log.error("{0}_{1}_target_counts.txt not found".format(self.args.Job_Name, library_index))
                continue
            first_line = True

            # Go through each target in the Library Control target counts file.
            for line in tmp_data_file:
                if first_line:
                    first_line = False
                    continue

                line_list = [x for x in line.strip("\n").split("\t")]
                target_name = line_list[0]

                # Count the number of reads for each target in each file for each mismatch.
                for i in range(int(self.args.Target_Mismatch) + 1):
                    tmp_bad_targets_dict[target_name] += int(line_list[i + 2])

            tmp_data_file.close()
            total_sublibrary_control_counts = sum(tmp_bad_targets_dict.values())
            total_guides = len(tmp_bad_targets_dict)
            expected_diversity_list.append(0.05*((total_sublibrary_control_counts/total_guides)/total_sublibrary_control_counts))
            for target_name, target_count in tmp_bad_targets_dict.items():
                bad_targets_dict[target_name].append((target_count/total_sublibrary_control_counts)+1.0e-10)
                percentile_list.append(target_count/total_sublibrary_control_counts)

            tmp_bad_targets_dict.clear()

        # Check for missing Library Control Targets.
        bad_target_cutoff = numpy.percentile(numpy.array(percentile_list), 2.5)
        # bad_target_cutoff = gmean(expected_diversity_list)
        percentile_list.clear()

        bad_target_outstring = "Targets excluded from analysis.\nFile Generated {}\nLibrary Control File: {}\n" \
                               "2.5 Percentile gTCnorm Cutoff: {}\n\nsgRNA Name\tgTCnorm\n"\
            .format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.args.Library_Control,
                    format(bad_target_cutoff, '.4g'))

        for target_name in bad_targets_dict:
            if target_name not in self.target_dict:
                self.log.error("{} not found in {}.  Confirm correct Target File is being used."
                               .format(target_name, self.args.Target_File))

            if gmean(bad_targets_dict[target_name]) <= bad_target_cutoff or target_name == "53BP1-intron":
                bad_targets_list.append(target_name)
                bad_value = format(gmean(bad_targets_dict[target_name])-1.0e-10, '.4g')

                # Because of a Python/CPU math quirk, this will never be 0.
                if float(bad_value) <= 1.0e-9:
                    bad_value = 0

                bad_target_outstring += "{}\t{}\n".format(target_name, bad_value)

                self.log.warning("Masking {} with a gTCnorm of {} in control {}"
                                 .format(target_name, bad_value, self.args.Library_Control))

        bad_guide_outfile = open("{}{}_Masked_Targets.txt"
                                 .format(self.args.Working_Folder, self.args.Job_Name), 'w')

        bad_guide_outfile.write(bad_target_outstring)
        bad_guide_outfile.close()

        # Process each library.
        for sample_name, library_index_list in self.sample_mapping_dict.items():
            if sample_name == "Unknown":
                continue

            for library_index in library_index_list:
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

                    target_name = line_list[0]
                    target_data_list.append(line_list)

                    # Count the number of reads for each target in each file for each mismatch.
                    for i in range(int(self.args.Target_Mismatch)+1):
                        library_index_target_counts[target_name] += int(line_list[i+2])

                    # Adding 1 to prevent errors from 0 counts.
                    library_index_target_counts[target_name] += 1
                tmp_data_file.close()

                library_index_total_count = sum(library_index_target_counts.values())
                # Normalize for reads in individual libraries
                for target_name in library_index_target_counts:
                    library_tc_norm_values[target_name]\
                        .append(library_index_target_counts[target_name]/library_index_total_count)
                library_index_target_counts.clear()

            # Determine the gTC_norm for each sample
            for target_name in natsort.natsorted(library_tc_norm_values):
                # Skip bad guides
                if target_name in bad_targets_list:
                    continue

                gene_name = target_name.split("_")[0]
                self.tc_norm_dict[sample_name][gene_name].extend(library_tc_norm_values[target_name])
                gtc_norm = statistics.mean(library_tc_norm_values[target_name])

                self.gtc_norm_dict[sample_name][target_name] = gtc_norm

                if gtc_norm == 0:
                    self.log.warning("gtc_norm=0|{} {} {}"
                                     .format(sample_name, target_name, sample_tc_data[target_name]))

            library_tc_norm_values.clear()

    def __summary_output(self, multiprocessor_tmp_data_list):
        """
        Process temporary data file into final output.
        """

        self.log.debug("Begin Writing Data Output File")
        args = self.args

        index_header = "Index Name\tSample Name\tSample Replica\tTotal Indexed\tFull Indexed Reads\tShort Indexed Reads"

        for i in range(int(args.Index_Mismatch)+1):
            index_header += "\t{}_mismatches".format(i)
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
                index_seq = data[1]
            except TypeError:
                continue

            index_name = self.index_dict[index_seq][4]
            # Index Name, Sample Name, Sample Replica
            param_string += "{}\t{}\t{}\t{}\t{}\t{}"\
                .format(index_name, self.index_dict[index_seq][5], self.index_dict[index_seq][6],
                        self.index_dict[index_seq][1], self.index_dict[index_seq][2], self.index_dict[index_seq][3])

            # Mismatched counts
            for i in range(len(self.index_dict[index_seq][0])):
                param_string += "\t{0}".format(self.index_dict[index_seq][0][i])

            # Targeted, Not Targeted, Fraction Targeted.
            param_string += "\t{}\t{}\t{}\n"\
                .format(data[0], int(self.index_dict[index_seq][1]) - int(data[0]),
                        int(data[0])/int(self.index_dict[index_seq][1]))

        run_param_out.write(param_string)
        run_param_out.close()

        self.log.debug("Data Summary File Written.")

    def dictionary_builds(self):
        """
        Build dictionaries, initialize output files and capture list of output file names.
        :return:
        """

        # Read master index file into a dictionary.
        fastq_out_list = []
        master_index_dict = {}
        fastq_file_dict = collections.defaultdict(object)
        index_dict = {"Unknown": [[0] * (int(self.args.Index_Mismatch) + 1), 0, 0, 0, "Unknown", "Unknown", "Unknown"]}

        with open(self.args.Master_Index_File) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                master_index_dict[l_list[0]] = l_list[1]

        for sample in self.index_list:
            index_sequence = master_index_dict[sample[0]]
            if index_sequence in index_dict:
                self.log.error("The index {0} is duplicated.  Correct the error in {1} and try again."
                               .format(sample[0], self.args.Index_File))
                raise SystemExit(1)

            # for each sample name append a list of all index ID's
            self.sample_mapping_dict[sample[1]].append(sample[0])
            index_dict[index_sequence] = \
                [[0]*(int(self.args.Index_Mismatch)+1), 0, 0, 0, sample[0], sample[1], sample[2]]

            if self.args.Target_Search:
                fastq_file_dict[index_sequence] = \
                    FASTQ_Tools.Writer(self.log, "{0}{1}_{2}.fastq"
                                       .format(self.args.Working_Folder, self.args.Job_Name, index_sequence))

                fastq_out_list.append("{0}{1}_{2}.fastq"
                                      .format(self.args.Working_Folder, self.args.Job_Name, index_sequence))

        # This is for no index found.
        self.index_list.append(("Unknown", "Unknown", "Unknown"))

        if self.args.Analyze_Unknowns:
            master_index_dict["Unknown"] = "Unknown"
            fastq_file_dict["Unknown"] = \
                FASTQ_Tools.Writer(self.log, "{0}{1}_Unknown.fastq"
                                   .format(self.args.Working_Folder, self.args.Job_Name))

            fastq_out_list.append("{0}{1}_Unknown.fastq".format(self.args.Working_Folder, self.args.Job_Name))

        # Fill target list and dictionary.  Do initial quality check on target file for duplicates.
        target_list = []

        for target in Tool_Box.FileParser.indices(self.log, self.args.Target_File):
            try:
                target_seq = target[2][int(self.args.Target_Start):][:int(self.args.Target_Length)]
            except ValueError:
                target_seq = target[2]

            target_name = target[1]
            target_list.append((target_seq, target_name))
            self.targets[len(target_seq)].append(target_seq)

            if target_seq in self.target_dict:
                self.log.error("The target sequence in {0} is duplicated.  Correct the error in {1} and try again."
                               .format(target, self.args.Target_File))
                raise SystemExit(1)
            elif target_seq in self.target_dict[target_name]:
                self.log.error("Target name in {0} is duplicated.  Correct the error in {1} and try again."
                               .format(target, self.args.Target_File))
                raise SystemExit(1)

            self.target_dict[target_name] = target_seq

        similarity_count = 0
        off_target_count = 0
        # Do a fine scale quality analysis of targets looking for similar sequences and off target guides.
        for target_name in self.target_dict:
            for target in target_list:
                mismatch_index = Sequence_Magic.match_maker(target[0], self.target_dict[target_name])
                if 0 < mismatch_index <= 3:
                    target_gene_name = target[1].split("_")[0]  # self.target_dict[target[1]].split("_")[0]
                    query_gene_name = target_name.split("_")[0]  # self.target_dict[target_name].split("_")[0]

                    if target_gene_name != query_gene_name:
                        off_target_count += 1
                        self.log.debug("!!!POTENTIAL OFF TARGET!!! {0} differs from {1} by {2}"
                                       .format(target_name, target[1], mismatch_index))
                    else:
                        similarity_count += 1
                        self.log.debug("{} differs from {} by only {}"
                                       .format(target_name, target[1], mismatch_index))
        if similarity_count > 0:
            self.log.info("{0} targets similar to each other".format(similarity_count))
            self.log.info("{0} potential off target guides".format(off_target_count))

        self.log.info("Dictionaries Built")
        return index_dict, fastq_file_dict, fastq_out_list, master_index_dict

    def fastq_processor(self):
        """
        Demultiplex FASTQ file.  Writes new files every 1.5 million reads processed.
        :rtype: object
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
                # this block limits reads for debugging.
                Tool_Box.debug_messenger("Limiting reads here to 1.5 million")
                self.log.warning("Limiting reads here to 1.5 million")
                eof = True
            index_loop_count = 0
            query_sequence = ""

            if args.Platform == "Illumina":
                # The indices are after the last ":" in the header.
                query_sequence = fastq_read.name.split(":")[-1]

            elif not args.Platform == "Illumina" and not args.Platform == "Ion":
                Tool_Box.debug_messenger("--Platform must be Illumina or Ion")
                raise SystemExit(1)

            while not index_found:
                for index_seq in self.index_dict:
                    index_loop_count += 1

                    # No need to search for the "Unknown" key.
                    unknown = True
                    mismatch_index = int(args.Index_Mismatch)+1
                    index_key_length = len(index_seq)
                    if index_seq != "Unknown":
                        unknown = False

                        if args.Platform == "Ion":
                            query_sequence = fastq_read.seq[:index_key_length]

                        mismatch_index = Sequence_Magic.match_maker(index_seq, query_sequence)

                    if mismatch_index <= int(args.Index_Mismatch) and not unknown:
                        index_found = True
                        self.index_dict[index_seq][1] += 1  # Total reads with index.

                        # Check read length and only output if the reads are long enough.
                        if len(fastq_read.seq) > int(args.Min_Length):
                            # if args.Platform == "Ion":
                            #     FASTQ_Tools.read_trim(fastq_read, trim5=index_key_length)

                            fastq_read.name = "{0}|{1}".format(index_seq, fastq_read.name)
                            Read = collections.namedtuple('Read', 'name, seq, index, qual')
                            new_read = Read(fastq_read.name, fastq_read.seq, fastq_read.index, fastq_read.qual)
                            output_dict[index_seq].append(new_read)
                            self.index_dict[index_seq][2] += 1  # Full reads with index.
                            self.index_dict[index_seq][0][mismatch_index] += 1  # Mismatch count.
                            break
                        else:
                            self.index_dict[index_seq][3] += 1  # Short reads with index.
                        # Break the loop if nothing matches.
                    elif index_loop_count == index_count:
                        index_found = True
                        self.fastq_read_counts[1] += 1  # Total Unknown index reads in FASTQ file.
                        self.index_dict["Unknown"][1] += 1  # Total reads without index.

                        if len(fastq_read.seq) > int(args.Min_Length):
                            self.fastq_read_counts[3] += 1  # Full reads with no index.
                            self.index_dict["Unknown"][2] += 1  # Full reads with no index.
                            if self.args.Analyze_Unknowns:
                                head_seq = fastq_read.seq[:index_key_length]
                                fastq_read.name = "Unknown|{0}|{1}".format(head_seq, fastq_read.name)
                                FASTQ_Tools.read_trim(fastq_read, trim5=index_key_length)
                                Read = collections.namedtuple('Read', 'name, seq, index, qual')
                                new_read = Read(fastq_read.name, fastq_read.seq, fastq_read.index, fastq_read.qual)
                                output_dict["Unknown"].append(new_read)
                        else:
                            self.fastq_read_counts[2] += 1  # Short reads with no index.
                            self.index_dict["Unknown"][3] += 1  # Short reads with no index.

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
            return []

        # If the FASTQ file is empty remove it and get out of here.
        elif os.stat(fq_file).st_size < 50:
            log.warning("{0} is empty; File removed." .format(fq_file))
            os.remove(str(fq_file))
            return []

        t0 = clock()

        # Retrieve the index sequence and index name for the file we are processing.
        index_key = ntpath.basename(fq_file).split(".")[0].split("_")[-1]
        index_name = index_dict[index_key][4]
        multiprocessor_tmp_data_list = [0, index_key]
        index_key_length = len(index_key)

        if not args.Analyze_Unknowns and index_key == "Unknown":
            log.info("\033[1;31mNotice:\033[m  Skipping {0} at user request.".format(ntpath.basename(fq_file)))
            return multiprocessor_tmp_data_list

        target_data_outstring = "Target\tTarget_Key"
        for i in range(int(args.Target_Mismatch)+1):
            target_data_outstring += "\t{0}_mismatches".format(i)

        # Iterate the FASTQ file; extract the target region; reverse-compliment it; check it against the target list
        # for matches; tabulate results.
        target_found_pos_list = []
        no_target_pos_list = []
        fastq_read_count = 0
        target_count = 0
        anchor_dict = {"index_key": index_key, "total_target_pos_list": [], "no_anchor_count": 0, "anchor_count": 0}
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
            # If the read length is too short skip it and go onto the next one.
            if len(fastq_read.seq) <= int(args.Min_Length):
                continue

            # Find the first position of the sgRNA
            unknown_seq_start, anchor_dict = SyntheticLethal.__anchor_search(args, fastq_read, anchor_dict)

            # Compare the sgRNA sequence to the targets.
            target_seq, mismatch_index = SyntheticLethal.__target_match(targets_dict, fastq_read,
                                                                        unknown_seq_start, args)
            # Count our targets or no targets.
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
        # Total Anchors data
        freq_pos_outstring = "Total_Anchors\tPosition\tCount\tFrequency"
        freq_pos_outstring = \
            SyntheticLethal.__frequency_outstring(freq_pos_outstring, anchor_dict["total_target_pos_list"],
                                                  index_key_length, anchor_dict)
        # Total Targets Data
        freq_pos_outstring += "\nTargets_Found\tPosition\tCount\tFrequency"
        freq_pos_outstring = \
            SyntheticLethal.__frequency_outstring(freq_pos_outstring, target_found_pos_list, index_key_length,
                                                  anchor_dict)
        # No Target Data
        freq_pos_outstring += "\nNo_Targets_Found\tPosition\tCount\tFrequency"
        freq_pos_outstring = \
            SyntheticLethal.__frequency_outstring(freq_pos_outstring, no_target_pos_list, index_key_length,
                                                  anchor_dict)

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
    def __target_match(targets_dict, fastq_read, unknown_seq_start, args):
        """
        Does a Levenshtein search for items in target list.
        :return:
        :param targets_dict:
        :param fastq_read:
        :param unknown_seq_start:
        :param args:
        """

        target_mismatch = int(args.Target_Mismatch)
        targets_found_dict = collections.defaultdict(list)

        # Go through targets based on size.
        for target_length in targets_dict:

            if args.RevComp:
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

    @staticmethod
    def __anchor_search(args, fastq_read, anchor_dict):
        """
        Looks for anchor sequence and returns the start position of the sgRNA.
        :param args:
        :param fastq_read:
        :param anchor_dict:
        :return:
        """
        anchor_found = False
        start_pos = int(args.AnchorStart) - len(anchor_dict["index_key"])
        while not anchor_found:
            mismatch_index = Sequence_Magic.match_maker(
                args.AnchorSeq, fastq_read.seq[start_pos:][:len(args.AnchorSeq)])

            # If we do not find the anchor sequence exit the loop and go to the next read.
            if start_pos > int(args.AnchorStop):
                unknown_seq_start = int(args.Expected_Position) - int(args.Target_Padding)
                anchor_dict["no_anchor_count"] += 1
                break

            elif mismatch_index <= int(args.AnchorMismatch):
                anchor_found = True
                unknown_seq_start = start_pos + len(args.AnchorSeq)
                anchor_dict["anchor_count"] += 1
                anchor_dict["total_target_pos_list"].append(unknown_seq_start)

            start_pos += 1
        return unknown_seq_start, anchor_dict

    @staticmethod
    def __frequency_outstring(freq_pos_outstring, data_list, index_key_length, anchor_dict):
        """
        This processess the data for the frequency position data file.
        :param freq_pos_outstring:
        :param data_list:
        :param index_key_length:
        :param anchor_dict:
        :return:
        """
        total_target_pos_counter = collections.Counter(data_list)
        for k in total_target_pos_counter.items():
            freq_pos_outstring += \
                "\n{0}\t{1}\t{2}".format(k[0]+index_key_length, k[1], round((k[1]/anchor_dict["anchor_count"]), 4))
        return freq_pos_outstring
