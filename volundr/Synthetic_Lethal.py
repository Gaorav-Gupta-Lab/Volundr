"""
Synthetic_Lethal.py 3.0.0
    October 20, 2021
    FASTQ demultiplexing done in parallel.

Synthetic_Lethal.py 2.0.0
    August 30, 2019
    Added multiple sample p-value correction.  Added percentile output.  Added output file for masked sgRNA sequences.
    Added a library control and sample control option.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2021
"""

import csv
import gc
import ntpath
import re
import time
from time import clock
import datetime
import collections
import itertools
import os
import statistics
import gzip
import statsmodels.stats.multitest as stats
import math
import numpy
from scipy.stats import gmean, ks_2samp, norm, combine_pvalues
import natsort
import pathos
from Valkyries import FASTQ_Tools, FASTQReader, Tool_Box, Sequence_Magic


__author__ = 'Dennis A. Simpson'
__version__ = '3.5.2'
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
        self.SampleManifest = Tool_Box.FileParser.indices(log, self.args.SampleManifest)
        self.fastq_read_counts = [0, 0, 0, 0]  # Tracks number of reads in input file.
        self.log = log
        self.permuted_null_dict = None
        self.sample_data_dict, self.fastq_file_dict, self.fastq_out_list, self.master_index_dict = \
            SyntheticLethal.dictionary_builds(self)

    def fastq_analysis(self):
        """
        This will send the FASTQ file off to be demultiplexed and quantified.  When that is done this will spawn
        parallel jobs to search for target sequences in the demultiplexed FASTQ files.  Each parallel job processes a
        single FASTQ file.
        """

        self.fastq_processing()

        self.log.info("Spawning \033[96m{0}\033[m parallel job(s) to search \033[96m{1}\033[m FASTQ files for targets"
                      .format(self.args.Spawn, len(self.fastq_out_list)))

        multiprocessor_tmp_data_list = []
        p = pathos.multiprocessing.Pool(self.args.Spawn)
        multiprocessor_tmp_data_list.append(
            p.starmap(self.target_search,
                      zip(self.fastq_out_list,
                          itertools.repeat((self.args, self.targets, self.log, self.sample_data_dict)))))

        self.log.info(" ***All Parallel Jobs Complete.***")

        # Process summary data
        self.__summary_output(multiprocessor_tmp_data_list)

        if self.args.Delete_Demultiplexed_FASTQ:
            self.log.debug("Deleting modified FASTQ files from system.")
            Tool_Box.delete(self.fastq_out_list)

    def statistics(self):
        """
        Runs through the methods to do the analysis in the correct order.
        """
        self.log.info("\033[93mBegin TCnorm calculations\033[m")
        self.tc_norm()
        self.log.info("TCnorm calculations complete.")
        return

        self.log.info("\033[93mBegin TDnorm calculations and Log2 transformation\033[m")
        self.td_norm()
        self.log.info("TDnorm calculations and Log2 transformation complete.")

        self.control_permutation()

        self.log.info("\033[93mBegin collapsing gene target groups into individual genes\033[m")
        self.gene_group()
        self.log.info("Gene data manipulations complete.")

        self.kolmogorov_smirnov()

    def control_permutation(self):
        """
        This bad boy does a permutation analysis using the control targets.  Outputs several data files.
        :param self:
        :return:
        """

        self.log.info("\033[93mBegin Control Permutation Analysis.\033[m")
        self.permuted_null_dict = collections.defaultdict(list)
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
        if self.args.Write_TDnorm_Log2_sgRNA_Control_File:
            for control_target_name in selection_data_dict:
                raw_data_string += "\n{0}".format(control_target_name)
                for v in selection_data_dict[control_target_name]:
                    raw_data_string += "\t{0}".format(math.log2(float(v)))

            # Write Log2 control target data to a file by itself.
            raw_data_out_string = "{0}{1}".format(raw_data_header, raw_data_string)
            raw_data_outfile = \
                open("{0}{1}_TDnorm_Log2_Control_Targets.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
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
                            if td_norm_control_ratio > 0:
                                working_dict[group_key][sample_name].append(td_norm_control_ratio)

                        except ZeroDivisionError:
                            pass

                            # self.log.error("Cannot complete Permutation Analysis.  {} {} missing from indices."
                            #                .format(self.args.Control_Sample, control_target_key))
                            # return

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
                str(numpy.percentile(numpy.array(percentile_dict[sample_name]), float(self.args.LowerPercentile),
                                     interpolation='linear'))
            upper_limit = \
                str(numpy.percentile(numpy.array(percentile_dict[sample_name]), float(self.args.UpperPercentile),
                                     interpolation='linear'))
            sample_name_list.append(sample_name)
            lower_limit_list.append(lower_limit)
            upper_limit_list.append(upper_limit)
            self.permuted_null_dict[sample_name] = [upper_limit, lower_limit]

        if self.args.Write_Permuted_Log2_Data_File:
            log2_out_string += "Sample:\t{}\n".format("\t".join(sample_name_list))
            log2_out_string += "Upper Limit:\t{}\n".format("\t".join(upper_limit_list))
            log2_out_string += "Lower Limit:\t{}\n\n{}".format("\t".join(lower_limit_list), log2_perm_data_string)

            log2_outfile = \
                open("{0}{1}_Permuted_Log2_GMeans.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
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

        for sample_name in self.tc_norm_dict:
            if sample_name in ["Unknown", self.args.Library_Control]:
                continue

            working_library_dict = collections.defaultdict(list)

            lower_limit = round(float(self.permuted_null_dict[sample_name][0]), 3)
            upper_limit = round(float(self.permuted_null_dict[sample_name][1]), 3)

            run_date = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
            out_string = "{}\nSample:\t{}\nControl Sample:\t{}\nLibrary Control:\t{}\nAlpha:\t{}\n" \
                         "Upper Null Set Limit:\t{}\nLower Null Set Limit:\t{}\n\nGene\tLog2\t" \
                         "Original p-value\t-Log10(pVal)\tCorrected p-value\tReject Null Hypothesis"\
                .format(run_date, sample_name, self.args.Control_Sample, self.args.Library_Control, self.args.Alpha,
                        lower_limit, upper_limit)

            p_value_list = []
            out_string_list = []
            null_set = []
            for target_name in self.sample_td_norm_dict[sample_name]:
                gene_name = target_name.split("_")[0]

                sample_lg2 = math.log2(self.sample_td_norm_dict[sample_name][target_name])
                try:
                    ctrl_lg2 = math.log2(self.sample_td_norm_dict[self.args.Control_Sample][target_name])

                except ValueError:
                    continue

                delta_value = sample_lg2 - ctrl_lg2
                working_library_dict[gene_name].append(delta_value)
                if gene_name == self.args.Species:
                    null_set.append(delta_value)

            for gene in working_library_dict:
                try:
                    v = ks_2samp(null_set, working_library_dict[gene])
                except RuntimeWarning:
                    v = [1, 1]

                p_value_list.append(v[1])
                gene_value = \
                    round(self.gene_data_dict[gene][sample_name]-self.gene_data_dict[gene][self.args.Control_Sample], 4)
                neg_log = round(-1*math.log(v[1], 10), 4)
                p_val = round(v[1], 4)
                out_string_list.append(["\n{}\t{}\t{}\t{}".format(gene, gene_value, p_val, neg_log)])

            fdr_data = stats.fdrcorrection_twostage(p_value_list, alpha=float(self.args.Alpha), method="bky")
            for v1, corrected_p, null_rejection in zip(out_string_list, fdr_data[1], fdr_data[0]):
                out_string += "{}\t{}\t{}".format(v1[0], round(corrected_p, 4), null_rejection)

            out_file = open("{0}{1}_{2}_KS_Log2_Delta_Genes.txt"
                            .format(self.args.WorkingFolder, self.args.Job_Name, sample_name), "w")
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
                    .append((td_norm + 1.0e-10) / (self.sample_td_norm_dict[self.args.Control_Sample][target_name]
                                                   + 1.0e-10))
                
                tmp_dict[gene_name].append(td_norm + 1.0e-10)

            for gene in natsort.natsorted(tmp_dict):
                # Sometimes all the guides are 0.
                try:
                    gene_value = math.log2(gmean(tmp_dict[gene]) - 1.0e-10)
                except ValueError:
                    gene_value = math.log2(1.0e-10)

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

        if self.args.Verbose == "DEBUG":
            delta_out_file = \
                open("{0}{1}_Log2_Delta_{2}_Genes.txt".format(self.args.WorkingFolder, self.args.Job_Name,
                                                          self.args.Control_Sample), "w")
            delta_out_file.write(delta_out_string)
            delta_out_file.close()

        if self.args.Write_Log2_sgRNA_File:
            log_out_file = open("{0}{1}_Log2_Genes.txt".format(self.args.WorkingFolder, self.args.Job_Name), "w")
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

            if self.args.Write_TDnorm_Log2_sgRNA_Sample_File:
                out_file = open("{0}{1}_{2}_TD_norm.txt"
                                .format(self.args.WorkingFolder, self.args.Job_Name, sample_name), "w")
                out_file.write(out_string)
                out_file.close()

    def bad_correlation(self, library_tc_norm_values):
        """

        :rtype: object
        """
        tmp_percentile_data_list = []
        sample_bad_targets_dict = collections.defaultdict(list)

        for target_name in library_tc_norm_values:
            avg_guide_tc_norm = statistics.mean(library_tc_norm_values[target_name])
            tmp_guide_log2 = []
            for guide_tc_norm in library_tc_norm_values[target_name]:
                try:
                    tmp_guide_log2.append(math.log2(guide_tc_norm / avg_guide_tc_norm))
                    sample_bad_targets_dict[target_name].append(math.log2(guide_tc_norm / avg_guide_tc_norm))
                    tmp_percentile_data_list.append(math.log2(guide_tc_norm / avg_guide_tc_norm))
                except ValueError:
                    pass
                    # tmp_guide_log2.append(1.0)
                    # sample_bad_targets_dict[target_name].append(1.0)
            # tmp_percentile_data_list.extend(tmp_guide_log2)

        upper_limit = \
            (numpy.percentile(numpy.array(tmp_percentile_data_list), self.args.UpperGuideLimit,
                              interpolation='linear'))
        lower_limit = \
            (numpy.percentile(numpy.array(tmp_percentile_data_list), self.args.LowerGuideLimit,
                              interpolation='linear'))

        return sample_bad_targets_dict, upper_limit, lower_limit

    def file_read(self, sample, control_file=True):
        tmp_bad_targets_dict = collections.defaultdict(int)
        tmp_tc_norm_dict = collections.defaultdict(list)
        bad_targets_dict = collections.defaultdict(list)
        percentile_list = []
        sample_control_dict = collections.defaultdict(list)

        for library_index in self.sample_mapping_dict[sample]:
            re.sub('[\s]', "", library_index)

            try:
                tmp_data_file = open("{0}{1}_{2}_target_counts.txt"
                                     .format(self.args.DataFiles, self.args.Job_Name, library_index))
            except FileNotFoundError:
                self.log.error("{0}_{1}_target_counts.txt not found".format(self.args.Job_Name, library_index))
                raise SystemExit(1)

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

            control_key = "{}_{} controls".format(sample, library_index)
            guide_key = "{}_{} guides".format(sample, library_index)
            count_key = "{}_{} counts".format(sample, library_index)
            if not control_file:
                sample_control_dict[count_key].append(total_sublibrary_control_counts)
            else:
                sample_control_dict["Library Reads"].append(total_sublibrary_control_counts)

            for target_name, target_count in tmp_bad_targets_dict.items():
                gene_name = target_name.split("_")[0]

                if not control_file and gene_name == self.args.Species:
                    sample_control_dict[control_key]\
                        .append([target_name, target_count/total_sublibrary_control_counts,
                                 total_sublibrary_control_counts])
                elif not control_file:
                    sample_control_dict[guide_key]\
                        .append([target_name, target_count/total_sublibrary_control_counts,
                                 total_sublibrary_control_counts])
                
                if target_count > 0:
                    tmp_tc_norm_dict[target_name].append(target_count/total_sublibrary_control_counts)

                bad_targets_dict[target_name].append((target_count/total_sublibrary_control_counts)+1.0e-10)
                percentile_list.append(target_count/total_sublibrary_control_counts)

            tmp_bad_targets_dict.clear()
            tmp_data_file.close()

        total_sublibrary_control_counts = sum(tmp_bad_targets_dict.values())

        for target_name, target_count in tmp_bad_targets_dict.items():
            if target_count > 0:
                tmp_tc_norm_dict[target_name].append(target_count / total_sublibrary_control_counts)
            bad_targets_dict[target_name].append((target_count / total_sublibrary_control_counts))
            percentile_list.append(target_count / total_sublibrary_control_counts)

        tmp_bad_targets_dict.clear()
        # Check for missing or vastly under represented Library Control Targets.
        percentile1 = self.args.Bad_sgRNA_Lower_Percentile
        percentile2 = self.args.Bad_sgRNA_Upper_Percentile

        upper_limit = \
            (numpy.percentile(numpy.array(percentile_list), percentile2, interpolation='linear'))
        lower_limit = \
            (numpy.percentile(numpy.array(percentile_list), percentile1, interpolation='linear'))

        return upper_limit, lower_limit, tmp_tc_norm_dict, bad_targets_dict, sample_control_dict

    def tc_norm(self):
        """
        This does the calculations to normalize the raw CRISPR sgRNA counts to the total counts for the library.
        :rtype: object
        """
        library_index_target_counts = collections.defaultdict(int)
        sample_tc_data = collections.defaultdict(float)

        sample_control_dict = collections.defaultdict(float)
        bad_targets_list = []

        # Process library control
        library_upper_limit, library_lower_limit, tmp_tc_norm_dict, bad_targets_dict, sample_control_data = \
            self.file_read(self.args.Library_Control)

        bad_target_outstring = "sgRNA Targets excluded from analysis.\nFile Generated {}\nLibrary Control File: {}\n" \
                               "{} Lower Percentile gTCnorm Lower Cutoff: {}\n" \
                               "{} Lower Percentile gTCnorm Upper Cutoff: {}\n\nsgRNA Name\tgTCnorm\n"\
            .format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.args.Library_Control,
                    self.args.Bad_sgRNA_Lower_Percentile, format(library_lower_limit, '.4g'),
                    self.args.Bad_sgRNA_Upper_Percentile, format(library_upper_limit, '.4g'))

        # sample_bad_targets_dict, upper_limit, lower_limit = self.bad_correlation(tmp_tc_norm_dict)
        # bad_correlation_outstring = "\n"
        for target_name in bad_targets_dict:
            if self.args.TargetSearch and target_name not in self.target_dict:
                self.log.error("{} not found in {}.  Confirm correct Target File is being used."
                               .format(target_name, self.args.Target_File))
            target_value = gmean(bad_targets_dict[target_name])
            if target_value <= library_lower_limit or target_value >= library_upper_limit:
                bad_targets_list.append(target_name)
                bad_value = format(target_value, '.4g')

                # Because of a Python/CPU math quirk, this will never be 0.
                if float(bad_value) <= 1.0e-9:
                    bad_value = 0

                bad_target_outstring += "{}\t{}\n".format(target_name, bad_value)

                self.log.warning("Masking {} with a gTCnorm of {} in control {}"
                                 .format(target_name, bad_value, self.args.Library_Control))

        bad_targets_list = list(set(bad_targets_list))

        bad_guide_outfile = open("{}{}_{}_Masked_Targets.txt"
                                 .format(self.args.WorkingFolder, self.args.Job_Name, self.args.Library_Control), 'w')
        bad_guide_outfile.write(bad_target_outstring)
        bad_guide_outfile.close()

        # Process sample control
        upper_limit, lower_limit, tmp_tc_norm_dict, targets_dict, control_dict = \
            self.file_read(self.args.Control_Sample)

        for target_name in targets_dict:

            if target_name in bad_targets_list:
                continue

            try:
                sample_control_dict[target_name] = gmean(targets_dict[target_name])
            except ValueError:
                sample_count = len(targets_dict[target_name])
                tmp_list = []
                for i in range(sample_count):
                    target_value = targets_dict[target_name][i]
                    read_count = control_dict["Library Reads"][i]
                    if target_value == 0:
                        target_value = 1/read_count
                    tmp_list.append(target_value)
                sample_control_dict[target_name] = gmean(tmp_list)

        # Process each library.
        for sample_name, library_index_list in self.sample_mapping_dict.items():

            if sample_name == "Unknown" or sample_name == self.args.Control_Sample:
                continue

            sample_control_guide_dict = collections.defaultdict(list)
            upper_limit, lower_limit, tmp_tc_norm_dict, sample_targets_dict, sample_dict = \
                self.file_read(sample_name, control_file=False)

            single_guide_dict = collections.defaultdict(lambda: collections.defaultdict(list))
            sample_data_dict = collections.defaultdict(lambda: collections.defaultdict(float))
            guide_counts_dict = collections.defaultdict(list)
            guide_counts_dict2 = collections.defaultdict(list)
            gene_pval_dict = collections.defaultdict(list)
            pval_used = []
            log_delta_control_guide_list = []
            index_key_list = []

            for library_index in library_index_list:
                index_key_list.append(library_index)
                control_key = "{}_{} controls".format(sample_name, library_index)
                guide_key = "{}_{} guides".format(sample_name, library_index)
                count_key = "{}_{} counts".format(sample_name, library_index)
                library_read_count = sample_dict[count_key][0]

                # Process the control sgRNA data
                control_guide_delta_list = []
                for control_target_data in sample_dict[control_key]:
                    target_name = control_target_data[0]
                    target_value = control_target_data[1]

                    if target_value == 0:
                        target_value = 1/library_read_count

                    if target_name not in bad_targets_list:
                        gene_name = target_name.split("_")[0]
                        sample_control_guide_dict[target_name].append(target_value)

                        delta_val = target_value/sample_control_dict[target_name]
                        log_delta_val = math.log2(delta_val)
                        single_guide_dict[gene_name][library_index].append([target_name, log_delta_val, delta_val])
                        log_delta_control_guide_list.append(log_delta_val)
                        control_guide_delta_list.append(delta_val)
                        sample_data_dict[library_index][target_name] = log_delta_val

                # Process the sample sgRNA data
                for guide_data in sample_dict[guide_key]:
                    target_name = guide_data[0]
                    target_value = guide_data[1]

                    if not target_value > 0:
                        target_value = 1/library_read_count
                        # target_value = library_lower_limit

                    if target_name not in bad_targets_list:
                        gene_name = target_name.split("_")[0]
                        delta_val = target_value/sample_control_dict[target_name]
                        log_delta_val = math.log2(delta_val)
                        sample_data_dict[target_name] = log_delta_val
                        single_guide_dict[gene_name][library_index].append([target_name, log_delta_val, delta_val])
                        guide_counts_dict2[target_name].append(target_value)

                sample_data_dict[library_index]['percentile_upperlimit'] = \
                    (numpy.percentile(numpy.array(log_delta_control_guide_list), self.args.UpperGuideLimit,
                     interpolation='linear'))

                sample_data_dict[library_index]['percentile_lowerlimit'] = \
                    (numpy.percentile(numpy.array(log_delta_control_guide_list), self.args.LowerGuideLimit,
                     interpolation='linear'))

            for index_key in library_index_list:
                percentile_upperlimit = sample_data_dict[index_key]['percentile_upperlimit']
                percentile_lowerlimit = sample_data_dict[index_key]['percentile_lowerlimit']
                file_run = datetime.datetime.today().strftime(self.date_format)
                outdata = "Running:\t{} Synthetic_Lethal v{}\nFile Generated:\t{}\nLibrary Control:\t{}\n" \
                          "Sample Control:\t{}\nSample:\t{}_{}\nLower Limit:\t{}\nUpper Limit:\t{}\nTarget\t" \
                          "Log2\tLower pVal\tUpper pVal\n"\
                    .format(__package__, __version__, file_run, self.args.Library_Control, self.args.Control_Sample,
                            sample_name, index_key, round(percentile_lowerlimit, 3), round(percentile_upperlimit, 3))

                for gene in single_guide_dict:
                    guide_count = len(single_guide_dict[gene][index_key])

                    if gene in gene_pval_dict:
                        gene_pval_dict[gene][0] += guide_count
                    else:
                        gene_pval_dict[gene] = [guide_count, 0, 0]

                    depleted_count = 0
                    enriched_count = 0
                    for target_data in single_guide_dict[gene][index_key]:
                        # target_name = target_data[0]
                        target_value = target_data[1]
                        guide_counts_dict[gene].append(target_value)

                        if target_value <= percentile_lowerlimit:
                            depleted_count += 1

                        if target_value >= percentile_upperlimit:
                            enriched_count += 1
                    gene_pval_dict[gene][1] += depleted_count
                    gene_pval_dict[gene][2] += enriched_count
                    enriched_sig_value = (self.args.LowerGuideLimit/100)**enriched_count
                    depleted_sig_value = (self.args.LowerGuideLimit/100)**depleted_count
                    upper_pval = "undefined"
                    lower_pval = "undefined"

                    # if not gene == self.args.Species:
                    # FixMe: This is a mess from before
                    if gene == "Dragons":
                        upper_pval = \
                            enriched_sig_value*((math.factorial(guide_count)/math.factorial(enriched_count))/math.factorial(guide_count-enriched_count))
                        upper_pval = round(upper_pval, 4)

                        lower_pval = \
                            depleted_sig_value*(math.factorial(guide_count)/math.factorial(depleted_count))/math.factorial(guide_count-depleted_count)
                        lower_pval = round(lower_pval, 4)

                    for target_data in single_guide_dict[gene][index_key]:
                        target_name = target_data[0]
                        target_value = target_data[1]
                        if target_value <= percentile_lowerlimit:
                            outdata += "{}\t{}\t{}\t\n"\
                                .format(target_name, round(target_value, 3), lower_pval)

                        elif target_value >= percentile_upperlimit:
                            outdata += "{}\t{}\t\t{}\n"\
                                .format(target_name, round(target_value, 3), upper_pval)

                outdatafile = \
                    open("{}{}_{}_{}_Filtered_Targets.txt"
                         .format(self.args.WorkingFolder, self.args.Job_Name, sample_name, index_key), 'w')

                outdatafile.write(outdata)
                outdatafile.close()

            log_delta_control_guide_list = []
            for target_name in sample_control_dict:
                gene_name = target_name.split("_")[0]

                if guide_counts_dict2[target_name]:
                    rep_avg = gmean(guide_counts_dict2[target_name])
                    guide_counts_dict2[gene_name].append(math.log2(rep_avg/sample_control_dict[target_name]))

                if self.args.Species in target_name:
                    replicate_avg = gmean(sample_control_guide_dict[target_name])
                    delta_val = replicate_avg / sample_control_dict[target_name]
                    log_delta_control_guide_list.append(math.log2(delta_val))

            avg_delta_controls = statistics.mean(log_delta_control_guide_list)
            stdev_delta_controls = statistics.stdev(log_delta_control_guide_list)

            genedata_list = []
            upper_lower_definition = []
            gene_abundance_score = []
            file_run = datetime.datetime.today().strftime(self.date_format)
            z_outdata = "Running:\t{} Synthetic_Lethal v{}\nFile Generated:\t{}\nLibrary Control:\t{}\n" \
                        "Sample Control:\t{}\nSample:\t{}\nAvg Delta Controls:\t{}\n" \
                        "stDev Delta Controls:\t{}\nGene\tLog2\tZ pVal\tNeg Log10(Z pVal)\tKS pVal\t" \
                        "Heatmap Data\tScorable Guides\tKS Test Vals\tZ Excluded Vals\tCall\n" \
                .format(__package__, __version__, file_run, self.args.Library_Control, self.args.Control_Sample,
                        sample_name, round(avg_delta_controls, 3), round(stdev_delta_controls, 3), )

            for gene in guide_counts_dict:

                if gene == self.args.Species:
                    continue

                gene_vals = []
                excluded_guide_vals = []
                heatmap_data = 0
                depleted_pval_list = []
                enriched_pval_list = []
                up = 0
                down = 0
                vals = []
                for val in guide_counts_dict[gene]:
                    t0_pval = norm(avg_delta_controls, stdev_delta_controls).cdf(val)
                    t1_pval = 1 - t0_pval
                    depleted_pval_list.append(t0_pval)
                    enriched_pval_list.append(t1_pval)
                    vals.append(val)

                    if t0_pval <= 0.05:
                        down += 1
                    elif t1_pval <= 0.05:
                        up += 1
                    if t0_pval <= 0.05 or t1_pval <= 0.05:
                        gene_vals.append(val)
                    else:
                        excluded_guide_vals.append(round(val, 3))
                if not gene_vals:
                    gene_vals = [0]

                epval = combine_pvalues(enriched_pval_list, method='fisher', weights=None)
                dpval = combine_pvalues(depleted_pval_list, method='fisher', weights=None)
                Tool_Box.debug_messenger([sample_name, gene, down, up, dpval[1], epval[1]])
                avg_delta = statistics.mean(gene_vals)

                t0_pval = norm(avg_delta_controls, stdev_delta_controls).cdf(avg_delta)
                t1_pval = 1 - t0_pval

                try:
                    v, ks_pval = ks_2samp(log_delta_control_guide_list, guide_counts_dict2[gene])
                except RuntimeWarning:
                    ks_pval = 1

                ks_vals = []
                for v in guide_counts_dict2[gene]:
                    ks_vals.append(round(v, 3))

                choosen_z_pval = round(t0_pval, 3)
                try:
                    log10_z_pval = abs(round(math.log10(t0_pval), 3))
                except ValueError:
                    log10_z_pval = 1e-17

                if t0_pval > t1_pval:
                    choosen_z_pval = round(t1_pval, 3)
                    try:
                        log10_z_pval = abs(round(math.log10(t1_pval), 3))
                    except ValueError:
                        log10_z_pval = 1e-17

                call = False
                alpha = float(self.args.Alpha)
                if choosen_z_pval <= alpha and ks_pval <= 0.1:
                    call = True
                    heatmap_data = round(avg_delta, 3)
                guides_scored = len(gene_vals)

                z_outdata += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                    .format(gene, round(avg_delta, 3), choosen_z_pval, log10_z_pval, round(ks_pval, 3), heatmap_data,
                            guides_scored, ks_vals, excluded_guide_vals, call)

                gene_total_guides = gene_pval_dict[gene][0]
                gene_depleted_guides = gene_pval_dict[gene][1]
                gene_enriched_guides = gene_pval_dict[gene][2]
                enriched_sig_value = (self.args.LowerGuideLimit / 100) ** gene_enriched_guides
                depleted_sig_value = (self.args.LowerGuideLimit / 100) ** gene_depleted_guides
                upper_pval = "undefined"
                lower_pval = "undefined"
                gene_pval = "undefined"

                # if not gene == self.args.Species:
                # FixMe: Part of mess trying to deal with sparse data
                if gene == "Dragons":
                    upper_pval = \
                        enriched_sig_value * (
                                    (math.factorial(gene_total_guides) / math.factorial(gene_enriched_guides)) / math.factorial(
                                     gene_total_guides - gene_enriched_guides))

                    lower_pval = \
                        depleted_sig_value * (
                                    math.factorial(gene_total_guides) / math.factorial(gene_depleted_guides)) / math.factorial(
                            gene_total_guides - gene_depleted_guides)

                    if lower_pval > 1:
                        lower_pval = 1
                    if upper_pval > 1:
                        upper_pval = 1

                avg = statistics.mean(guide_counts_dict[gene])
                mdn = statistics.median(guide_counts_dict[gene])
                gene_abundance_score.append(avg)
                neg_log10 = "undefined"
                heatmap_val = "undefined"

                if avg <= 0 and lower_pval != "undefined":
                    neg_log10 = round(-1*math.log(lower_pval, 10), 4)
                    gene_pval = round(lower_pval, 4)
                    pval_used.append(lower_pval)
                    heatmap_val = -1*neg_log10
                    upper_lower_definition.append("Depleted")
                elif avg > 0 and upper_pval != "undefined":
                    neg_log10 = round(-1*math.log(upper_pval, 10), 4)
                    gene_pval = round(upper_pval, 4)
                    pval_used.append(upper_pval)
                    heatmap_val = neg_log10
                    upper_lower_definition.append("Enriched")

                if not gene == self.args.Species:
                    upper_pval = round(upper_pval, 4)
                    lower_pval = round(lower_pval, 4)

                genedata_list.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                          .format(gene, round(avg, 3), round(mdn, 3), lower_pval, upper_pval, gene_pval, neg_log10,
                                  heatmap_val))
            z_outfile = open("{}{}_{}_Z_Gene_Targets.txt"
                             .format(self.args.WorkingFolder, self.args.Job_Name, sample_name), 'w')
            z_outfile.write(z_outdata)
            z_outfile.close()

            fdr_data = stats.fdrcorrection_twostage(pval_used, alpha=float(self.args.Alpha), method="bky")

            genedata = ""

            for v1, corrected_p, null_rejection, gene_abundance, depleted_or_enriched in zip(genedata_list, fdr_data[1], fdr_data[0], gene_abundance_score, upper_lower_definition):
                if corrected_p > 1:
                    corrected_p = 1
                corrected_neg_log10 = round(-1 * math.log(corrected_p, 10), 4)

                if corrected_neg_log10 < 0:
                    corrected_neg_log10 = -1*corrected_neg_log10
                    
                corrected_heatmap_value = 0
                if corrected_p <= float(self.args.Alpha):

                    if depleted_or_enriched == "Depleted":
                        corrected_heatmap_value = -1*corrected_neg_log10
                    elif depleted_or_enriched == "Enriched":
                        corrected_heatmap_value = corrected_neg_log10

                genedata += "{}\t{}\t{}\t{}\t{}\n"\
                    .format(v1, round(corrected_p, 4), corrected_neg_log10, corrected_heatmap_value, null_rejection)

            run_stop = datetime.datetime.today().strftime(self.date_format)
            header = "Running:\t{} Synthetic_Lethal v{}\nProcess Date:\t{}\nOriginal Sig pVal (Alpha):\t{}\n" \
                     "Library Control:\t{}\nSample Control:\t{}\nSample:\t{}\n\nGene\tMean\tMedian\tDepleted pVal\t" \
                     "Enriched pVal\tpVal Used\tNeg. Log10(pVal)\t" \
                     "Heatmap Values\tCorrected pVal\tNeg. Log10(Corrected pVal)\tCorrected Heatmap Value\t" \
                     "Reject Null Hypothesis\n" \
                     .format(__package__, __version__, run_stop, self.args.Alpha, self.args.Library_Control,
                             self.args.Control_Sample, sample_name)

            genedata_outfile = open("{}{}_{}_Gene_Targets.txt"
                                    .format(self.args.WorkingFolder, self.args.Job_Name, sample_name), 'w')
            genedata_outfile.write(header+genedata)
            genedata_outfile.close()

        return

        """
            for library_index in library_index_list:
                re.sub('[\s]', "", library_index)
                self.log.debug("processing {} sample {}".format(sample_name, library_index))

                try:
                    tmp_data_file = open("{0}{1}_{2}_target_counts.txt"
                                         .format(self.args.DataFiles, self.args.Job_Name, library_index))
                except FileNotFoundError:
                    self.log.error("{0}_{1}_target_counts.txt not found".format(self.args.Job_Name, library_index))
                    raise SystemExit(1)

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

                    # Skip bad guides
                    if target_name in bad_targets_list:
                        continue

                    library_tc_norm_values[target_name]\
                        .append(library_index_target_counts[target_name]/library_index_total_count)

                    gene_name = target_name.split("_")[0]
                    if gene_name == self.args.Species:
                        key = "{}_{}".format(self.args.DataFiles, library_index)
                        sample_control_dict[key].append(library_index_target_counts[target_name] / library_index_total_count)
                library_index_target_counts.clear()

            # sample_bad_targets_dict, upper_limit, lower_limit = self.bad_correlation(library_tc_norm_values)

            # Determine the gTC_norm for each sample
            for target_name in natsort.natsorted(library_tc_norm_values):

                # Screen individual sgRNA for bad correlation between replicates
                '''
                bad_correlation = False
                for check_value in sample_bad_targets_dict[target_name]:
                    if check_value >= upper_limit or check_value <= lower_limit:
                        bad_correlation = True

                if bad_correlation and sample_name != self.args.Control_Sample:
                    tmp = []
                    for x in sample_bad_targets_dict[target_name]:
                        tmp.append("{}".format(round(x, 3)))
                    v = ";".join(tmp)
                    bad_correlation_outstring += "{}\t{}\t{}\t{}\t{}\n"\
                        .format(sample_name, target_name, v, round(upper_limit, 3), round(lower_limit, 3))
                    continue
                '''
                gene_name = target_name.split("_")[0]
                self.tc_norm_dict[sample_name][gene_name].extend(library_tc_norm_values[target_name])
                gtc_norm = statistics.mean(library_tc_norm_values[target_name])

                self.gtc_norm_dict[sample_name][target_name] = gtc_norm

                if gtc_norm == 0:
                    self.log.warning("gtc_norm=0|{} {} {}"
                                     .format(sample_name, target_name, sample_tc_data[target_name]))

            library_tc_norm_values.clear()
        '''
        # bad_correlation_outfile = open("{}{}_Bad_Correlation.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
        # bad_correlation_outfile.write(bad_correlation_outstring)
        # bad_correlation_outfile.close()
        """

    def __summary_output(self, multiprocessor_tmp_data_list):
        """
        Process temporary data file into final output.
        """

        self.log.debug("Begin Writing Data Output File")
        args = self.args

        # Data captured from the multiprocessor calls is a list of lists.  This is to remove the outer level without
        # the need to modify my existing code below.
        multiprocessor_data_list = []

        for b in multiprocessor_tmp_data_list:
            for c in b:
                multiprocessor_data_list.append(c)

        param_string = ""
        unknown_count = 0
        total_indexed_reads = 0
        file_delete_list = []
        for data in multiprocessor_data_list:

            # If we are not analyzing the unknowns there will be a NoneType entry in the list.
            try:
                index_seq = data[1]
            except TypeError:
                continue

            index_name = self.sample_data_dict[index_seq][1]
            temp_file = open("{}{}_counts.tmp".format(args.WorkingFolder, index_name), "r")
            file_delete_list.append("{}{}_counts.tmp".format(args.WorkingFolder, index_name))
            sample_indexed_reads = 0

            for line in temp_file:
                if index_name == "Unknown" or index_name == "GhostIndex":
                    unknown_count += int(line)
                else:
                    sample_indexed_reads = int(line)
                    total_indexed_reads += int(line)

            temp_file.close()

            re.sub('[\s]', "", index_name)

            # Index Name, Sample Name, Sample Replica
            param_string += "{}\t{}\t{}"\
                .format(index_name, self.sample_data_dict[index_seq][2], self.sample_data_dict[index_seq][3])

            # Index Mismatched counts
            for i in range(len(self.sample_data_dict[index_seq][0])):
                param_string += "\t{}".format(self.sample_data_dict[index_seq][0][i])

            # Targeted, Not Targeted, Fraction Targeted.
            good_reads = \
                self.sample_data_dict[index_seq][0][0]+self.sample_data_dict[index_seq][0][1]

            targeted_reads = data[0]
            if good_reads == 0:
                param_string += "\t0\t0\t0\n"
            else:
                param_string += "\t{}\t{}\t{}\n"\
                    .format(targeted_reads, sample_indexed_reads-targeted_reads, targeted_reads/sample_indexed_reads)

        # Create and populate Summary File
        index_header = "Index Name\tSample Name\tSample Replica\tIndex Mismatch 0\tIndex Mismatch 1\tFiltered Reads"

        for i in range(3):
            index_header += "\t{}_mismatches".format(i)
        index_header += "\tTargeted\tNot Targeted\tFraction Targeted\n"

        run_param_out = open("{0}{1}_summary.txt".format(args.WorkingFolder, args.Job_Name), "w")
        run_stop = datetime.datetime.today().strftime(self.date_format)
        run_param_out.write("Running:\t{} Synthetic_Lethal v{}\nStart_Time:\t{}\nStop_Time\t{}\nFASTQ_File:\t{}\n"
                            "SampleManifest:\t{}\nTarget_File:\t{}\nIndex_Mismatches\t1\nTarget_Mismatches\t{}\n"
                            "Target_Padding\t{}\nExpected_Position\t{}\nMin_Read_length\t{}\nTarget_Start\t{}\n"
                            "Target_Length\t{}\nTotal_Reads:\t{}\nIndexed_Reads:\t{}\nUnknown_Count:\t{}\n\n{}"
                            .format(__package__, __version__, self.run_start, run_stop, args.FASTQ1,
                                    args.SampleManifest, args.Target_File, args.Target_Mismatch, args.Target_Padding,
                                    args.Expected_Position, args.MinimumReadLength, args.Target_Start,
                                    args.Target_Length, self.fastq_read_counts[0], total_indexed_reads, unknown_count,
                                    index_header))

        run_param_out.write(param_string)
        run_param_out.close()

        Tool_Box.delete(file_delete_list)

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
        sample_data_dict = {}
        sample_index_mismatch = 1

        '''
        if self.args.Statistics:
            sample_index_mismatch = 1
        else:
            sample_index_mismatch = self.args.Index_Mismatch
        '''

        sample_data_dict["Unknown"] = [[0] * (sample_index_mismatch + 2), "Unknown", "Unknown", "Unknown"]

        sample_data_dict["GhostIndex"] = \
            [[0] * (sample_index_mismatch + 2), "GhostIndex", "GhostIndex", "GhostIndex"]

        with open(self.args.Master_Index_File) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                master_index_dict[l_list[0]] = "{}+{}".format(l_list[1], l_list[2])

        for sample in self.SampleManifest:
            index_name = re.sub('[\s]', '', sample[0])
            # index_name = sample[0].strip()

            if index_name in sample_data_dict:
                self.log.error("The index {0} is duplicated.  Correct the error in {1} and try again."
                               .format(sample[0], self.args.SampleManifest))
                raise SystemExit(1)

            # for each sample name append a list of all index ID's
            self.sample_mapping_dict[sample[1]].append(sample[0])
            sample_data_dict[index_name] = \
                [[0]*(sample_index_mismatch+2), sample[0], sample[1], sample[2]]

            if self.args.TargetSearch:
                fastq_file_dict[index_name] = \
                    FASTQ_Tools.Writer(self.log, "{0}{1}_{2}.fq.gz"
                                       .format(self.args.WorkingFolder, self.args.Job_Name, index_name))

                fastq_out_list.append("{0}{1}_{2}.fq.gz"
                                      .format(self.args.WorkingFolder, self.args.Job_Name, index_name))

        # This is for no index found.
        self.SampleManifest.append(("Unknown", "Unknown", "Unknown"))
        self.SampleManifest.append(("GhostIndex", "GhostIndex", "GhostIndex"))

        # If doing Statistics there is no need to run the sgRNA check.
        if self.args.Statistics:
            return sample_data_dict, fastq_file_dict, fastq_out_list, master_index_dict

        if self.args.Analyze_Unknowns:
            master_index_dict["Unknown"] = "Unknown"
            master_index_dict["GhostIndex"] = "GhostIndex"
            fastq_file_dict["Unknown"] = \
                FASTQ_Tools.Writer(self.log, "{0}{1}_Unknown.fq.gz"
                                   .format(self.args.WorkingFolder, self.args.Job_Name))

            fastq_out_list.append("{0}{1}_Unknown.fq.gz".format(self.args.WorkingFolder, self.args.Job_Name))

            fastq_file_dict["GhostIndex"] = \
                FASTQ_Tools.Writer(self.log, "{0}{1}_GhostIndex.fq.gz"
                                   .format(self.args.WorkingFolder, self.args.Job_Name))

            fastq_out_list.append("{0}{1}_GhostIndex.fq.gz".format(self.args.WorkingFolder, self.args.Job_Name))

        # Fill target list and dictionary.  Do initial quality check on target file for duplicates.
        target_list = []
        for target in Tool_Box.FileParser.indices(self.log, self.args.Target_File):
            try:
                target_seq = target[1][self.args.Target_Start:][:self.args.Target_Length]
            except ValueError:
                target_seq = target[1]

            if self.args.RevComp:
                target_seq = Sequence_Magic.rcomp(target_seq)

            target_name = target[0]
            target_list.append((target_seq, target_name))
            self.targets[len(target_seq)].append(target_seq)

            if target_seq in self.target_dict:
                self.log.error("The sgRNA target sequence in {} is duplicated.  Correct the error in {} and try "
                               "again.".format(target, self.args.Target_File))
                raise SystemExit(1)

            elif target_seq in self.target_dict[target_name]:
                self.log.error("sgRNA target name {0} is duplicated.  Correct the error in {1} and try again."
                               .format(target_name, self.args.Target_File))
                raise SystemExit(1)

            self.target_dict[target_name] = target_seq

        similarity_count = 0
        off_target_count = 0

        # Do a fine scale quality analysis of targets looking for similar sequences and off target guides.
        for target_name in self.target_dict:
            if self.args.Verbose == "INFO":
                break

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
        return sample_data_dict, fastq_file_dict, fastq_out_list, master_index_dict

    def fastq_processing(self):
        self.log.info("\033[96mDemultiplexing Input FASTQ File.\033[m")
        fastq1 = FASTQReader.Reader(self.args.FASTQ1, self.args.BatchSize)
        fastq_data = FASTQ_Tools.FastqProcessing(self.args, self.log)
        fastq_data.dataframe_build()

        # Setup Multiprocessor Workers
        p = pathos.multiprocessing.Pool(self.args.Spawn)
        worker = []
        counter = []

        for i in range(self.args.Spawn):
            worker.append(i)
            counter.append(i)

            for sample_index in self.sample_data_dict:
                # Delete tmp files if they exist
                Tool_Box.delete(["{}{}_{}.tmp".format(self.args.WorkingFolder, i, sample_index),
                                 "{}{}_GhostIndex.tmp".format(self.args.WorkingFolder, i),
                                 "{}{}_Unknown.tmp".format(self.args.WorkingFolder, i)])

        count_increment = 0
        previous_count = 0
        avg_time = []
        eof = False
        run_start = time.time()
        while not eof:
            fq1_group = []

            for i in range(self.args.Spawn):
                try:
                    read_group = next(fastq1.grouper())
                    fq1_group.append(read_group)
                    # Count total reads
                    self.fastq_read_counts[0] += len(read_group)
                except StopIteration:
                    eof = True

            count_increment += self.args.Spawn
            p.starmap(fastq_data.file_writer, zip(worker, fq1_group))

            if count_increment % 100 == 0:
                gc.collect()
                elapsed_time = int(time.time() - run_start)
                avg_time.append(elapsed_time)

                increment_completed = count_increment - previous_count
                self.log.info("{} Batches completed in {} seconds.  Avg Elapsed Time {}.  Total Reads Completed: {}".
                              format(increment_completed, elapsed_time, round(statistics.mean(avg_time), 2),
                                     count_increment*self.args.BatchSize))

                previous_count = count_increment
                run_start = time.time()

        # Delete processed FASTQ files, if they exist, so we don't add data to old ones.
        for sample_index in self.sample_data_dict:
            Tool_Box.delete(["{}{}_{}.fq.gz".format(self.args.WorkingFolder, self.args.Job_Name, sample_index),
                             "{}{}_GhostIndex.fq.gz".format(self.args.WorkingFolder, self.args.Job_Name)])

        self.log.info("FASTQ Processing Done.  Begin combining temporary files")
        file_merge_parameters = (self.args, self.log, self.sample_data_dict, worker)
        p.starmap(Tool_Box.file_merge, zip(self.sample_data_dict, itertools.repeat(file_merge_parameters)))

        for worker_id in worker:

            for sample_index in self.sample_data_dict:
                try:
                    r1_tmp_file = open("{}{}_{}.tmp".format(self.args.WorkingFolder, worker_id, sample_index), "r")
                except FileNotFoundError:
                    continue
                outstring = ""

                count_tmp_file = \
                    list(csv.reader(open("{}{}_{}_mismatch.tmp"
                                         .format(self.args.WorkingFolder, worker_id, sample_index)), delimiter='\t'))

                # Total reads for index 0 mismatch, 1 mismatch, filtered reads
                for line in count_tmp_file:
                    self.sample_data_dict[sample_index][0][0] += int(line[0])
                    self.sample_data_dict[sample_index][0][1] += int(line[1])
                    self.sample_data_dict[sample_index][0][2] += int(line[2])

                self.log.debug("Reading Temp FASTQ {}{}_{}.tmp"
                               .format(self.args.WorkingFolder, worker_id, sample_index))
                line_count = 0
                target_count = 0
                for line in r1_tmp_file:
                    outstring += line
                    line_count += 1

                if sample_index is not "Unknown":
                    self.fastq_read_counts[1] += int(line_count*0.25)
                r1_out = \
                    gzip.open("{}{}_{}.fq.gz".format(self.args.WorkingFolder, self.args.Job_Name, sample_index), "a")
                r1_out.write(outstring.encode())
                r1_out.close()

                # Close the open temp files and delete them
                r1_tmp_file.close()
                Tool_Box.delete(["{}{}_{}.tmp".format(self.args.WorkingFolder, worker_id, sample_index),
                                 "{}{}_{}_mismatch.tmp".format(self.args.WorkingFolder, worker_id, sample_index)])

    @staticmethod
    def target_search(fq_file, argvs):
        """
        I intend this to be called using multiprocessor.  Searches a demultiplexed FASTQ file for target sequences.
        :param fq_file:
        :param argvs:
        :return:
        """
        def frequency_position():
            # Total Anchors data
            freq_pos_outstring = "Position\tTotal_Anchors\tFrequency"
            freq_pos_outstring = \
                SyntheticLethal.__frequency_outstring(freq_pos_outstring, anchor_dict["total_target_pos_list"],
                                                      index_key_length, anchor_dict)
            # Total Targets Data
            freq_pos_outstring += "\nPosition\tTargets_Found\tFrequency"
            freq_pos_outstring = \
                SyntheticLethal.__frequency_outstring(freq_pos_outstring, target_found_pos_list, index_key_length,
                                                      anchor_dict)
            # No Target Data
            freq_pos_outstring += "\nPosition\tNo_Targets_Found\tFrequency"
            freq_pos_outstring = \
                SyntheticLethal.__frequency_outstring(freq_pos_outstring, no_target_pos_list, index_key_length,
                                                      anchor_dict)

            target_position_freq_outfile = open("{0}{1}_{2}_Target_Position_Freq.txt"
                                                .format(args.WorkingFolder, args.Job_Name, index_name), "w")
            target_position_freq_outfile.write(freq_pos_outstring)
            target_position_freq_outfile.close()

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
        index_name = index_dict[index_key][1]
        re.sub('[\s]', "", index_name)
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
            if len(fastq_read.seq) <= args.MinimumReadLength or fastq_read.seq.count("T") > len(fastq_read.seq)/2:
                continue

            # Find the first position of the sgRNA
            anchor_found, unknown_seq_start, anchor_dict = \
                SyntheticLethal.__anchor_search(args, fastq_read, anchor_dict)

            target_seq = False
            # Compare the sgRNA sequence to the targets.
            if anchor_found:
                target_seq, mismatch_index = \
                    SyntheticLethal.__target_match(targets_dict, fastq_read, unknown_seq_start, args)

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

        if args.Verbose == "DEBUG":
            frequency_position()

        # Format target count data for output file and write data to file.
        for line in target_file:
            target_name = line[0]

            if args.Target_Length == 'Variable':
                sgrna = line[1]
            else:
                sgrna = line[1][int(args.Target_Start):][:int(args.Target_Length)]

            target_key = sgrna
            if args.RevComp:
                target_key = Sequence_Magic.rcomp(sgrna)

            target_data_outstring += "\n{0}\t{1}".format(target_name, sgrna)
            for i in range(int(args.Target_Mismatch)+1):
                target_data_outstring += "\t{}".format(target_count_dict[target_key][i])

        target_data_file_name = "{0}{1}_{2}_target_counts.txt".format(args.WorkingFolder, args.Job_Name, index_name)
        target_data_out = open(target_data_file_name, "w")
        target_data_out.write(target_data_outstring)
        target_data_out.close()

        log.info("{} written".format(target_data_file_name))
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

        target_mismatch = args.Target_Mismatch
        targets_found_dict = collections.defaultdict(list)

        # Go through targets based on size.
        for target_length in targets_dict:

            # if args.RevComp:
            #     unknown_seq = Sequence_Magic.rcomp(fastq_read.seq[unknown_seq_start:][:target_length])
            # else:
            #     unknown_seq = fastq_read.seq[unknown_seq_start:][:target_length]

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
        # start_pos = args.AnchorStart - len(anchor_dict["index_key"])
        start_pos = args.AnchorStart

        while not anchor_found:
            unknown_seq_start = start_pos + len(args.AnchorSeq)
            mismatch_index = Sequence_Magic.match_maker(
                args.AnchorSeq, fastq_read.seq[start_pos:][:len(args.AnchorSeq)])

            # If we do not find the anchor sequence exit the loop and go to the next read.
            if start_pos > args.AnchorStop:
                unknown_seq_start = args.Expected_Position - args.Target_Padding
                anchor_dict["no_anchor_count"] += 1
                break

            elif mismatch_index <= args.AnchorMismatch:
                anchor_found = True
                anchor_dict["anchor_count"] += 1
                anchor_dict["total_target_pos_list"].append(unknown_seq_start)

            start_pos += 1
        return anchor_found, unknown_seq_start, anchor_dict

    @staticmethod
    def __frequency_outstring(freq_pos_outstring, data_list, index_key_length, anchor_dict):
        """
        This processes the data for the frequency position data file.
        :param freq_pos_outstring:
        :param data_list:
        :param index_key_length:
        :param anchor_dict:
        :return:
        """
        total_target_pos_counter = collections.Counter(data_list)

        for k in natsort.natsorted(total_target_pos_counter.items()):
            freq_pos_outstring += \
                "\n{0}\t{1}\t{2}".format(k[0]+index_key_length, k[1], round((k[1]/anchor_dict["anchor_count"]), 4))
        return freq_pos_outstring
