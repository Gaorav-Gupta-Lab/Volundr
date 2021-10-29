#!usr/bin/env python3

"""
Volundr.py v 3.0.0
    Entry point for the Volundr bioinformatics package.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2021
"""
import ntpath
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import time
import pathlib
import magic
from distutils.util import strtobool
import volundr.Synthetic_Lethal as Synthetic_Lethal
import Valkyries.Version_Dependencies as VersionDependencies
import Valkyries.Tool_Box as Tool_Box

from scipy import stats
import math
from numpy import log as ln

__author__ = 'Dennis A. Simpson'
__version__ = '3.0.0'
__package__ = 'Völundr'


def main():
    """

    """
    VersionDependencies.python_check()

    parser = argparse.ArgumentParser(description="A package to process Synthetic Lethal Data.\n {0} v{1}"
                                     .format(__package__, __version__), formatter_class=RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    # Convert strings to int, float, boolean, check file names and paths for errors
    args, log = error_checking(parser)
    start_time = time.time()

    # Initialize program
    synthetic_lethal = Synthetic_Lethal.SyntheticLethal(log, args)

    if args.TargetSearch:
        module_name = "Target Search"
        log.info("{} v{}; Module: {} v{} Beginning"
                 .format(__package__, __version__, module_name, Synthetic_Lethal.__version__))

        synthetic_lethal.fastq_analysis()

    elif args.Statistics:
        module_name = "Statistical Analysis"
        log.info("{} v{}; Module: {} v{} Beginning"
                 .format(__package__, __version__, module_name, Synthetic_Lethal.__version__))
        synthetic_lethal.statistics()

    else:
        module_name = "No module selected"
        log.error('No module selected to run.')

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Völundr {0} complete ({1} seconds, {2} Mb peak memory).****\n{3}"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


def error_checking(parser):
    """
    Check parameter file for errors and return parser object.
    :param parser:
    :return:
    """
    def string_conversions(parser):
        """
        Convert True/False statements in parameter file to boolean
        :param parser:
        :return:
        """
        options_parser = Tool_Box.options_file(parser)
        initial_args = options_parser.parse_args()

        options_parser.set_defaults(TargetSearch=bool(strtobool(initial_args.TargetSearch)))
        options_parser.set_defaults(Statistics=bool(strtobool(initial_args.Statistics)))

        options_parser.set_defaults(Verbose=initial_args.Verbose.upper())

        if initial_args.Statistics == "False":
            options_parser.set_defaults(AnchorSeq=initial_args.AnchorSeq.upper())
            options_parser.set_defaults(Analyze_Unknowns=bool(strtobool(initial_args.Analyze_Unknowns)))
            options_parser.set_defaults(Delete_Demultiplexed_FASTQ=bool(strtobool(initial_args.Delete_Demultiplexed_FASTQ)))
            options_parser.set_defaults(RevComp=bool(strtobool(initial_args.RevComp)))
            options_parser.set_defaults(BatchSize=int(initial_args.BatchSize))
            options_parser.set_defaults(Target_Mismatch=int(initial_args.Target_Mismatch))
            options_parser.set_defaults(MinimumReadLength=int(initial_args.MinimumReadLength))
            options_parser.set_defaults(N_Limit=10)
            options_parser.set_defaults(Target_Length=int(initial_args.Target_Length))
            options_parser.set_defaults(Target_Start=int(initial_args.Target_Start))
            # options_parser.set_defaults(Index_Mismatch=int(initial_args.Index_Mismatch))
            options_parser.set_defaults(Spawn=int(initial_args.Spawn))
            options_parser.set_defaults(Target_Padding=int(initial_args.Target_Padding))
            options_parser.set_defaults(Expected_Position=int(initial_args.Expected_Position))
            options_parser.set_defaults(AnchorMismatch=int(initial_args.AnchorMismatch))
            options_parser.set_defaults(AnchorStart=int(initial_args.AnchorStart))
            options_parser.set_defaults(AnchorStop=int(initial_args.AnchorStop))
        else:
            options_parser.set_defaults(Write_TDnorm_Log2_sgRNA_Control_File=
                                        bool(strtobool(initial_args.Write_TDnorm_Log2_sgRNA_Control_File)))
            options_parser.set_defaults(Write_TDnorm_Log2_sgRNA_Sample_File=
                                        bool(strtobool(initial_args.Write_TDnorm_Log2_sgRNA_Sample_File)))
            options_parser.set_defaults(Write_Log2_sgRNA_File=
                                        bool(strtobool(initial_args.Write_Log2_sgRNA_File)))
            options_parser.set_defaults(Write_Permuted_Log2_Data_File=
                                        bool(strtobool(initial_args.Write_Permuted_Log2_Data_File)))
            options_parser.set_defaults(Bad_sgRNA_Lower_Percentile=float(initial_args.Bad_sgRNA_Lower_Percentile))
            options_parser.set_defaults(Bad_sgRNA_Upper_Percentile=float(initial_args.Bad_sgRNA_Upper_Percentile))
            options_parser.set_defaults(UpperPercentile=float(initial_args.UpperPercentile))
            options_parser.set_defaults(LowerPercentile=float(initial_args.LowerPercentile))
            options_parser.set_defaults(PermutationCount=int(initial_args.PermutationCount))
            options_parser.set_defaults(Alpha=float(initial_args.Alpha))
            options_parser.set_defaults(Target_Mismatch=float(initial_args.Target_Mismatch))
            options_parser.set_defaults(UpperGuideLimit=float(initial_args.UpperGuideLimit))
            options_parser.set_defaults(LowerGuideLimit=float(initial_args.LowerGuideLimit))

        initial_args = options_parser.parse_args()

        return initial_args

    args = string_conversions(parser)
    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, sys.argv)

    if not pathlib.Path(args.WorkingFolder).exists():
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Parameter File."
              .format(args.WorkingFolder))
        raise SystemExit(1)

    if args.Statistics:
        if not pathlib.Path(args.DataFiles).exists():
            print("\033[1;31mERROR:\n\t--DataFiles Folder Path: {} Not Found.  Check Parameter File."
                  .format(args.DataFiles))
            raise SystemExit(1)

    if not pathlib.Path(args.SampleManifest).exists():
        print("\033[1;31mERROR:\n\t--SampleManifest: {} Not Found.  Check Parameter File."
              .format(args.SampleManifest))
        raise SystemExit(1)

    if not pathlib.Path(args.Master_Index_File).exists():
        print("\033[1;31mERROR:\n\t--Master_Index_File: {} Not Found.  Check Parameter File."
              .format(args.Master_Index_File))
        raise SystemExit(1)

    if not pathlib.Path(args.Target_File).exists():
        print("\033[1;31mERROR:\n\t--Target_File: {} Not Found.  Check Parameter File."
              .format(args.Target_File))
        raise SystemExit(1)

    if args.TargetSearch:
        if getattr(args, "FASTQ1", False) and not pathlib.Path(args.FASTQ1).exists():
            print("\033[1;31mERROR:\n\t--FASTQ1: {} Not Found.  Check Parameter File."
                  .format(args.FASTQ1))
            raise SystemExit(1)

        try:
            mime_type1 = magic.from_file(args.FASTQ1, mime=True).decode()

        except AttributeError:
            mime_type1 = magic.from_file(args.FASTQ1, mime=True)

        if "text" in mime_type1 or "gzip" in mime_type1:
            pass
        else:
            log.error("Unsupported FASTQ file-type.  Only TEXT or GZIP Allowed.")
            raise SystemExit(1)

    return args, log


if __name__ == '__main__':
    '''
    # This is some code we have been using to test different stats approaches for the analysis when the data is sparse.
    # chi_sq = np.sum(np.divide(np.square(observation - expectation), expectation))
    # stats.chisquare(f_obs=observation, f_exp=expectation, ddof=0)
    v = [0.031576546525932093, 0.02067054709516096, 0.9990222022222689, 0.10449493367883328, 0.019424285019606614, 0.9990222022222689, 0.34240059947767854, 0.10172969102060508, 0.2673538627084865, 0.03968862219705106]
    v0 = [0.968423453, 0.979329453, 0.000977798, 0.895505066, 0.980575715, 0.000977798, 0.657599401, 0.898270309, 0.732646137, 0.960311378]
    v1 = -2*ln(0.05)
    v2 = -2*ln(0.01)
    v3 = -2*math.log(0.05, math.e)
    v4 = -2*math.log(0.01, math.e)
    print(v1+v2, stats.chi2.cdf(v1+v2, 4), v3+v4, stats.chi2.cdf(v3+v4, 3))

    fstat = stats.combine_pvalues([v1, v2], method='fisher')
    print(fstat, stats.combine_pvalues(v0, method='fisher', weights=None))
    # [1 - stats.chi2.cdf(chi_sq, len(observation) - 1),
    # stats.chisquare(f_obs=observation, ddof=0)[1]]
    '''
    main()
