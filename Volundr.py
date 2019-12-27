#!usr/bin/env python3
# coding=utf-8
"""
Volundr.py v 1.0.0
    Entry point for the Volundr bioinformatics package.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import time
from distutils.util import strtobool
import volundr.Synthetic_Lethal as Synthetic_Lethal
import Valkyries.Version_Dependencies as VersionDependencies
import Valkyries.Tool_Box as Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = '2.x.x'
__package__ = 'Völundr'


def main(command_line_args=None):
    """

    :param command_line_args:
    """
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv

    parser = argparse.ArgumentParser(description="A package to process Synthetic Lethal Data.\n {0} v{1}"
                                     .format(__package__, __version__), formatter_class=RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    # Convert universal variables intended as boolean from string to boolean.
    args, options_parser = string_to_boolean(Tool_Box.options_file(parser))

    # Check file names and paths for errors
    error_checking(args)

    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)

    start_time = time.time()
    module_name = "Synthetic_Lethal"

    log.info("{0} v{1}; Module: Synthetic Lethal Analysis v{2} Beginning"
             .format(__package__, __version__, Synthetic_Lethal.__version__))

    synthetic_lethal = Synthetic_Lethal.SyntheticLethal(log, args)

    if args.TargetSearch:
        synthetic_lethal.fastq_analysis()
    elif args.Statistics:
        synthetic_lethal.statistics()
    else:
        log.error('No module selected to run.')

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Völundr {0} complete ({1} seconds, {2} Mb peak memory).****\n{3}"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


def error_checking(args):
    """
    Make sure all paths and files exist.
    :param args:
    """
    if not os.path.exists(args.Working_Folder):
        print("\033[1;31mERROR:\n\t--Working_Folder: {} Not Found.  Check Options File."
              .format(args.Working_Folder))
        raise SystemExit(1)

    if not os.path.isfile(args.Target_File):
        print("\033[1;31mERROR:\n\t--Target_File: {} Not Found.  Check Options File."
              .format(args.Target_File))
        raise SystemExit(1)

    if not os.path.isfile(args.Master_Index_File):
        print("\033[1;31mERROR:\n\t--Master_Index_File: {} Not Found.  Check Options File."
              .format(args.Master_Index_File))
        raise SystemExit(1)

    if not os.path.isfile(args.SampleManifest):
        print("\033[1;31mERROR:\n\t--Index_File: {} Not Found.  Check Options File."
              .format(args.Index_File))
        raise SystemExit(1)


def string_to_boolean(options_parser):
    """
    Converts strings to boolean.  Done to keep the eval() function out of the code.
    :param options_parser:
    :return:
    """
    args = options_parser.parse_args()

    options_parser.set_defaults(TargetSearch=bool(strtobool(args.TargetSearch)))
    options_parser.set_defaults(Statistics=bool(strtobool(args.Statistics)))

    if not getattr(args, "Index_Mismatch", False):
        options_parser.add_argument("--Index_Mismatch", dest="Index_Mismatch", default=0)
        options_parser.add_argument("--Analyze_Unknowns", dest="Analyze_Unknowns", default=False)
        options_parser.set_defaults(Write_TDnorm_Log2_sgRNA_Control_File=
                                    bool(strtobool(args.Write_TDnorm_Log2_sgRNA_Control_File)))
        options_parser.set_defaults(Write_TDnorm_Log2_sgRNA_Sample_File=
                                    bool(strtobool(args.Write_TDnorm_Log2_sgRNA_Sample_File)))
        options_parser.set_defaults(Write_Log2_sgRNA_File=
                                    bool(strtobool(args.Write_Log2_sgRNA_File)))
        options_parser.set_defaults(Write_Permuted_Log2_Data_File=
                                    bool(strtobool(args.Write_Permuted_Log2_Data_File)))

    if args.TargetSearch == "True":
        options_parser.set_defaults(Analyze_Unknowns=bool(strtobool(args.Analyze_Unknowns)))
        options_parser.set_defaults(RevComp=bool(strtobool(args.RevComp)))
        options_parser.set_defaults(Delete_Demultiplexed_FASTQ=bool(strtobool(args.Delete_Demultiplexed_FASTQ)))
        options_parser.set_defaults(Compress=bool(strtobool(args.Compress)))

    args = options_parser.parse_args()

    return args, options_parser


if __name__ == '__main__':
    main()
