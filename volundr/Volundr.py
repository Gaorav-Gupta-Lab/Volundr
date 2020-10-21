#!usr/bin/env python3
"""
Volundr.py v 1.0.0
    Entry point for the Volundr bioinformatics package.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""

import sys
import argparse
from argparse import RawTextHelpFormatter
import time
import volundr.Synthetic_Lethal as Synthetic_Lethal
import Valkyries.Version_Dependencies as VersionDependencies
import Valkyries.Tool_Box as Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = '1.0.1'
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

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()

    # If we are doing statistical analysis the user will not input an Index_Mismatch value
    if not getattr(args, "Index_Mismatch", False):
        options_parser.add_argument("--Index_Mismatch", dest="Index_Mismatch", default=0)
        options_parser.add_argument("--Analyze_Unknowns", dest="Analyze_Unknowns", default="False")
        args = options_parser.parse_args()

    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)
    start_time = time.time()
    module_name = "Synthetic_Lethal"

    log.info("{0} v{1}; Module: Synthetic Lethal Analysis v{2} Beginning"
             .format(__package__, __version__, Synthetic_Lethal.__version__))

    # Convert universal variables intended as boolean from string to boolean.
    # ToDo: Should be a cleaner method to do this.
    if args.Target_Search == "True":
        options_parser.set_defaults(Target_Search=True)
        if args.RevComp == "True":
            options_parser.set_defaults(RevComp=True)
        else:
            options_parser.set_defaults(RevComp=False)
        if args.Delete_Demultiplexed_FASTQ == "True":
            options_parser.set_defaults(Delete_Demultiplexed_FASTQ=True)
        else:
            options_parser.set_defaults(Delete_Demultiplexed_FASTQ=False)
        if args.compress == "True":
            options_parser.set_defaults(compress=True)
        else:
            options_parser.set_defaults(compress=False)
    else:
        options_parser.set_defaults(Target_Search=False)

    if args.Statistics == "True":
        options_parser.set_defaults(Statistics=True)
    else:
        options_parser.set_defaults(Statistics=False)

    if args.Analyze_Unknowns == "True":
        options_parser.set_defaults(Analyze_Unknowns=True)
    else:
        options_parser.set_defaults(Analyze_Unknowns=False)

    args = options_parser.parse_args()

    synthetic_lethal = Synthetic_Lethal.SyntheticLethal(log, args)

    # Add some parameters to our options parser object.
    args = options_parser.parse_args()

    if args.Target_Search:
        synthetic_lethal.fastq_analysis()
    elif args.Statistics:
        synthetic_lethal.statistics()
    else:
        log.error('No module selected to run.')

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Volundr {0} complete ({1} seconds, {2} Mb peak memory).****\n{3}"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


if __name__ == '__main__':
    main()
