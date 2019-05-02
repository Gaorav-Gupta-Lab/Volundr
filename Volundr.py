#!usr/bin/env python3
"""
Volundr.py v 0.5.6
    Aug. 29, 2018
    Dennis A. Simpson
    Removed the Single Cell Processing calls.
Volundr.py v 0.5.6
    March 24, 2018
    Dennis A. Simpson
    Removed Heatmap Generator.  This is now stand alone.  Added few coding changes for clarity.
Volundr.py v 0.5.3
    July 25, 2017
    Dennis A. Simpson
    Added some code to deal with Killdevil.
Volundr.py v 0.5.0
    February 13, 2017
    Dennis A. Simpson
    Alignment launcher integrated.  Ability to do FASTQ quality and alignment using a single call.
Volundr.py v 0.3.0
    December 26, 2016
    Dennis A. Simpson
    Added call to FASTQ quality analysis.
Volundr.py v 0.2.1
    November 30, 2016
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).  Code cleaning.  Made more compact.
Volundr.py v 0.2
    August 27, 2016
    Dennis A. Simpson
    Added links to more modules and submodules.
Volundr.py v 0.15
    September 30, 2015
    Dennis A. Simpson
    Entry point for the Volundr bioinformatics package.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2016
"""

import argparse
from argparse import RawTextHelpFormatter
import time
import volundr.Synthetic_Lethal as Synthetic_Lethal
import Valkyries.Version_Dependencies as VersionDependencies
import Valkyries.Tool_Box as Tool_Box
import sys

__author__ = 'Dennis A. Simpson'
__version__ = '0.7.0'
__package__ = 'Volundr'


def main(command_line_args=None):
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

    sl = Synthetic_Lethal.SyntheticLethal(log, args)

    if eval(args.Target_Search):
        sl.fastq_analysis()
    elif eval(args.Statistics):
        sl.statistics()
    elif eval(args.Combine_Replicates):
        sl.combo_test()
    else:
        log.error('No module selected to run.')

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Volundr {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


if __name__ == '__main__':
    main()
