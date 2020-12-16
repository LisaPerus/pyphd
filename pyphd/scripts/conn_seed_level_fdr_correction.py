#!/usr/bin/env python
# Script used to run Conn seed level correction
# Copyright (C) 2020  Lisa Perus
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# long with this program.  If not, see <https://www.gnu.org/licenses/>.

# System import
import os
import argparse
import textwrap
import json
import glob
from datetime import datetime
from argparse import RawTextHelpFormatter
from pprint import pprint

# Pyphd imports
from pyphd.fc_analyses.conn import (
    parse_conn_roi_to_roi_output_textfile, conn_seed_level_fdr_correction)


# Script documentation
DOC = """

Script used to run Conn seed level correction.
----------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/conn_seed_level_fdr_correction.py \
    -i /tmp/conn_results.txt \
    -s F \
    -o /tmp \
    -E 35

"""


def is_file(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError("File does not exist: %s" % filepath)
    return filepath


def is_directory(dirarg):
    """ Type for argparse - checks that directory exists.
    """
    if not os.path.isdir(dirarg):
        raise argparse.ArgumentError(
            "The directory '{0}' does not exist!".format(dirarg))
    return dirarg


# Parse input arguments
def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python conn_seed_level_fdr_correction.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-conn-all-connections", type=is_file, required=True,
        help="File exported from Conn ROI-to-ROI analysis gui containing all "
             "connections test pval.")
    required.add_argument(
        "-s", "--test-statistic", type=str, required=True,
        help="Test statistic used (T or F).")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-C", "--conn-version", type=str, default="Conn19c",
        help="Conn version used to get input file.")
    parser.add_argument(
        "-E", "--expected-nb-conns-per-seed", type=int,
        help="Expected number of connections found for each seed.")
    parser.add_argument(
        "-V", "--verbose",
        type=int, choices=[0, 1], default=1,
        help="Increase the verbosity level: 0 silent, 1 verbose.")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")
    return kwargs, verbose


"""
Parse the command line.
"""
inputs, verbose = get_cmd_line_args()
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "conn_seed_level_fdr_correction.py"
}
outputs = {}


"""
Reorganize input data with all connections
"""
reorganized_input = os.path.join(
    inputs["outdir"], "all_connections_reorganized.txt")
parse_conn_roi_to_roi_output_textfile(
    conn_textfile=inputs["input_conn_all_connections"],
    method="0_ALL_CONNECTIONS",
    statistic=inputs["test_statistic"],
    conn_version=inputs["conn_version"],
    outfile=reorganized_input)
outputs["Reorganized all connections"] = reorganized_input


"""
Run seed level FDR correction
"""
output_file = os.path.join(
    inputs["outdir"], "all_conns_seed_fdr_corrected.txt")
output_file_sig_only = os.path.join(
    inputs["outdir"], "significant_only_conns_seed_fdr_corrected.txt")
seed_corr_conn_file, sig_seed_corr_conn_file = conn_seed_level_fdr_correction(
    conn_parsed_textfile=reorganized_input,
    outfile=output_file,
    outfile_threshold_alpha=output_file_sig_only,
    alpha=0.05,
    expected_nb_conns=inputs["expected_nb_conns_per_seed"])
outputs["All FDR seed corrected conn file"] = seed_corr_conn_file
outputs["Significant FDR seed corrected conn file"] = sig_seed_corr_conn_file


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)

output_basename = "conn_seed_level_fdr_correction"
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
