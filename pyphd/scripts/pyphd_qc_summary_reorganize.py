#!/usr/bin/env python
# Script to reorganize QC Summary data
# Copyright (C) 2019  Lisa Perus
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

# System imports
import os
import re
from pprint import pprint
import argparse
from argparse import RawTextHelpFormatter
import textwrap
from datetime import datetime

# Pyphd imports
from pyphd.utils.excel import excel_to_csv
from pyphd.utils.excel import compute_excel_index

# Script documentation
DOC = """Convert Cati xls QC Summary file to csv
------------------------------------------------

Example on MEMENTO QC Summary file:

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_qc_summary_reorganize.py \
    -q QC_summary_MEMENTO_M024_ongoing.xls \
    -c A AU DS\
    -o $HOME/PHD/NOTES/INFOS_MAPT_CATI/memento \
    -e Sid T1_note Func_note \
    -S 9 \
    -V 2
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
        prog="python pyphd_qc_summary_reorganize.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-q", "--qc-summary-file",
        type=is_file, required=True, metavar="<path>",
        help="Path to the QC summary file.")
    required.add_argument(
        "-c", "--columns", type=str, nargs='+', required=True,
        help="Excel indexes of the columns to keep in QC summary file "
             "e.g : AU")
    parser.add_argument(
        "-e", "--output-header", type=str, nargs='+', required=True,
        help="Output header elements.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-H", "--intermediate-header", type=str, nargs='+', default=[0],
        help="Original header in xls file.")
    parser.add_argument(
        "-S", "--skip-row", type=int, default=0,
        help="Skip rows under this value.")
    parser.add_argument(
        "-V", "--verbose",
        type=int, choices=[0, 1, 2], default=2,
        help="Increase the verbosity level: 0 silent, [1, 2] verbose.")

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
    "tool": "pyphd_qc_summary_reorganize.py",
    "timestamp": datetime.now().isoformat()}
outputs = {}

if verbose > 0:
    pprint("[info] Starting reorganizing QC Summary data...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Step 1 : Convert to csv files
"""
qc_summary_dir = os.path.join(inputs["outdir"], "QC_summary")
if not os.path.isdir(qc_summary_dir):
    os.mkdir(qc_summary_dir)
csv_files = excel_to_csv(
    xls_file=inputs["qc_summary_file"], outdir=qc_summary_dir,
    skip_row=inputs["skip_row"], header=inputs["intermediate_header"])


"""
Step 2 : Concatenate csv files
"""
cols_indexes = []
for col in inputs["columns"]:
    idx, excel_columns = compute_excel_index(col)
    cols_indexes.append(idx)
max_idx = 0
for elt in cols_indexes:
    if elt > max_idx:
        max_idx = elt

data = []
for idx, fid in enumerate(csv_files):
    print("Reading file : '{0}'...".format(fid))
    with open(fid, "rt") as open_file:
        lines = open_file.readlines()
    for line in lines[inputs["skip_row"]:]:
        line = line.split(",")
        if len(line) <= max_idx:
            print(
                "Line {0} ( {1} ) did not split correctly, add it "
                "manually.".format(idx, line))
            continue
        data.append([line[idx_col] for idx_col in cols_indexes])

"""
Step 3 : Write output
"""
output_basename = os.path.basename(inputs["qc_summary_file"])
output_basename = re.sub("\\..*$", "", output_basename)
output_csv = os.path.join(inputs["outdir"], output_basename + ".csv")
with open(output_csv, "wt") as open_file:
    open_file.write(",".join(inputs["output_header"]))
    open_file.write("\n")
    for line in data:
        open_file.write(",".join(line))
        open_file.write("\n")
outputs["QC_summary_csv"] = output_csv

if verbose > 0:
    pprint("Outputs")
    pprint(outputs)
