#!/usr/bin/env python
# Script used to extract Conn 1rst level analysis connectivity scores.
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

# System import
import os
import argparse
import textwrap
import json
from datetime import datetime
from argparse import RawTextHelpFormatter
from datetime import datetime
from pprint import pprint

# Third party imports
import pandas as pd
from numpy import tanh
from scipy.io import loadmat

# Script documentation
DOC = """
Extract Conn 1rst level analysis subjects connectivity scores.
--------------------------------------------------------------

1) Extract Z-score matrix
2) Compute r-scores from Z-score matrix

Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_get_conn_first_level_scores.py \
    -m resultsROI_Subject001_Condition001.mat \
    -o /tmp/test_conn \
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
        prog="python pyphd_get_conn_first_level_scores.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-m", "--conn-mat-file", type=is_file, required=True,
        help="Conn mat files from 1rst level analysis results. "
             "Must be subject-condition specific file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
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
    "timestamp": datetime.now().isoformat(),
    "tool": "pyphd_get_conn_first_level_scores"
}
outputs = {}
if verbose > 0:
    pprint("[info] Get connectivity scores...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Load data
conn_data = loadmat(inputs["conn_mat_file"])

# Create z-score dataframe
# > Get row and column names
source_names = [x[0] for x in conn_data["names"][0]]
target_names = [x[0] for x in conn_data["names2"][0]]
zscores_df = pd.DataFrame(
    conn_data["Z"], index=source_names, columns=target_names)

# Compute r-scores
rscores_df = zscores_df.apply(tanh)

# Save outputs
output_basename = os.path.basename(inputs["conn_mat_file"]).replace(".mat", "")
out_zscores = os.path.join(
    inputs["outdir"], output_basename + "_zscores" + ".txt")
zscores_df.to_csv(out_zscores, sep=",")
out_rscores = os.path.join(
    inputs["outdir"], output_basename + "_rscores" + ".txt")
rscores_df.to_csv(out_rscores, sep=",")
outputs["Zscore file"] = out_zscores
outputs["Rscore file"] = out_rscores


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(
        logdir,
        "pyphd_get_conn_first_level_scores_{0}_{1}.json".format(
            output_basename, name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 1:
    print("[final]")
    pprint(outputs)
