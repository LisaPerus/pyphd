#!/usr/bin/env python
# Script used to substract two ROI-to-ROI Conn 1rst level analysis connectivity
# scores mat files.
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
import progressbar

# Third party imports
import pandas as pd
import numpy as np
from numpy import tanh, arctanh
from scipy.io import loadmat, savemat

# Script documentation
DOC = """
Script used to substract two ROI-to-ROI Conn 1rst level analysis connectivity
scores.
-----------------------------------------------------------------------------

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/substract_conn_ROI_mat.py \
    -f resultsROI_Subject044_Condition001.mat \
    -s resultsROI_Subject044_Condition002.mat \
    -o /tmp \
    -V 1
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
        prog="python pyphd_conn_rtr_to_palm.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-f", "--first-mat-file", type=is_file, required=True,
        help="First conn ROI mat file. Substraction is Second - First.")
    required.add_argument(
        "-s", "--second-mat-file", type=is_file, required=True,
        help="Conn ROI mat file. Substraction is Second - First.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
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
    "tool": "substract_conn_ROI_mat.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Substract conn ROI mat files...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Load data
first_roi_data = loadmat(inputs["first_mat_file"])
second_roi_data = loadmat(inputs["second_mat_file"])
first_roi_source_names = [x[0] for x in first_roi_data["names"][0]]
first_roi_target_names = [x[0] for x in first_roi_data["names2"][0]]
second_roi_source_names = [x[0] for x in second_roi_data["names"][0]]
second_roi_target_names = [x[0] for x in second_roi_data["names2"][0]]
first_roi_basename = os.path.basename(
    inputs["first_mat_file"]).replace(".mat", "")
second_roi_basename = os.path.basename(
    inputs["second_mat_file"]).replace(".mat", "")

if first_roi_source_names != second_roi_source_names:
    raise ValueError(
        "Different source names for {0} and {1}!".format(
            inputs["first_mat_file"], inputs["second_mat_file"]))
if first_roi_target_names != second_roi_target_names:
    raise ValueError(
        "Different target names for {0} and {1}!".format(
            inputs["first_mat_file"], inputs["second_mat_file"]))

# Get z-scores
first_roi_zscores_df = pd.DataFrame(
    first_roi_data["Z"], index=first_roi_source_names,
    columns=first_roi_target_names)
second_roi_zscores_df = pd.DataFrame(
    second_roi_data["Z"], index=second_roi_source_names,
    columns=second_roi_target_names)

# Get r-scores
first_roi_rscores_df = first_roi_zscores_df.apply(tanh)
second_roi_rscores_df = second_roi_zscores_df.apply(tanh)

# Get difference
diff_rscores_df = second_roi_rscores_df - first_roi_rscores_df

# > Save diff r-scores
out_csv = os.path.join(
    inputs["outdir"], "rscores_{0}-{1}.csv".format(
        first_roi_basename, second_roi_basename))
diff_rscores_df.to_csv(out_csv)
outputs["Rscores diff csv file"] = out_csv

# Z-score the difference
diff_zscores_df = diff_rscores_df.apply(arctanh)

# Save data
diff_roi_mat_data = first_roi_data.copy()
diff_roi_mat_data["Z"] = diff_zscores_df.as_matrix()

out_mat = os.path.join(
    inputs["outdir"], "zscores_{0}-{1}.mat".format(
        first_roi_basename, second_roi_basename))
savemat(out_mat, diff_roi_mat_data)
outputs["Diff mat file"] = out_mat

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
        "substract_conn_ROI_mat_{0}_{1}_{2}.json".format(
            first_roi_basename, second_roi_basename, name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 1:
    print("[final]")
    pprint(outputs)
