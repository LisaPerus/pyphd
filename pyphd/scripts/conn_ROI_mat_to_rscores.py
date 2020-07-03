#!/usr/bin/env python
# Script used to convert ROI-to-ROI Conn 1rst level analysis connectivity
# zscores to rscores
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
Script used to convert ROI-to-ROI Conn 1rst level analysis connectivity zscores
to rscores
-------------------------------------------------------------------------------

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/conn_ROI_mat_to_rscores.py \
    -i resultsROI_Subject001_Condition001.mat \
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
        prog="python conn_ROI_mat_to_rscores.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-mat-file", type=is_file, required=True,
        help="Input conn ROI mat file.")
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
    "tool": "conn_ROI_mat_to_rscores.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Convert conn ROI mat zscores to rscores...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Load data
input_roi_data = loadmat(inputs["input_mat_file"])
input_roi_source_names = [x[0] for x in input_roi_data["names"][0]]
input_roi_target_names = [x[0] for x in input_roi_data["names2"][0]]
input_roi_basename = os.path.basename(
    inputs["input_mat_file"]).replace(".mat", "")

# Get z-scores
input_roi_zscores_df = pd.DataFrame(
    input_roi_data["Z"], index=input_roi_source_names,
    columns=input_roi_target_names)

# Z-score the difference of R-scores
# > Get r-scores
input_roi_rscores_df = input_roi_zscores_df.apply(tanh)

# Save outputs
roi_mat_data = input_roi_data.copy()
roi_mat_data["Z"] = input_roi_rscores_df.as_matrix()
out_mat = os.path.join(
    inputs["outdir"], "rscores_{0}.mat".format(input_roi_basename))
savemat(out_mat, roi_mat_data)
outputs["Rscores mat file"] = out_mat


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
        "conn_ROI_mat_to_rscores_{0}_{1}.json".format(
            input_roi_basename, name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
