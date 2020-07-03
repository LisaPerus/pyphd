#!/usr/bin/env python
# Script used to create substract one image to another using FSL
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
from datetime import datetime
from pprint import pprint
from collections import OrderedDict
from pyconnectome.wrapper import FSLWrapper

# Script documentation
DOC = """

Script used to substract one image to another image
---------------------------------------------------

This script is used to substract one image to another image. It
is used for instance to have a delta scalar on TBSS skeleton for two
timepoints.

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/create_diff_image.py \
    -i ${MEDIA_SCRIPT}/MAPT_DTI/MAPT_RERUN/tbss/M0/FA_individ/01310009GEU/stats/01310009GEU_masked_FAskel.nii.gz \
    -s ${MEDIA_SCRIPT}/MAPT_DTI/MAPT_RERUN/tbss/M36/FA_individ/01310009GEU/stats/01310009GEU_masked_FAskel.nii.gz \
    -f ~/fsl_init.sh \
    -o /tmp \
    -M ${MEDIA_SCRIPT}/MAPT_DTI/MAPT_RERUN/tbss/common_mask_M0_M36_M60_mean_FA_skel.nii.gz \
    -N 01310009GEU_masked_FAskel.nii.gz \
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
        prog="python create_diff_image.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--first-im", type=is_file, required=True,
        help="First image. This image will be substracted to the second one.")
    required.add_argument(
        "-s", "--second-im", type=is_file, required=True,
        help="Second image. The first image will be substracted to this one.")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-M", "--mask", type=is_file,
        help="File to mask final image with. Used if images at different "
             "timepoints do not have the same shape. E.g : tbss skeleton "
             "at different timepoints need to be masked by common mask.")
    parser.add_argument(
        "-N", "--outname", type=str, default="out.nii",
        help="File outname.")
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
    "tool": "create_diff_image.py"
}
outputs = {}

# Substract first image to the second
out_nii = os.path.join(inputs["outdir"], inputs["outname"])
cmd = ["fslmaths", inputs["second_im"], "-sub", inputs["first_im"], out_nii]
fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
fslprocess()

# Mask image
if inputs["mask"]:
    cmd = ["fslmaths", out_nii, "-mas", inputs["mask"], out_nii]
    fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
    fslprocess()
outputs["Diff image"] = out_nii

"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)

output_basename = "create_diff_image"
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
