#!/usr/bin/env python
# Script generate MRI image ReHo map.
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

# Third-Party imports
# TODO : Resolve why importing CPAC only once does not work when importing
# it twice works!
try:
    import CPAC
except ValueError:
    print("Failed attempt at importing CPAC.")

import numpy as np
import nibabel as nb
from CPAC import version as CPAC_version
from CPAC.reho.utils import compute_reho

# Pyconnectome imports
from pyconnectome.wrapper import FSLWrapper

# Script documentation
DOC = """Run reho analysis
--------------------------

/!\ Only available for python2

# TODO : fix C-PAC installation
cd $SCRIPT_DIR/GIT_REPOS/C-PAC
python2 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_create_reho_map.py \
    -i /tmp/wrasub-03990129GPA_ses-M0_task-rest_bold_orig_ac.nii \
    -o /tmp \
    -C 27 \
    -F $HOME/fsl_init.sh \
    -V 1
or
python2 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_create_reho_map.py \
    -i $HOME/PHD/DATA/MAPT/0399/derivatives/rsfmri_preproc/sub-03990129GPA/ses-M0/wrasub-03990129GPA_ses-M0_task-rest_bold_orig_ac.nii \
    -o /tmp/test_reho \
    -F $HOME/fsl_init.sh \
    -M $HOME/PHD/DATA/MAPT/0399/derivatives/reho/sub-03990129GPA/ses-M0/wrasub-03990129GPA_ses-M0_task-rest_bold_orig_ac_mean_brain_mask.nii \
    -R $HOME/PHD/DATA/MAPT/0399/derivatives/reho/sub-03990129GPA/ses-M0/wrasub-03990129GPA_ses-M0_task-rest_bold_orig_ac_raw_reho.nii \
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
        prog="python pyphd_create_reho_map.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-im", type=is_file, required=True, metavar="<path>",
        help="Path to input image.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-C", "--cluster-size", type=int, default=27,
        help="Cluster size.")
    parser.add_argument(
        "-R", "--raw-reho-map", type=is_file, metavar="<path>",
        help="Raw reho map if it already exists.")
    parser.add_argument(
        "-M", "--mask", type=is_file, metavar="<path>",
        help="Image brain mask.")
    parser.add_argument(
        "-F", "--fsl-config",
        type=is_file, metavar="<path>",
        help="Bash script initializing FSL's environment.")
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
    "tool": "pyphd_create_reho_map.py",
    "CPAC version": CPAC_version,
    "fsl_version": FSLWrapper([], shfile=inputs["fsl_config"]).version,
    "timestamp": datetime.now().isoformat()}
outputs = {}

if verbose > 0:
    pprint("[info] Create ReHo map...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Extract brain mask if not provided
im_basename = os.path.basename(
    inputs["input_im"]).replace(".nii", "").replace(".gz", "")
if inputs["mask"] is None and inputs["raw_reho_map"] is None:
    im_brain = os.path.join(inputs["outdir"], im_basename + "_brain")
    cmd = ["bet2", inputs["input_im"], im_brain, "-m"]
    fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
    fslprocess()
    im_brain_mask = os.path.join(
        inputs["outdir"], im_basename + "_brain_mask.nii.gz")
    if not os.path.isfile(im_brain_mask):
        raise ValueError("Could not find file : {0}".format(im_brain_mask))
else:
    im_brain_mask = inputs["mask"]

# Create reho map
if inputs["raw_reho_map"] is None:
    raw_reho_map = compute_reho(
        in_file=inputs["input_im"],
        mask_file=im_brain_mask,
        cluster_size=inputs["cluster_size"])
    outputs["Raw reho map"] = raw_reho_map
else:
    raw_reho_map = inputs["raw_reho_map"]

# Create Z-score reho map
# > Compute mean image
reho_im = nb.load(raw_reho_map)
mask_im = nb.load(im_brain_mask)
mean_reho = np.mean(reho_im.get_data()[mask_im.get_data() == 1])

# > Compute std image
std_reho = np.std(reho_im.get_data()[mask_im.get_data() == 1])

# > Substract mean and divide by std
z_scored_reho_map = os.path.join(
    inputs["outdir"], inputs["input_im"].replace(".nii", "") + "_zscored_reho")
cmd = ["fslmaths", raw_reho_map, "-sub", mean_reho, "-div", std_reho,
       z_scored_reho_map]
fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
fslprocess()

if verbose > 0:
    print("[Outputs]:")
    print(outputs)
