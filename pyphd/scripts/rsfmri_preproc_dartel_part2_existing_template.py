# Script reproducing Clara Manesco matlab script for rsfMRI data preprocessing
# on the MAPT cohort with Dartel registration and an existing Dartel template.
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
import glob
import re
import json
import shutil
import subprocess
import time
from datetime import datetime
from argparse import RawTextHelpFormatter
from datetime import datetime
from pprint import pprint

# Pyphd imports
from pyphd.spm.utils import run_dartel_from_existing_template

# Script documentation
DOC = """
Execute Clara Manesco fMRI preprocessing pipeline initially designed to be run
on MAPT data with DARTEL registration and an existing Dartel template
(from Template 1 to 6).
------------------------------------------------------------------------------

Example:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/rsfmri_preproc_dartel_part2_existing_template.py \
    -g /tmp/sub-03990228BJA/rc1sub-03990228BJA_ses-M36_T1w_orig_ac.nii \
    -w /tmp/sub-03990228BJA/ses-M36/rc2sub-03990228BJA_ses-M36_T1w_orig_ac.nii \
    -t /tmp/sub-03990228BJA/Template_1.nii \
       /tmp/sub-03990228BJA/Template_2.nii \
       /tmp/sub-03990228BJA/Template_3.nii \
       /tmp/sub-03990228BJA/Template_4.nii \
       /tmp/sub-03990228BJA/Template_5.nii \
       /tmp/sub-03990228BJA/Template_6.nii \
    -o /tmp/sub-03990228BJA/ses-M36/ \
    -S $HOME/opt/spm12/spm12/run_spm12.sh \
    -M $HOME/opt/spm12/mcr/v713 \
    -V 2
"""


def spm_get_version(spm_sh, spm_mcr):
    """Returns SPM version
    """
    cmd = [spm_sh, spm_mcr, "--version"]
    version = subprocess.check_output(cmd).decode("utf-8")
    version = version.split("\n")[0:2]
    version = " : ".join(version)
    return version


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
        prog="python rsfmri_preproc_dartel_existing_template.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-g", "--gm-seg", type=is_file, required=True,
        help="Path to subject GM segmentation FOR DARTEL (rc1 file).")
    required.add_argument(
        "-w", "--wm-seg", type=is_file, required=True,
        help="Path to subject WM segmentation FOR DARTEL (rc2 file).")
    required.add_argument(
        "-t", "--templates", type=is_file, required=True, nargs="+",
        help="Dartel template files (Template_1..6.nii).")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-S", "--spm-sh", type=is_file, help="Path to SPM sh file.")
    parser.add_argument(
        "-M", "--spm-mcr", type=is_directory, help="Path to SPM MCR directory")
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
spm_version = spm_get_version(inputs["spm_sh"], inputs["spm_mcr"])
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "rsfmri_preproc_dartel_part2_existing_template",
    "spm_version": spm_version
}
outputs = {}
if verbose > 0:
    pprint("[info] Starting rsfMRI preprocessing from existing template...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Run dartel
flow_field, m_file = run_dartel_from_existing_template(
    templates=inputs["templates"],
    rc1_file=inputs["gm_seg"],
    rc2_file=inputs["wm_seg"],
    outdir=inputs["outdir"],
    spm_sh=inputs["spm_sh"],
    spm_mcr=inputs["spm_mcr"],
    inner_its=[3, 3, 3, 3, 3, 3],
    reg_form="Linear Elastic Energy",
    timesteps=[1, 1, 2, 4, 16, 64],
    reg_params=[["4", "2", "1e-06"], ["2", "1", "1e-06"],
                ["1", "0.5", "1e-06"], ["0.5", "0.25", "1e-06"],
                ["0.25", "0.125", "1e-06"], ["0.25", "0.125", "1e-06"]],
    lm_reg=0.01,
    nb_cycles=3,
    its=3)
outputs["Flow field"] = flow_field
outputs["Matlab script"] = m_file


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
script_name = "rsfmri_preproc_dartel_part2_existing_template"
for name, final_struct in [("{0}_inputs".format(script_name), inputs),
                           ("{0}_outputs".format(script_name), outputs),
                           ("{0}_runtime".format(script_name), runtime)]:
    log_file = os.path.join(logdir, "{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 1:
    print("[info] Outputs:")
    pprint(outputs)
