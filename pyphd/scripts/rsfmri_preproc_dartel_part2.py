# Script reproducing Clara Manesco matlab script for rsfMRI data preprocessing
# on the MAPT cohort with Dartel registration.
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

# Third party imports
from traits import trait_base
import nipype
from nipype.interfaces import spm

# Script documentation
DOC = """
Execute Clara Manesco fMRI preprocessing pipeline initially designed to be run
on MAPT data with DARTEL registration.
------------------------------------------------------------------------------

/!\ This script is part of a group of three scripts
    rsfmri_preproc_dartel_part*.py to execute preprocessing with dartel
    registration /!\

Create the DARTEL template from a list of subjects segmentation files.
rsfmri_preproc_dartel_part1.py must be run first.

/!\ IMPORTANT /!\ 
Run DARTEL (existing template) is not available through nipype.
If template1...6.nii has already been created, go through the regular SPM
interface to register rc*.nii files.

Steps:
0) Set environment and copy all DARTEL GM/WM seg files in a single directory
   (nipype oblige...)
1) Generate DARTEL template.

Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/rsfmri_preproc_dartel_part2.py \
    -g /tmp/test_fmri/sub-03990425CHE/ses-M0/rc1sub-03990425CHE_ses-M0_T1w_orig_ac.nii \
       /tmp/test_fmri/sub-03990171AGE/ses-M0/rc1sub-03990171AGE_ses-M0_T1w_orig_ac.nii \
    -w /tmp/test_fmri/sub-03990425CHE/ses-M0/rc2sub-03990425CHE_ses-M0_T1w_orig_ac.nii \
       /tmp/test_fmri/sub-03990171AGE/ses-M0/rc2sub-03990171AGE_ses-M0_T1w_orig_ac.nii \
    -o /tmp/test_fmri \
    -S $SPM12_DIR/spm12/run_spm12.sh \
    -M $SPM12_DIR/mcr/v713/ \
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
        prog="python rsfmri_preproc_dartel_part2.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-g", "--gm-seg", type=is_file, required=True, nargs="+",
        help="Path to subjects GM segmentations FOR DARTEL "
             "(look at rc1 files).")
    required.add_argument(
        "-w", "--wm-seg", type=is_file, required=True, nargs="+",
        help="Path to subjects WM segmentations FOR DARTEL "
             "(look at rc2 files). Must be in the same order as GM seg files.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-C", "--center", type=str, help="Center name.")
    parser.add_argument(
        "-T", "--timepoint", type=str, help="Timepoint.")
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
    "tool": "rsfmri_preproc_dartel_part2",
    "nipype_version": nipype.__version__,
    "nipype_file": nipype.__file__,
    "spm_version": spm_version
}
outputs = {}
if verbose > 0:
    pprint("[info] Starting rsfMRI preprocessing...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Step 0 - Set up SPM environment
"""

# Force to use SPM standalone
matlab_cmd = '{0} {1} script'.format(inputs["spm_sh"], inputs["spm_mcr"])
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

# Create dartel directory
dartel_dir = os.path.join(inputs["outdir"], "dartel")
if not os.path.isdir(dartel_dir):
    os.mkdir(dartel_dir)

# Copy all files to dartel directory
# TODO : find a way to specify moving files
gm_segs = []
wm_segs = []
to_delete = True  # Delete files after copy
for idx, gm_seg in enumerate(inputs["gm_seg"]):
    wm_seg = inputs["wm_seg"][idx]
    if os.path.dirname(gm_seg) == dartel_dir:
        to_delete = False
    if os.path.dirname(wm_seg) == dartel_dir:
        to_delete = False
    shutil.copy(gm_seg, dartel_dir)
    gm_segs.append(os.path.join(dartel_dir, os.path.basename(gm_seg)))
    shutil.copy(wm_seg, dartel_dir)
    wm_segs.append(os.path.join(dartel_dir, os.path.basename(wm_seg)))
os.chdir(dartel_dir)


"""
Step 1 - Create dartel template
"""

# Note : c1=gm, c2=wm, c3=csf. Here we only care for gm and wm segmentations.
dartel = spm.DARTEL()
dartel.inputs.image_files = [gm_segs, wm_segs]
dartel_results = dartel.run()
dartel_results = dartel_results.outputs
outputs["Dartel flow fields"] = dartel_results.dartel_flow_fields
outputs["Final template file"] = dartel_results.final_template_file
outputs["Intermediate template files"] = dartel_results.template_files


# Delete all files copied
if to_delete:
    to_delete_files = gm_segs + wm_segs
    for fid in to_delete_files:
        print("Deleting {0}...".format(fid))
        os.remove(fid)


"""
Sanitize outputs
"""
for key, data in outputs.items():
    if type(data) == type(trait_base.Undefined):
        outputs[key] = None


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
script_name = "rsfmri_preproc_dartel_part2"
if inputs["center"] is not None:
    script_name += "_" + inputs["center"]
if inputs["timepoint"] is not None:
    script_name += "_" + inputs["timepoint"]
for name, final_struct in [("{0}_inputs".format(script_name), inputs),
                           ("{0}_outputs".format(script_name), outputs),
                           ("{0}_runtime".format(script_name), runtime)]:
    log_file = os.path.join(logdir, "{0}_{1}.json".format(
        name, inputs["timepoint"]))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 1:
    print("[info] Outputs:")
    pprint(outputs)
