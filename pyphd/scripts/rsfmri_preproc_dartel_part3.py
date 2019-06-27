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
from collections import OrderedDict

# Third party imports
from traits import trait_base
import nipype
from nipype.interfaces import spm

# PyPHD imports
from pyphd.constants import RSFMRI_PREPROC_CLARA_MANESCO
from pyphd.spm.utils import dartel_normalize_to_mni

# Script documentation
DOC = """
Execute Clara Manesco fMRI preprocessing pipeline initially designed to be run
on MAPT data with DARTEL registration.
------------------------------------------------------------------------------

/!\ This script is part of a group of three scripts
    rsfmri_preproc_dartel_part*.py to execute preprocessing with dartel
    registration /!\

Steps:
1) Normalize subject anatomical and functional data from DARTEL generated
template.
2) Normalize subject functional data from DARTEL generated template
with minimal smoothing (fwhm at 3)
3) Normalize subject functional data from DARTEL generated template WITH
smoothing (Smoothing is included in the normalization function).
4) Normalize segmentations files (wm, gm, csf only).

Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/rsfmri_preproc_dartel_part3.py \
    -s sub-03990171AGE \
    -t ses-M0 \
    -f /tmp/test_fmri/sub-03990171AGE/ses-M0/rasub-03990171AGE_ses-M0_task-rest_bold_orig_ac.nii \
    -a /tmp/test_fmri/sub-03990171AGE/ses-M0/sub-03990171AGE_ses-M0_T1w_orig_ac.nii \
    -d /tmp/test_fmri/dartel/Template_6.nii \
    -w /tmp/test_fmri/dartel/u_rc1sub-03990171AGE_ses-M0_T1w_orig_ac_Template.nii \
    -o /tmp/test_fmri/ \
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
        prog="python rsfmri_preproc_dartel_part3.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-s", "--sid", type=str, required=True,
        help="Subject ID.")
    required.add_argument(
        "-t", "--timepoint", type=str, required=True,
        help="Timepoint.")
    required.add_argument(
        "-f", "--func-im", type=is_file, required=True,
        help="Path to functional image already preprocessed by "
             "rsfmri_preproc_dartel_part1.py.")
    required.add_argument(
        "-a", "--anat-im", type=is_file, required=True,
        help="Path to anatomical image already preprocessed by "
             "rsfmri_preproc_dartel_part1.py.")
    required.add_argument(
        "-e", "--seg-im", type=is_file, required=True, nargs="+",
        help="Path to segmentation image c1*.nii, c2*.nii and c3*.nii.")
    required.add_argument(
        "-d", "--dartel-template", type=is_file, required=True,
        help="Path to dartel template.")
    required.add_argument(
        "-w", "--flow-field", type=is_file, required=True,
        help="Path to flow field file.")
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
    "tool": "rsfmri_preproc_dartel_part3",
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

# Create outdir + go to directory
subdir = os.path.join(inputs["outdir"], inputs["sid"], inputs["timepoint"])

# > If directory does not exist, create it
if not os.path.isdir(subdir):
    os.makedirs(subdir)

# > Go to directory
os.chdir(subdir)

# > Copy Template to subject directory (later delete copy)
# Else nipype cannot find its output spec.
delete_template = False
if os.path.dirname(inputs["dartel_template"]) != subdir:
    copy_template = os.path.join(
        subdir, os.path.basename(inputs["dartel_template"]))
    shutil.copy(inputs["dartel_template"], subdir)
    delete_template = True
    template = copy_template
else:
    template = inputs["dartel_template"]

# Get parameters
parameters = RSFMRI_PREPROC_CLARA_MANESCO
np_parameters = parameters["DartelNormalize2MNI"]


"""
Step 1 : Normalize subject anatomical data
"""
print("Normalize anatomical data...")
registered_input_file, template_mat_file, m_file = dartel_normalize_to_mni(
    template=template,
    flowfield=inputs["flow_field"],
    infile=inputs["anat_im"],
    fwhm=np_parameters["fwhm"]["anat"],
    voxel_size=np_parameters["voxel_size"],
    outdir=subdir,
    spm_sh=inputs["spm_sh"],
    spm_mcr=inputs["spm_mcr"],
    bb=np_parameters["bounding_box"],
    modulate=np_parameters["modulate"]
)
outputs["Anat normalization file"] = template_mat_file
outputs["Normalized anat file"] = registered_input_file


"""
Step 2 : Normalize subject functional data with minimal smoothing
"""
print("Normalize functional data with minimal smoothing...")
registered_input_file, template_mat_file, m_file = dartel_normalize_to_mni(
    template=template,
    flowfield=inputs["flow_field"],
    infile=inputs["func_im"],
    fwhm=[3., 3., 3.],
    voxel_size=np_parameters["voxel_size"],
    outdir=subdir,
    spm_sh=inputs["spm_sh"],
    spm_mcr=inputs["spm_mcr"],
    bb=np_parameters["bounding_box"],
    modulate=np_parameters["modulate"]
)

# Change output name
new_registered_input_file = registered_input_file.replace(
    ".nii", "_fwhm3iso.nii")
shutil.move(registered_input_file, new_registered_input_file)
outputs["Func with minimal smoothing normalization file"] = template_mat_file
outputs["Normalized minimal smoothing func file"] = new_registered_input_file


"""
Step 3 : Normalize subject functional data with smoothing
"""
print("Normalize functional data with smoothing...")
registered_input_file, template_mat_file, m_file = dartel_normalize_to_mni(
    template=template,
    flowfield=inputs["flow_field"],
    infile=inputs["func_im"],
    fwhm=np_parameters["fwhm"]["func"],
    voxel_size=np_parameters["voxel_size"],
    outdir=subdir,
    spm_sh=inputs["spm_sh"],
    spm_mcr=inputs["spm_mcr"],
    bb=np_parameters["bounding_box"],
    modulate=np_parameters["modulate"]
)
outputs["Func with smoothing normalization file"] = template_mat_file
outputs["Normalized smoothing func file"] = registered_input_file


"""
Step 4 : Normalize segmentation
"""
print("Normalize segmentation data...")
for seg in inputs["seg_im"]:
    registered_input_file, template_mat_file, m_file = dartel_normalize_to_mni(
        template=template,
        flowfield=inputs["flow_field"],
        infile=seg,
        fwhm=[0.0, 0.0, 0.0],
        voxel_size=np_parameters["voxel_size"],
        outdir=subdir,
        spm_sh=inputs["spm_sh"],
        spm_mcr=inputs["spm_mcr"],
        bb=np_parameters["bounding_box"],
        modulate=np_parameters["modulate"]
    )
    basename_seg = ""
    if "c1" in os.path.basename(seg):
        basename_seg = "GM"
    if "c2" in os.path.basename(seg):
        basename_seg = "WM"
    if "c3" in os.path.basename(seg):
        basename_seg = "CSF"

    outputs["{0} Seg normalization file".format(
        basename_seg)] = template_mat_file
    outputs["Normalized {0} seg file".format(
        basename_seg)] = registered_input_file


# If template has been copied to subdir, delete it
if delete_template:
    os.remove(copy_template)

"""
Sanitize outputs
"""
for key, data in outputs.items():
    if type(data) == type(trait_base.Undefined):
        outputs[key] = None


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(subdir, "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
script_name = "rsfmri_preproc_dartel_part3"
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
