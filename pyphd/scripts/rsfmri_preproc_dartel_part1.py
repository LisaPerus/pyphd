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
import nibabel
from traits import trait_base
import nipype
from nipype.interfaces import spm

# Pypreprocess + own code imports
from pypreprocess import __version__ as pypreprocess_version
from pypreprocess.io_utils import niigz2nii

# PyPHD imports
from pyphd.constants import RSFMRI_PREPROC_CLARA_MANESCO
from pyphd.spm.utils import (spm_standalone_reorient, get_slice_order,
                             st_get_ref_slice)


# Script documentation
DOC = """
Execute Clara Manesco fMRI preprocessing pipeline initially designed to be run
on MAPT data with DARTEL registration.
------------------------------------------------------------------------------

Note : it is controversial whether to do slice timing correction after or
before motion correction. As nothing clear was specified from SPM, it was
chosen to do stc first.

/!\ This script is part of a group of three scripts
    rsfmri_preproc_dartel_part*.py to execute preprocessing with dartel
    registration /!\

Executes the preprocessing pipeline as follows:
1) Reorient functional images by centering on AC.
2) Correct slice timing for functional image.
3) Correct functional image for motion.
4) Coregistration of functional to anatomical images.
5) Segment anatomical images.


Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/rsfmri_preproc_dartel_part1.py \
    -s sub-03990171AGE \
    -t ses-M0 \
    -f /home/lp259104/PHD/DATA/MAPT/0399/sub-03990171AGE/ses-M0/func/sub-03990171AGE_ses-M0_task-rest_bold.nii.gz \
    -a /home/lp259104/PHD/DATA/MAPT/0399/sub-03990171AGE/ses-M0/anat/sub-03990171AGE_ses-M0_T1w.nii.gz \
    -l "Ascending interleaved" \
    -n SIEMENS \
    -r Middle \
    -o /tmp/test_fmri/ \
    -C 4.6 34.9 13.2 \
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
        prog="python rsfmri_preproc_dartel_part1.py",
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
        help="Path to functional image.")
    required.add_argument(
        "-a", "--anat-im", type=is_file, required=True,
        help="Path to anatomical image.")
    required.add_argument(
        "-l", "--slice-order", type=str, required=True,
        help="Image slice order.")
    required.add_argument(
        "-r", "--st-ref-slice", type=str, required=True,
        help="Slice timing reference slice. Can be First, Middle or Last.")
    required.add_argument(
        "-n", "--scanner", type=str, required=True,
        help="Scanner.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-C", "--commissure", type=float, nargs="+",
        help="Commissure coordinates in scanner space (mm).")
    parser.add_argument(
        "-E", "--clean", action="store_true",
        help="Clean subject directory if it already exists.")
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
    "tool": "rsfmri_preproc_dartel_part1",
    "nipype_version": nipype.__version__,
    "nipype_file": nipype.__file__,
    "spm_version": spm_version,
    "pypreprocess_version": pypreprocess_version
}
outputs = {}
if verbose > 0:
    pprint("[info] Starting rsfMRI preprocessing...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Step 0 -
1) Set up SPM environment
2) Create outdir + go to outdir (so that m file are not created anywhere)
TODO: find if a path can be specified for the creation of the intermediate
      m files. Have a look to pypreprocess.
3) Load parameters
4) Transform .nii.gz to .nii
"""

# Force to use SPM standalone
matlab_cmd = '{0} {1} script'.format(inputs["spm_sh"], inputs["spm_mcr"])
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

# Create outdir + go to directory
subdir = os.path.join(inputs["outdir"], inputs["sid"], inputs["timepoint"])

# > Clean directory if it exists and if option clean is set to True
if os.path.isdir(subdir) and inputs["clean"]:
    shutil.rmtree(subdir)

# > If directory does not exist, create it
if not os.path.isdir(subdir):
    os.makedirs(subdir)

# > Go to directory
os.chdir(subdir)

# Get parameters
parameters = RSFMRI_PREPROC_CLARA_MANESCO

# Unzip nii.gz to nii
func_im = niigz2nii(inputs["func_im"], output_dir=subdir)
anat_im = niigz2nii(inputs["anat_im"], output_dir=subdir)


"""
Step 1 - Reorient subject by setting origin to Anterior Commissure coordinates
"""
print("Reorienting to anterior commissure...")
if inputs["commissure"] is not None:

    # Copy original images
    reoriented_func = os.path.join(
        subdir, "{0}_{1}_task-rest_bold_orig_ac.nii".format(
            inputs["sid"], inputs["timepoint"]))
    shutil.copy(func_im, reoriented_func)
    reoriented_anat = os.path.join(
        subdir, "{0}_{1}_T1w_orig_ac.nii".format(
            inputs["sid"], inputs["timepoint"]))
    shutil.copy(anat_im, reoriented_anat)

    # Reorient func
    print(
        "Centering on {0} for func and anat images...".format(
            " ".join([str(x) for x in inputs["commissure"]])))
    func_im = spm_standalone_reorient(
        im=reoriented_func,
        sid=inputs["sid"],
        origin_coords=inputs["commissure"],
        outdir=subdir,
        spm_sh=inputs["spm_sh"],
        spm_mcr=inputs["spm_mcr"],
        delete_mfile=False,
        delete_mat_file=True)
    anat_im = spm_standalone_reorient(
        im=reoriented_anat,
        sid=inputs["sid"],
        origin_coords=inputs["commissure"],
        outdir=subdir,
        spm_sh=inputs["spm_sh"],
        spm_mcr=inputs["spm_mcr"],
        delete_mfile=False,
        delete_mat_file=True)
else:
    reoriented_func = None
    reoriented_anat = None
outputs["Reoriented images"] = [reoriented_func, reoriented_anat]


"""
Step 2 - Slice Timing Correction
"""

# Get number of slices in volumes
# > Compute ta and slice order
st_parameters = parameters["SliceTiming"]
im = nibabel.load(func_im)
nslices = im.shape[2]
tr = im.header.get_zooms()[3]
if tr == 0 or tr == 1:
    raise ValueError("Suspicious TR value in NIFTI header : {0}".format(tr))
ta = tr - (tr / nslices)
slice_order_list = get_slice_order(
    nslices, inputs["slice_order"], inputs["scanner"])

# > Compute slice of reference index
ref_slice_idx = st_get_ref_slice(
    ref_slice=inputs["st_ref_slice"],
    slice_order=slice_order_list)

# Run Slice timing correction
print("Slice timing correction...")
st = spm.SliceTiming()
st.inputs.in_files = func_im
st.inputs.num_slices = nslices
st.inputs.time_repetition = tr
st.inputs.time_acquisition = ta
st.inputs.slice_order = slice_order_list
st.inputs.ref_slice = ref_slice_idx
st.inputs.out_prefix = st_parameters["out_prefix"]
st_results = st.run()
timecorrected_func_im = st_results.outputs.timecorrected_files
outputs["TA"] = ta
outputs["Slice order"] = slice_order_list
outputs["Ref slice"] = ref_slice_idx
outputs["Slice timing corrected func image"] = timecorrected_func_im


"""
Step 3 - Motion correction (Realignement)
"""
print("Func data motion correction...")
motion_corr_parameters = parameters["Realign"]
rl = spm.Realign()
rl.inputs.in_files = timecorrected_func_im
rl.inputs.quality = motion_corr_parameters["quality"]
rl.inputs.separation = motion_corr_parameters["separation"]
rl.inputs.fwhm = motion_corr_parameters["fwhm"]
rl.inputs.register_to_mean = motion_corr_parameters["register_to_mean"]
rl.inputs.interp = motion_corr_parameters["interp"]
rl.inputs.wrap = motion_corr_parameters["wrap"]
if motion_corr_parameters["weight_img"] is not None:
    rl.inputs.weight_img = motion_corr_parameters["weight_img"]
rl.inputs.write_which = motion_corr_parameters["write_which"]
rl.inputs.write_interp = motion_corr_parameters["write_interp"]
rl.inputs.write_wrap = motion_corr_parameters["write_wrap"]
rl.inputs.write_mask = motion_corr_parameters["write_mask"]
rl.inputs.out_prefix = motion_corr_parameters["out_prefix"]
rl.inputs.jobtype = motion_corr_parameters["jobtype"]
realign_results = rl.run()
motion_corr_func_im = realign_results.outputs.realigned_files
realignment_parameters = realign_results.outputs.realignment_parameters
mean_image = realign_results.outputs.mean_image
modified_in_file = realign_results.outputs.modified_in_files
outputs["St + realigned func images"] = motion_corr_func_im
outputs["Realignement motion parameters"] = realignment_parameters


"""
Step 4 - Coregistration with T1
"""

# Coregister func and anat - estimate only
print("Func and anat coregistration estimation...")
coreg_parameters = parameters["Coregister"]
coreg = spm.Coregister()
coreg.inputs.target = anat_im
coreg.inputs.source = mean_image
coreg.inputs.apply_to_files = motion_corr_func_im
coreg.inputs.cost_function = coreg_parameters["cost_function"]
coreg.inputs.separation = coreg_parameters["separation"]
coreg.inputs.tolerance = coreg_parameters["tolerance"]
coreg.inputs.fwhm = coreg_parameters["fwhm"]
coreg.inputs.jobtype = "estimate"
coreg_results = coreg.run()
coregistered_estimate_source = coreg_results.outputs.coregistered_source
coregistered_estimate_files = coreg_results.outputs.coregistered_files

# Save coregistered mean func to anat image, to check later if
# coregistration worked
coreg = spm.Coregister()
coreg.inputs.target = anat_im
coreg.inputs.source = coregistered_estimate_source
coreg.inputs.write_interp = coreg_parameters["write_interp"]
coreg.inputs.write_mask = coreg_parameters["write_mask"]
coreg.inputs.write_wrap = coreg_parameters["write_wrap"]
coreg.inputs.jobtype = "write"
coreg_results = coreg.run()
mean_func_coregistered = coreg_results.outputs.coregistered_source

"""
Step 5 - Segment with Dartel
"""
print("Segment anatomical data...")
seg_parameters = parameters["NewSegmentDartel"]
seg = spm.NewSegment()

# Json file cannot contain tuple type : change array to tuples when necessary
tissues = []
for tissue in seg_parameters["tissues"]:
    tuple_tissue = []
    for idx, elt in enumerate(tissue):
        if type(elt) == list:
            tuple_tissue.append(tuple(elt))
        else:
            tuple_tissue.append(elt)
    tissues.append(tuple(tuple_tissue))
channel_info_tuple = []
for idx, info in enumerate(seg_parameters["channel_info"]):
    if type(info) == list:
        channel_info_tuple.append(tuple(info))
    else:
        channel_info_tuple.append(info)
channel_info_tuple = tuple(channel_info_tuple)
seg.inputs.channel_files = anat_im
seg.inputs.tissues = tissues
seg.inputs.affine_regularization = seg_parameters["affine_regularization"]
seg.inputs.sampling_distance = seg_parameters["sampling_distance"]
seg.inputs.warping_regularization = seg_parameters["warping_regularization"]
seg.inputs.channel_info = channel_info_tuple
seg.inputs.use_mcr = True

# Run segmentation
seg_results = seg.run()
seg_results = seg_results.outputs
outputs["Bias corrected image"] = seg_results.bias_corrected_images
outputs["Bias field images"] = seg_results.bias_field_images
outputs["Dartel inputs image"] = seg_results.dartel_input_images
outputs["Forward deformation field"] = seg_results.forward_deformation_field
outputs["Inverse Deformation field"] = seg_results.inverse_deformation_field
outputs["Modulated probability maps images"] = (
    seg_results.modulated_class_images)
outputs["Native probability maps images"] = seg_results.native_class_images
outputs["Normalized probability maps"] = seg_results.normalized_class_images
outputs["Transformation matrix"] = seg_results.transformation_mat


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
script_name = "rsfmri_preproc_dartel_part1"
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
