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
import progressbar
from datetime import datetime
from argparse import RawTextHelpFormatter
from datetime import datetime
from pprint import pprint

# Third party imports
from pyfreesurfer.wrapper import FSWrapper

# Script documentation
DOC = """
Get cortical thickness for one subject and one ROI
--------------------------------------------------

Script to get cortical thickness from a subject anatomical data processed
with Freesurfer (recon-all) over a ROI volume projected onto fsaverage surface.

Before running this script make sure that:
1) The subject has been processed with recon-all
2) You have a transformation (.dat) from your ROI to the fsaverage subject
space
3) You have mapped your ROI volume to the fsaverage surface
4) You have mapped your subject thickness to the fsaverage subject

For 2) and 3) see script map_nifti_roi_to_freesurfer_fsaverage_surface.py
For 4) see script map_freesurfer_subject_thickness_to_fsaverage.py

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/get_freesurfer_subject_roi_thickness.py \
    -s 03990171AGE_M0 \
    -d $NFS_CATI/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/ \
    -f $HOME/setup_freesurfer_i2bm.sh \
    -l /nfs/neurospin/cati/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/roi_to_fsaverage_surface/lh.fsaverage.resliced_to_fsl6.0.0_MNI152_T1_2mm_braintemplate_cambridge_basc_multiscale_asym_scale036_cluster_35.0.mgh \
    -r /nfs/neurospin/cati/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/roi_to_fsaverage_surface/rh.fsaverage.resliced_to_fsl6.0.0_MNI152_T1_2mm_braintemplate_cambridge_basc_multiscale_asym_scale036_cluster_35.0.mgh \
    -n ct_cambridge_basc_multiscale_asym_scale036_cluster_35 \
    -o /tmp/test \
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
        prog="python map_nifti_roi_to_freesurfer_fsaverage_surface.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-s", "--sid", type=str, required=True,
        help="Subject ID.")
    required.add_argument(
        "-d", "--fs-sub-dir", type=is_directory,
        help="Path to Freesurfer subjects directory.")
    required.add_argument(
        "-f", "--fs-sh", type=is_file, help="Path to Freesurfer sh file.")
    required.add_argument(
        "-l", "--lh-fsaverage-roi-surface-overlay", type=is_file, nargs="+",
        help="Path to Freesurfer fsaverage ROI overlay files for left hemi.")
    required.add_argument(
        "-r", "--rh-fsaverage-roi-surface-overlay", type=is_file, nargs="+",
        help="Path to Freesurfer fsaverage ROI overlay files for right hemi.")
    required.add_argument(
        "-n", "--output-basename", type=str, nargs="+",
        help="Output files basenames. This basename will be added after sid "
             "and left/right hemi.")
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
outputs = {}
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "get_freesurfer_subject_roi_thickness.py"
}

# Set freesurfer subjects directory path
os.environ["SUBJECTS_DIR"] = inputs["fs_sub_dir"]

# Go to subject directory
sid_surf_dir = os.path.join(inputs["fs_sub_dir"], inputs["sid"], "surf")
if not os.path.isdir(sid_surf_dir):
    raise ValueError(
        "Cannot find surf directory for subject {0}".format(inputs["sid"]))
os.chdir(sid_surf_dir)

# Run segstats to get subject thickness on ROI
for idx_roi, roi_name in enumerate(inputs["output_basename"]):
    for side in ["lh", "rh"]:
        outfile = os.path.join(
            inputs["outdir"],
            inputs["sid"] + ".{0}.".format(side) + roi_name + ".txt")
        cmd = ["mri_segstats", "--seg",
               inputs["{0}_fsaverage_roi_surface_overlay".format(side)][
                idx_roi],
               "--in", "{0}.thickness.fsaverage.mgh".format(side),
               "--sum", outfile]
        fscmd = FSWrapper(cmd, inputs["fs_sh"])
        fscmd()
        if not os.path.isfile(outfile):
            raise ValueError(
                "Cannot find stat summary file {0}...".format(outfile))
        outputs["{0} {1} stats".format(roi_name, side)] = outfile


"""
Write outputs
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(
        logdir,
        "get_freesurfer_subject_roi_thickness_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)

if verbose > 0:
    pprint("[Outputs]:")
    pprint(outputs)
