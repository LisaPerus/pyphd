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
Map a subject freesurfer cortical thickness to fsaverage
--------------------------------------------------------

/!\ Freesurfer must be activated in shell first /!\

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/map_freesurfer_subject_thickness_to_fsaverage.py \
    -s 03990171AGE_M36 \
    -d $NFS_CATI/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/03990171AGE_M36/surf \
    -f ~/setup_freesurfer_i2bm.sh \
    -i $NFS_CATI/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03 \
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
        prog="python map_freesurfer_subject_thickness_to_fsaverage.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-s", "--sid", type=str, required=True,
        help="Subject ID.")
    required.add_argument(
        "-d", "--sid-surf-dir", type=is_directory, required=True,
        help="Path to freesurfer subject surf directory.")
    required.add_argument(
        "-f", "--fs-sh", type=is_file, help="Path to Freesurfer sh file.")
    required.add_argument(
        "-i", "--fs-sub-dir", type=is_directory,
        help="Path to Freesurfer subjects directory.")

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
    "tool": "map_freesurfer_subject_thickness_to_fsaverage.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Map subject {0} thickness to fsaverage...".format(
           inputs["sid"]))
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Set freesurfer subjects directory path
os.environ["SUBJECTS_DIR"] = inputs["fs_sub_dir"]

# Go to subject directory
os.chdir(inputs["sid_surf_dir"])

# Run mapping
# > For left hemisphere
cmd = ["mri_surf2surf", "--s", inputs["sid"], "--trgsubject", "fsaverage",
       "--hemi", "lh", "--sval", "lh.thickness", "--tval",
       "lh.thickness.fsaverage.mgh"]
fscmd = FSWrapper(cmd, inputs["fs_sh"])
fscmd()

# > For right hemisphere
cmd = ["mri_surf2surf", "--s", inputs["sid"], "--trgsubject", "fsaverage",
       "--hemi", "rh", "--sval", "rh.thickness", "--tval",
       "rh.thickness.fsaverage.mgh"]
fscmd = FSWrapper(cmd, inputs["fs_sh"])
fscmd()

# Check that output has been created
output_lh_thick_fsaverage = os.path.join(
    inputs["sid_surf_dir"], "lh.thickness.fsaverage.mgh")
output_rh_thick_fsaverage = os.path.join(
    inputs["sid_surf_dir"], "rh.thickness.fsaverage.mgh")
if not os.path.isfile(output_lh_thick_fsaverage):
    raise ValueError(
        "Cannot find subject {0} lh thickness on fsaverage".format(
            inputs["sid"]))
outputs["lh thickness fsaverage"] = output_lh_thick_fsaverage
outputs["rh thickness fsaverage"] = output_rh_thick_fsaverage


"""
Write outputs
"""
logdir = os.path.join(inputs["sid_surf_dir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(
        logdir,
        "map_freesurfer_subject_thickness_to_fsaverage_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)

if verbose > 0:
    pprint("[Outputs]:")
    pprint(outputs)
