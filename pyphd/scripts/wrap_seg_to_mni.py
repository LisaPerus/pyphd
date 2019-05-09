# Script to wrap c1, c2, c3 segmentation to MNI
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
import nipype
from nipype.interfaces import spm


# Script documentation
DOC = """
Wraps c1, c2, c3 outputs of SPM segmentation step to MNI, provided a
deformation field is given.
------------------------------------------------------------------------------

Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/wrap_seg_to_mni.py \
    -c /tmp/test_warp/c1sub-03990366TRE_ses-M0_T1w_orig_ac.nii \
       /tmp/test_warp/c2sub-03990366TRE_ses-M0_T1w_orig_ac.nii \
    -y /tmp/test_warp/y_sub-03990366TRE_ses-M0_T1w_orig_ac.nii \
    -o /tmp/test_warp/ \
    -S $HOME/opt/spm12/spm12/run_spm12.sh \
    -M $HOME/opt/spm12/mcr/v713/ \
    -V 1
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
        prog="python wrap_seg_to_mni.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-c", "--seg-im", type=is_file, nargs='+', required=True,
        help="Path to segmentation outputs.")
    required.add_argument(
        "-y", "--deformation-field", type=is_file, required=True,
        help="Deformation field.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-F", "--m-file-dir", type=str, default="m_files",
        help="Additional directory at the root of output where to store m "
             "file.")
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
    "tool": "wrap_seg_to_mni.py",
    "nipype_version": nipype.__version__,
    "nipype_file": nipype.__file__,
    "spm_version": spm_version
}
outputs = {}
if verbose > 0:
    pprint("[info] Wrapping segmentation to MNI...")
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

# Go to output directory
os.chdir(inputs["outdir"])


"""
Step 1 - Use SPM normalize12 to wrap segmentation
"""
norm12 = spm.Normalize12()
norm12.inputs.deformation_file = inputs["deformation_field"]
norm12.inputs.apply_to_files = inputs["seg_im"]
norm12.inputs.jobtype = 'write'
norm12_results = norm12.run()
wrapped_files = norm12_results.outputs.normalized_files
outputs["Wrapped seg files"] = wrapped_files


"""
Step 2 - Move m file to outdir sub-directory
"""
if inputs["m_file_dir"] is not None:
    m_file = glob.glob(
        os.path.join(inputs["outdir"], "pyscript_normalize12.m"))
    if len(m_file) == 0:
        raise ValueError("No m file generated by nipype SPM normalize12.")
    m_file = m_file[0]
    m_file_dir = os.path.join(inputs["outdir"], inputs["m_file_dir"])
    if not os.path.isdir(m_file_dir):
        os.mkdir(m_file_dir)
    m_move_file = os.path.join(m_file_dir, "pyscript_seg_normalize12.m")
    shutil.move(m_file, m_move_file)

    # Also remove pyscript.m
    pyscript_m_file = glob.glob(
        os.path.join(inputs["outdir"], "pyscript.m"))
    if len(pyscript_m_file) == 0:
        raise ValueError(
            "No pyscript.m file generated by nipype SPM normalize12.")
    pyscript_m_file = pyscript_m_file[0]
    print("Removing {0}...".format(pyscript_m_file))
    os.remove(pyscript_m_file)


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(os.path.join(inputs["outdir"], "logs"))
if not os.path.isdir(logdir):
    os.mkdir(logdir)
script_name = "wrap_seg_to_mni"
for name, final_struct in [("{0}_inputs".format(script_name), inputs),
                           ("{0}_outputs".format(script_name), outputs),
                           ("{0}_runtime".format(script_name), runtime)]:

    log_file = os.path.join(logdir, "{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[info] Outputs:")
    pprint(outputs)
