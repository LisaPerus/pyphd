#!/usr/bin/env python
# Script to compare reho map evolution at two timepoints
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
import shutil

# Third party imports
from nipype import Node, Workflow
from nipype.interfaces import spm


# Script documentation
DOC = """Design t-tests from reho map data.
-------------------------------------------

1) Create test
2) Estimate model

Example on MAPT data:
python3 $SCRIPT_DIR/GIT_REPOST/pyphd/pyphd/design_test_reho_maps.py \
    -i /home/lp259104/PHD/SCRIPTS/fMRI_ANALYSIS/MAPT_reho_maps.txt \
    -o /tmp/test_second_level \
    -S $SPM12_DIR/spm12/run_spm12.sh \
    -M $SPM12_DIR/mcr/v713/ \
    -V 1

python3 $SCRIPT_DIR/GIT_REPOST/pyphd/pyphd/design_test_reho_maps.py \
    -i /home/lp259104/PHD/DATA/RESULTS/MAPT/REHO_ANALYSIS/Clara_n118/MAPT_reho_maps.txt \
    -o /tmp/test_second_level/PairedTTest \
    -a PairedTTestDesign \
    -S $SPM12_DIR/spm12/run_spm12.sh \
    -M $SPM12_DIR/mcr/v713/ \
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
        "-i", "--data-file", type=is_file, required=True, metavar="<path>",
        help="Path to input file listing reho maps at two timepoints."
             "Must be comma separated with a one-line header.")
    required.add_argument(
        "-a", "--analysis", type=str, required=True, metavar="<path>",
        help="Type of analysis, can be : OneSampleTTestDesign or "
             "PairedTTestDesign.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-T", "--timepoint", type=int, metavar="<path>",
        default=0, help="If design is a OneSampleTTestDesign, select column "
        "in data-file to use as sample.")
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
runtime = {
    "tool": "evolution_reho_maps.py",
    "timestamp": datetime.now().isoformat()}
outputs = {}

if verbose > 0:
    pprint("[info] Compare evolution of reho maps...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Get to outdir
os.chdir(inputs["outdir"])

# Set up spm env
# > Force to use SPM standalone
matlab_cmd = '{0} {1} script'.format(inputs["spm_sh"], inputs["spm_mcr"])
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

# Read data
timepoints_data = []
with open(inputs["data_file"], "rt") as open_file:
    lines = open_file.readlines()

# Copy files to outdir
imgs_per_timepoints = [[], []]
copied_files = []
for line in lines[1:]:
    line = line.strip("\n").strip(" ").split(",")
    shutil.copy(line[0].strip("\n"), inputs["outdir"])
    imgs_per_timepoints[0].append(
        os.path.join(inputs["outdir"], os.path.basename(line[0].strip("\n"))))
    if inputs["analysis"] == "OneSampleTTestDesign":
        continue
    shutil.copy(line[1].strip("\n"), inputs["outdir"])
    copied_file = os.path.join(
        inputs["outdir"], os.path.basename(line[1].strip("\n")))
    copied_files.append(copied_file)
    imgs_per_timepoints[1].append(
        os.path.join(inputs["outdir"], os.path.basename(line[1].strip("\n"))))
outputs["Copied files "] = copied_files

# Create workflow
analysis2nd = Workflow(
    name='work_2nd', base_dir=inputs["outdir"])

# Create Design
if inputs["analysis"] == "OneSampleTTestDesign":
    stat_test = Node(spm.model.OneSampleTTestDesign(), name="stattest")
    stat_test.inputs.in_files = imgs_per_timepoints[inputs["timepoint"]]
elif inputs["analysis"] == "PairedTTestDesign":
    stat_test = Node(spm.model.PairedTTestDesign(), name="stattest")
    paired_ims = []
    for idx, elt in enumerate(imgs_per_timepoints[0]):
        paired_ims.append(
            [imgs_per_timepoints[0][idx], imgs_per_timepoints[1][idx]])
    stat_test.inputs.paired_files = paired_ims
else:
    raise NotImplementedError(
        "{0} test has not been implemented yet.".inputs["analysis"])

# Estimate model
# > Will compute beta maps from the GLM
level2estimate = Node(spm.model.EstimateModel(
    estimation_method={'Classical': 1}), name="level2estimate")

# Connect notes
analysis2nd.connect([
    (stat_test, level2estimate, [('spm_mat_file', 'spm_mat_file')])

])

# Run workflow
res = analysis2nd.run('MultiProc', plugin_args={'n_procs': 8})
nodes_wf = list(res.nodes)
outputs["Outputs"] = []
for node in nodes_wf:
    outputs["Outputs"].extend(node.result.outputs)

# Save output logs
log_dir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(log_dir):
    os.mkdir(log_dir)

for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, "{0}_{1}.json".format(
        name, runtime["timestamp"]))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    pprint(outputs)
