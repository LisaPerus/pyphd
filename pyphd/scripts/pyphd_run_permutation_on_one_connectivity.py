#!/usr/bin/env python
# Script used to run FSL PALM on one connection after inputs have been prepared
# with pyphd_prepare_permutation_connectivities.py 
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

# Pyphd imports
from pyphd.fsl.utils import palm

# Script documentation
DOC = """

Script used to run FSL PALM on one connection after inputs have been prepared
with pyphd_prepare_permutation_on_one_connectivity
-----------------------------------------------------------------------------

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_run_permutation_on_one_connectivity.py \
    -i /nfs/neurospin/tmp/lperus2/MAPT_rsfmri/IM_vs_Placebo/M36/ttest/conn_connectivities_IM_1_Placebo_0_0399_rscores_all_subjects_at_POST_connections/ventral_DMN_2_to_post_Salience_4_connectivity_values.csv \
    -d /nfs/neurospin/tmp/lperus2/MAPT_rsfmri/IM_vs_Placebo/M36/ttest/design.mat \
    -c /nfs/neurospin/tmp/lperus2/MAPT_rsfmri/IM_vs_Placebo/M36/ttest/ttest.con \
    -n conn_connectivities_IM_1_Placebo_0_0399_rscores_all_subjects_at_POST \
    -f $HOME/fsl_init.sh \
    -o /nfs/neurospin/tmp/lperus2/MAPT_rsfmri/IM_vs_Placebo/M36/ttest/conn_connectivities_IM_1_Placebo_0_0399_rscores_all_subjects_at_POST_connections \
    -T
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
        prog="python pyphd_run_permutation_on_one_connectivity.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-conn-file", type=is_file, required=True,
        help="Csv file listing for connection values for each subject."
             "File for one connection.")
    required.add_argument(
        "-d", "--design-mat", type=is_file, required=True,
        help="Design mat file.")
    required.add_argument(
        "-c", "--contrast-con", type=is_file, required=True,
        help="Contrast con file.")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-n", "--outname",
        type=str, required=True, metavar="<path>",
        help="Outname for palm outdir.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-F", "--f-contrast", type=is_file,
        help="F contrast txt file. If set :  test the t contrasts in a single "
             "F contrast.")
    parser.add_argument(
        "-A", "--alternate-palm", type=is_file,
        help="Alternate palm bin file.")
    parser.add_argument(
        "-T", "--two-tail", action="store_true",
        help="If set, perform two-tailed test.")
    parser.add_argument(
        "-N", "--nb-permutations", type=int, default=10000,
        help="Number of permutation for permutation testing.")
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
    "tool": "pyphd_run_permutation_on_one_connectivity.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Prepare permutation tests on rsfmri connections...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

"""
Run Palm
"""
palm_outdir = os.path.join(
    inputs["outdir"], inputs["outname"] + "_palm_output")
if not os.path.isdir(palm_outdir):
    os.mkdir(palm_outdir)

conn_data_path, conn_data_ext = os.path.splitext(inputs["input_conn_file"])
conn_basename = os.path.basename(
    inputs["input_conn_file"]).replace(conn_data_ext, "")
palm_output_basename = os.path.join(
    palm_outdir, "{0}_palm".format(conn_basename))
if inputs["two_tail"]:
    palm_output_basename += "two_tail"
stat_val, pval_unc, p_fwe = palm(
    indata=inputs["input_conn_file"],
    design_file=inputs["design_mat"],
    contrast_file=inputs["contrast_con"],
    f_contrast=inputs["f_contrast"],
    output_basename=palm_output_basename,
    nb_permutations=inputs["nb_permutations"],
    twotail=inputs["two_tail"],
    alternate_palm_bin=inputs["alternate_palm"])
outputs["{0}_palm_stat_val".format(conn_basename)] = stat_val
outputs["{0}_palm_pval_unc".format(conn_basename)] = pval_unc
outputs["{0}_palm_p_fwe".format(conn_basename)] = p_fwe

"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(palm_outdir, "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)

output_basename = "palm" + conn_basename.replace(conn_data_ext, "")
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
