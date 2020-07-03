#!/usr/bin/env python
# Script used to prepare permutation tests on rsfmri connections using FSL PALM
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
import progressbar

# Third party imports
import pandas as pd
import numpy as np

# Pyphd imports
from pyphd.fsl.utils import text2vest

# Script documentation
DOC = """
Script used to prepare permutation tests on rsfmri connections using FSL PALM
-----------------------------------------------------------------------------

Outputs design and contrast file and connectivities input files.

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_prepare_permutation_connectivities.py \
    -i /home/lp259104/PHD/DATA/RESULTS/MAPT/CONNECTIVITY_ANALYSIS/SBC_01_Grecius_POST/CONN_FIRST_LEVEL_CONNECTIVITY_RESULTS/PARAMETRIC_TESTS/conn_connectivities_IM_1_Placebo_0_0399_rscores_all_subjects_at_POST.txt \
    -g IM_1_Placebo_0 \
    -s 0/0 1/1 \
    -c /tmp/contrast.txt \
    -f $HOME/fsl_init.sh \
    -o /tmp/IM_1_Placebo_0_0399_rscores_all_subjects_at_POST

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_prepare_permutation_connectivities.py \
    -i /home/lp259104/PHD/DATA/RESULTS/MAPT/CONNECTIVITY_ANALYSIS/SBC_01_Grecius_POST/CONN_FIRST_LEVEL_CONNECTIVITY_RESULTS/PARAMETRIC_TESTS/conn_connectivities_IM_1_Placebo_0_0399_rscores_all_subjects_at_POSTage_sexe_NIVSCOL_APOE4_intIRM36-V1.txt \
    -g IM_1_Placebo_0 \
    -s 0/0 1/1 \
    -c /tmp/contrast.txt \
    -f $HOME/fsl_init.sh \
    -C age sexe NIVSCOL APOE4 intIRM36_V1 \
    -D \
    -o /tmp/IM_1_Placebo_0_0399_rscores_all_subjects_at_POSTage_sexe_NIVSCOL_APOE4_intIRM36-V1
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
        prog="python pyphd_prepare_permutation_connectivities.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-file", type=is_file, required=True,
        help="Csv file listing for each subject its group, values for the "
             "different dependant variables, and eventual covariates. "
             "Must be comma-separated.")
    required.add_argument(
        "-g", "--gpe-col", type=str, required=True,
        help="Group column name.")
    required.add_argument(
        "-s", "--gpe-separate-cols", type=str, required=True, nargs="+",
        help="Indicates how to divide one col with groups "
             "of subject in 0 and 1. 3:IM/0 4:ctrl/1 indicates that subjects "
             "3:IM will be in the first column, and 4:ctrl in the second col. "
             "It is very important because in order to known which contrast "
             "is specified.")
    required.add_argument(
        "-c", "--contrast", type=is_file, help="Contrast file (text file).")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-C", "--covariates", type=str, nargs="+",
        help="Covariates names in input csv file header.")
    parser.add_argument(
        "-D", "--demean", action="store_true", help="Demean covariates.")
    parser.add_argument(
        "-F", "--f-contrast", type=is_file,
        help="F contrast txt file. If set :  test the t contrasts in a single "
             "F contrast.")
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
    "tool": "pyphd_prepare_permutation_connectivities.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Prepare permutation tests on rsfmri connections...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

"""
Load data
"""
input_data_path, input_data_ext = os.path.splitext(inputs["input_file"])
indata = pd.read_csv(inputs["input_file"], dtype=str)

# Get rid of NaN data and save NaN subjects
nan_values_df = indata[indata.isnull().values]
nan_subjects_csv = os.path.join(
    inputs["outdir"], os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "_nan_subjects" + input_data_ext))
nan_values_df.to_csv(nan_subjects_csv, index=False)
indata = indata.dropna()

"""
Create design matrix
"""
gpe_order = {}
for gpe_info in inputs["gpe_separate_cols"]:
    gpe_info = gpe_info.split("/")
    gpe = gpe_info[0]
    gpe_col = int(float(gpe_info[1]))
    gpe_order[gpe_col] = gpe

# Get dependant variables values
not_dependant_variables = [inputs["gpe_col"]]
if inputs["covariates"] is not None:
    not_dependant_variables += inputs["covariates"]
dependant_variables = [
    x for x in indata.columns if x not in not_dependant_variables]

# Create empty dataframe for design matrix
# > Create new cols for design matrix df
cols = []
for gpe_col in sorted(gpe_order.keys()):
    cols.append(gpe_order[gpe_col].replace(":", "_") + "_" + str(gpe_col))
if inputs["covariates"] is not None:
    cols += inputs["covariates"]

# > Create new index for design matrix df
index_design = [x for x in range(indata.shape[0])]

# > Create empty design matrix df
design_mat_df = pd.DataFrame(index=index_design, columns=cols)

# > Fill design matrix by columns
# >> with group columns
for gpe_col in sorted(gpe_order.keys()):
    gpe_col_name = gpe_order[gpe_col].replace(":", "_") + "_" + str(gpe_col)
    col_values = []
    for elt in indata[inputs["gpe_col"]]:
        if elt == gpe_order[gpe_col]:
            col_values.append("1")
        else:
            col_values.append("0")
    design_mat_df[gpe_col_name] = col_values

# >> with covariates
if inputs["covariates"] is not None:
    for cov in inputs["covariates"]:
        col_values = indata[cov]
        if inputs["demean"]:
            mean_val_col = np.mean([float(x) for x in col_values])
            col_values = [float(x) - mean_val_col for x in col_values]
        design_mat_df[cov] = col_values

# > Save intermediate desin matrix output
out_design_csv = os.path.join(inputs["outdir"], "design.csv")
design_mat_df.to_csv(out_design_csv, index=False)
outputs["Design file with colnames"] = out_design_csv

# > Save design mat file
design_txt_file = os.path.join(inputs["outdir"], "design.txt")
design_mat_df.to_csv(design_txt_file, index=False, header=False, sep=" ")


"""
Tranform design matrix/contrast files to valid PALM inputs
"""
design_mat_file = design_txt_file.replace(".txt", ".mat")
text2vest(
    indata=design_txt_file,
    outdata=design_mat_file,
    fsl_sh=inputs["fsl_config"])
outputs["Design mat file"] = design_mat_file

# If contrast is in text file, convert in .con file
if inputs["contrast"].endswith(".txt"):
    constrat_file = os.path.join(
        inputs["outdir"],
        os.path.basename(inputs["contrast"]).replace(".txt", ".con"))
    text2vest(
        indata=inputs["contrast"],
        outdata=constrat_file,
        fsl_sh=inputs["fsl_config"])
else:
    constrat_file = inputs["contrast"]
outputs["Contrast con file"] = constrat_file

# If f contrast text file is specified, convert it to fts file
if inputs["f_contrast"] is not None:
    f_contrast_file = os.path.join(
        inputs["outdir"],
        os.path.basename(inputs["f_contrast"]).replace(".txt", ".fts"))
    text2vest(
        indata=inputs["f_contrast"],
        outdata=f_contrast_file,
        fsl_sh=inputs["fsl_config"])
else:
    f_contrast_file = None
outputs["F Contrast fts file"] = f_contrast_file


"""
Write connections input files and perform palm
"""
connectivity_dir = os.path.join(
    inputs["outdir"], os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "") + "_connections")
if not os.path.isdir(connectivity_dir):
    os.mkdir(connectivity_dir)
outputs["connectivity_files"] = {}
with progressbar.ProgressBar(max_value=len(dependant_variables),
                             redirect_stdout=True) as bar:
    for idx_conn, connection in enumerate(dependant_variables):
        print("Processing {0}...".format(connection))
        connections_values = indata[connection]
        conn_name = connection.replace(":", "_to_")

        # Write connections file
        connection_file = os.path.join(
            connectivity_dir, "{0}_connectivity_values.csv".format(conn_name))
        with open(connection_file, "wt") as open_file:
            for val in connections_values:
                open_file.write(str(val))
                open_file.write("\n")
        outputs["connectivity_files"][
            "{0}_connectivity_values".format(conn_name)] = connection_file


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
output_basename = "pyphd_prepare_permutation_connectivities"
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
