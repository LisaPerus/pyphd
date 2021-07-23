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
from collections import OrderedDict

# pyphd import
import pandas as pd
import numpy as np
from pyphd.fc_analyses.conn import (extract_connectivities, extract_group,
                                    parse_conn_roi_to_roi_output_textfile)
from pyphd.constants import (SCRIPTS_STATS, MODELS_COVARIATES,
                             RSFMRI_TEMPORAL_COVARIATES, ANALYSES_DETAILS_JSON,
                             MAPT_RSFMRI_COV_TYPES, CONN_NETWORKS_ANALYSES)
from pyphd.display.report import generate_pdf_struct_file

# Pyconnectomist imports
from pyconnectomist.utils.pdftools import generate_pdf

# Script documentation
DOC = """
Compute for conn ROI-to-ROI results post-hoc tests
--------------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_conn_posthoc_tests.py \
    -a CDR_x_gpeMapt4c \
    -i /media/lp259104/10086bbb-e609-4274-a6b1-fcb6e9e0f9cb/tmp/MAPT_rsfmri/MONTPELLIER_ONLY_CONN/ROI_TO_ROI/R7/CDR_x_gpeMapt4c/M36/anova/3_TFCE/1000_perms/conn_results.txt \
    -n R7 \
    -t M36 \
    -o /tmp/test \
    -V 1

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_conn_posthoc_tests.py \
    -a CDR_x_gpeMapt4c \
    -i /media/lp259104/10086bbb-e609-4274-a6b1-fcb6e9e0f9cb/tmp/MAPT_rsfmri/MONTPELLIER_ONLY_CONN/ROI_TO_ROI/R7/CDR_x_gpeMapt4c/M36/ancova/m1/3_TFCE/1000_perms/conn_results.txt \
    -n R7 \
    -t M36 \
    -o /tmp/test \
    -C m1 \
    -V 1


python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_conn_posthoc_tests.py \
    -a CDR_x_gpeMapt4c \
    -i /media/lp259104/10086bbb-e609-4274-a6b1-fcb6e9e0f9cb/tmp/MAPT_rsfmri/MONTPELLIER_ONLY_CONN/ROI_TO_ROI/R36/CDR_x_gpeMapt4c/M36/ancova/m1/6_NBS/5000_perms/conn_results.txt \
    -n R36 \
    -t M36 \
    -o /tmp/test \
    -O 6_NBS \
    -C m1 \
    -s F \
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
        prog="python mapt_rsfmri_run_parametric_tests.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--analysis-name", type=str, required=True,
        help="Analysis name.")
    required.add_argument(
        "-i", "--input-conn-file", type=is_file, required=True,
        help="Input file with conn results created using the export table "
             "option in Conn GUI for ROI-to-ROI analysis.")
    required.add_argument(
        "-s", "--statistic", type=str, required=True,
        help="Type of statistic used. T or F.")
    required.add_argument(
        "-n", "--network-analysis", type=str, required=True,
        help="Name of Conn analysis to select in constant "
             "CONN_NETWORKS_ANALYSES.")
    required.add_argument(
        "-t", "--timepoint", type=str, required=True,
        help="Timepoint.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-J", "--group-json", type=is_file,
        help="Json with group extraction info", default=ANALYSES_DETAILS_JSON)
    parser.add_argument(
        "-R", "--rename-conn", type=is_file,
        help="File used to rename conn in extract_connectivities."
             "Must be semicolon separated")
    parser.add_argument(
        "-N", "--rename-conn-from-conn-of-interest-file", type=is_file,
        help="File used to rename conn saved from conn results."
             "Must be semicolon separated")
    parser.add_argument(
        "-K", "--tukeytest", type=str,
        default="two_way_anova", help="Type of posthoc test to apply.")
    parser.add_argument(
        "-C", "--covmodel", type=str,
        help="Specify model if covariate have to be added.")
    parser.add_argument(
        "-D", "--do-not-add-tp-covs", action="store_true",
        help="Do not add timepoint covariates.")
    parser.add_argument(
        "-I", "--conn-version", type=str, default="Conn19c",
        help="Conn version.")
    parser.add_argument(
        "-E", "--center", type=str, default="0399",
        help="Specify center.")
    parser.add_argument(
        "-O", "--conncorr", type=str, default="3_TFCE",
        help="Type of correction used for ROI-to-ROI analysis. "
             "Needed to know the structure of input conn file.")
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
    "tool": "mapt_rsfmri_conn_posthoc_tests.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Run post hoc test...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Charge group extraction info
"""
with open(inputs["group_json"], "rt") as open_file:
    group_extraction_info = json.load(open_file)
analysis_name = inputs["analysis_name"]
inputs["groups_info"] = group_extraction_info[analysis_name]["groups_info"]
inputs["rename_cols"] = group_extraction_info[analysis_name]["rename_cols"]
inputs["rename_file"] = group_extraction_info[analysis_name]["rename_file"]
inputs["erase_cols"] = group_extraction_info[analysis_name]["erase_cols"]
inputs["group_name"] = group_extraction_info[analysis_name]["group_name"]
inputs["extract_group"] = group_extraction_info[analysis_name]["extract_group"]
models = group_extraction_info[analysis_name]["models"]
covariates_info = group_extraction_info[analysis_name]["covariates"]
if not inputs["do_not_add_tp_covs"]:
    add_cov_spe_timepoints = group_extraction_info[
        analysis_name]["add_cov_spe_timepoints"]
else:
    add_cov_spe_timepoints = None
conn_file_additional_covariates = None
if "additional_conn_file_covariates" in group_extraction_info[
        analysis_name].keys():
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]
if "subgroup_two_way_anova" in group_extraction_info[
        analysis_name].keys():
    subgroup_two_way_anova = group_extraction_info[
        analysis_name]["subgroup_two_way_anova"]
else:
    subgroup_two_way_anova = None

# Get path to conn data
conn_inputs = CONN_NETWORKS_ANALYSES[
    inputs["network_analysis"]][inputs["timepoint"]]
conn_of_interest = []  # All connections of interest

# Reorganize conn output file
outfile = os.path.join(inputs["outdir"], "reorganized_conn_results.txt")
parse_conn_roi_to_roi_output_textfile(
        conn_textfile=inputs["input_conn_file"],
        method=inputs["conncorr"],
        statistic=inputs["statistic"],
        conn_version=inputs["conn_version"],
        outfile=outfile)
outputs["Reorganised input conn results file"] = outfile

# Get all connections that were statistically significant
# > Check if conn of interest have to be renamed
if inputs["rename_conn_from_conn_of_interest_file"] is not None:
    rename_conn_from_conn_of_interest_file = {}
    with open(inputs[
              "rename_conn_from_conn_of_interest_file"], "rt") as open_file:
        for line in open_file.readlines():
            line = line.strip("\n").split(";")
            rename_conn_from_conn_of_interest_file[line[0]] = line[1]
else:
    rename_conn_from_conn_of_interest_file = None

reorganized_conn_data = pd.read_csv(outfile)
for idx_row, row in reorganized_conn_data.iterrows():
    if rename_conn_from_conn_of_interest_file is None:
        conn_of_interest.append([row["ROI1"], row["ROI2"]])
    else:
        if row["ROI1"] in rename_conn_from_conn_of_interest_file.keys():
            new_roi1 = rename_conn_from_conn_of_interest_file[row["ROI1"]]
        else:
            new_roi1 = row["ROI1"]
        if row["ROI2"] in rename_conn_from_conn_of_interest_file.keys():
            new_roi2 = rename_conn_from_conn_of_interest_file[row["ROI2"]]
        else:
            new_roi1 = row["ROI2"]
        conn_of_interest.append([new_roi1, new_roi2])

# Extract group and subgroup for all connections
# > Check if covariates have to be added
covariates = None
list_covariates = None
if inputs["covmodel"] is not None:
    list_covariates = MODELS_COVARIATES[inputs["covmodel"]]

    # some covariates may be added depending on the timepoint for
    # the data (e.g : cov -> delay between mri and beginning of
    # the treatment)
    if add_cov_spe_timepoints is not None:
        list_covariates = list_covariates + [
            RSFMRI_TEMPORAL_COVARIATES[inputs["timepoint"]]]
if list_covariates is None:
    if conn_file_additional_covariates is not None:
        add_covariates = conn_file_additional_covariates
    else:
        add_covariates = None
else:
    if conn_file_additional_covariates is not None:
        add_covariates = list_covariates + conn_file_additional_covariates
    else:
        add_covariates = list_covariates

# Rename connectivities after extraction if needed
if inputs["rename_conn"] is not None:
    rename_conns_data = pd.read_csv(inputs["rename_conn"], index_col=False)
    rename_conns = {}
    with open(inputs["rename_conn"], "rt") as open_file:
        for line in open_file.readlines():
            line = line.strip("\n").split(";")
            if line[0] in rename_conns.keys():
                raise ValueError("Doublon in renaming connectivity")
            rename_conns[line[0]] = line[1]
else:
    rename_conns = None

conn_file = extract_connectivities(
    inputs["group_name"], tp=inputs["timepoint"], center_name=inputs["center"],
    covariates=add_covariates, network=None,
    conn_file_pattern=conn_inputs["conn_file_pattern"],
    datafile=conn_inputs["datafile"],
    outdir=inputs["outdir"],
    conn_datapath=conn_inputs["conn_datapath"],
    tp_name=conn_inputs["timepoint_name"],
    rename_conns=rename_conns)

group_colname = inputs["group_name"]
if inputs["extract_group"]:
    groups_info = {}
    for elt, elt_data in inputs["groups_info"].items():
        groups_info[int(elt)] = elt_data
    conn_file = extract_group(
        conn_file, groups_info=groups_info,
        rename_cols=inputs["rename_cols"],
        erase_cols=inputs["erase_cols"],
        rename_file=inputs["rename_file"])
    if group_colname in inputs["rename_cols"]:
        group_colname = inputs["rename_cols"][group_colname]

# Extract only connection of interest
conn_data = pd.read_csv(conn_file, dtype=str, header=0)
all_conns = conn_data.columns
keep_conn = [group_colname]
if subgroup_two_way_anova is not None:
    keep_conn += [subgroup_two_way_anova]
if list_covariates is not None:
    for col in conn_data.columns:
        if col in [x.replace("-", "_") for x in list_covariates]:
            keep_conn += [col]
for conn in all_conns:
    conn = conn.split(":")
    for conn_to_compare in conn_of_interest:
        if (conn[0] in conn_to_compare) and (conn[1] in conn_to_compare):
            keep_conn.append(":".join(conn))
conn_data = conn_data[keep_conn]
conn_file_with_selected_conns = os.path.join(
    inputs["outdir"], os.path.basename(conn_file).replace(
        "conn_connectivities", "reorganized_conn_connectivities"))
conn_data.to_csv(conn_file_with_selected_conns, index=0)
outputs["Input posthoc test"] = conn_file_with_selected_conns

# Do Tukey test
statistic_dir = os.path.join(inputs["outdir"], "statistics")
if not os.path.isdir(statistic_dir):
    os.mkdir(statistic_dir)
if inputs["tukeytest"] == "two_way_anova":
    cmd = ["Rscript", "--vanilla", SCRIPTS_STATS["posthoc_two_way_anova"],
           "-d", conn_file_with_selected_conns,
           "-g", group_colname,
           "-s", subgroup_two_way_anova,
           "-o", statistic_dir,
           "-v"]
    if list_covariates is not None:
        covariates_type = []
        for cov in list_covariates:
            covariates_type.append(MAPT_RSFMRI_COV_TYPES[cov])
        cmd += ["-c", " ".join([x.replace("-", "_") for x in list_covariates]),
                "-t", " ".join(covariates_type)]
    print(" ".join(cmd))
else:
    raise ValueError("Post hoc tests implement only for two-way anova.")

# Run cmd
proc = subprocess.Popen(cmd,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()
if proc.returncode == 1:
    raise ValueError(
        "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))
msg_out = stdout.decode(
    "utf-8").strip("\n").split("[1] \"Ouput : \"")
msg = msg_out[0]
output = msg_out[-1].strip("\n")
output = output.replace("[1] ", "").replace("\"", "")

# Assemble pngs created by R script
# > Get subgroups levels to filter png by subgroup level
subgroup_levels = list(np.unique(conn_data[subgroup_two_way_anova]))

# > Check if png files can be found for each connection
output_data = pd.read_csv(output)
output_outdir = os.path.dirname(output)
output_pngs = []
texts = []
nb_images_per_page = None
for idx, row in output_data.iterrows():
    var = row["Variables"]
    for level in subgroup_levels:
        var_pngs = glob.glob(
            os.path.join(
                output_outdir, "png_snaps", var + "_" + level + "_" + "*.png"))
        if len(var_pngs) > 0:
            output_pngs += sorted(var_pngs)
            texts += [var + "_" + level for x in range(len(var_pngs))]
            if nb_images_per_page is not None:
                if len(var_pngs) != nb_images_per_page:
                    raise ValueError(
                        "Different number of pngs for each variable.")
            nb_images_per_page = len(var_pngs)

# > Create pdf
json_struct = os.path.join(output_outdir, "struct.json")
generate_pdf_struct_file(
    images=output_pngs,
    out_file=json_struct,
    texts=texts,
    nb_im_per_pages=nb_images_per_page,
    style="TwoCol")
tic = datetime.now()
report_file = os.path.join(output_outdir, "summary_statistics.pdf")
generate_pdf(
    datapath=output_outdir,
    struct_file=json_struct,
    author="NA",
    client="NA",
    poweredby="NA",
    project="NA",
    timepoint="NA",
    subject="NA",
    date="{0}-{1}-{2}".format(tic.year, tic.month, tic.day),
    title="Statistics report",
    filename=report_file,
    pagesize=None,
    left_margin=10,
    right_margin=10,
    top_margin=20,
    bottom_margin=20,
    show_boundary=False)
outputs["Statistics report file"] = report_file


"""
Write outputs
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(
        logdir, "mapt_rsfmri_conn_posthoc_tests_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[Output]")
    pprint(outputs)
