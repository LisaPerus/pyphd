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

# pyphd import
import pandas as pd
import numpy as np
from pyphd.fc_analyses.conn import extract_connectivities, extract_group
from pyphd.constants import (SCRIPTS_STATS, COVARIATES_RSFMRI_ANALYSES,
                             ANALYSES_DETAILS_JSON, CONN_INPUTS)

# Script documentation
DOC = """
Compute MAPT fc analysis using parametric tests
-----------------------------------------------

NOTE: script previously named rerun_analysis_inter_intra_networks.py

Steps
1) Extract connectivities
2) Extract subgroups if needed
3) Run test
4) Save test results + nb of subjects in each groups

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_run_parametric_tests.py \
    -a CDR05_vs_CDR0 \
    -m ttest glm \
    -g scoreCDR1 \
    -o /tmp \
    -J ~/PHD/SCRIPTS/fMRI_ANALYSIS/trash_scripts/rerun_analysis_inter_intra_networks.json \
    -G \
    -V 2

###############################################
    
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_run_parametric_tests.py \
    -a Placebo_only_CDR0_vs_CDR05 \
    -m ttest glm \
    -g scoreCDR1 \
    -o /tmp \
    -G \
    -V 2

#    -J ~/PHD/SCRIPTS/fMRI_ANALYSIS/trash_scripts/rerun_analysis_inter_intra_networks.json \

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
        "-m", "--models", type=str, required=True, nargs="+",
        help="Models to utilize.")
    required.add_argument(
        "-g", "--group-name", type=str, required=True,
        help="Group name.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-G", "--extract-group", action="store_true",
        help="Option to extract subgroup.")
    parser.add_argument(
        "-J", "--group-json", type=is_file,
        help="Json with group extraction info", default=ANALYSES_DETAILS_JSON)
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        default=["M0", "M36", "M36-M0"], help="Timepoints.")
    parser.add_argument(
        "-E", "--expected-nb-conn", type=int,
        help="Expected number of connections. If set, check that there is the"
              "right number of connections.")
    parser.add_argument(
        "-N", "--center", type=str, help="Center.")
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
    "tool": "mapt_rsfmri_run_parametric_tests.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Run fc analysis...")
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
if group_extraction_info[analysis_name]["extract_group"]:
    inputs["groups_info"] = group_extraction_info[analysis_name]["groups_info"]
    inputs["rename_cols"] = group_extraction_info[analysis_name]["rename_cols"]
    inputs["rename_file"] = group_extraction_info[analysis_name]["rename_file"]
    inputs["erase_cols"] = group_extraction_info[analysis_name]["erase_cols"]
else:
    inputs["groups_info"] = {}
    inputs["rename_cols"] = {}
    inputs["rename_file"] = {}
    inputs["erase_cols"] = []

conn_file_additional_covariates = None
if "additional_conn_file_covariates" in group_extraction_info[
        analysis_name].keys():
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]

"""
Process without covariates
"""
subjects_info = {"without_covariates": {},
                 "with_covariates": {}}
commands_stats = {"without_covariates": {},
                  "with_covariates": {}}
outputs["without_covariates"] = {}
outputs["with_covariates"] = {}

# Create analysis dir
analysis_dir = os.path.join(inputs["outdir"], inputs["analysis_name"])
if not os.path.isdir(analysis_dir):
    os.mkdir(analysis_dir)

for tp in inputs["timepoints"]:

    # > Create tp/model dir and statistics dir
    # Note 06/08/2020 : previous code version could go wrong
    # if multiple processes with the same group_name and subgroup extraction,
    # with access to the SAME files were run in parallel.
    # To avoid this issue analysis are run in specific analysis dir.
    tp_dir = os.path.join(analysis_dir, tp)
    if not os.path.isdir(tp_dir):
        os.mkdir(tp_dir)
    tp_model_dir = os.path.join(tp_dir, inputs["models"][0])
    if not os.path.isdir(tp_model_dir):
        os.mkdir(tp_model_dir)
    statistics_dir = os.path.join(tp_model_dir, "statistics")
    if not os.path.isdir(statistics_dir):
        os.mkdir(statistics_dir)

    # > Extract connectivity file
    conn_file = extract_connectivities(
        inputs["group_name"], tp=tp, center_name=inputs["center"],
        covariates=conn_file_additional_covariates, network=None,
        conn_file_pattern=None, datafile=None, outdir=tp_model_dir,
        conn_datapath=None, tp_name=None)

    # > Extract subgroup if necessary
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

    # > Get number of subjects by group
    conn_data = pd.read_csv(conn_file).dropna()
    sub_col = conn_data[group_colname]
    sub_groups = np.unique(sub_col)
    subjects_info["without_covariates"][tp] = {}
    for gpe in sub_groups:
        cpt_gpe = 0
        for elt in sub_col:
            if elt == gpe:
                cpt_gpe += 1
        subjects_info["without_covariates"][tp][str(gpe)] = str(cpt_gpe)

    # > If defined check that the number of expected connections is good
    if inputs["expected_nb_conn"] is not None:
        conn_data_cols = conn_data.columns
        if len(conn_data_cols) - 1 != inputs["expected_nb_conn"]:
            raise ValueError(
                "Not expected number of connections in file {0}".format(
                    conn_file))

    # > Run statistical test
    cmd = ["Rscript", "--vanilla", SCRIPTS_STATS[inputs["models"][0]],
           "-d", conn_file,
           "-o", statistics_dir,
           "-g", group_colname,
           "-v"]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))
    msg_out = stdout.decode("utf-8").strip("\n").split("[1] \"Ouput : \"")
    msg = msg_out[0]
    output = msg_out[-1]
    output_test = output.replace("[1] ", "").replace("\"", "")
    outputs["without_covariates"][tp] = output_test
    commands_stats["without_covariates"][tp] = " ".join(cmd)
    print(output_test)

"""
Process with covariates
"""
for tp in inputs["timepoints"]:

    # > Create tp/model dir and statistics dir
    tp_dir = os.path.join(analysis_dir, tp)
    if not os.path.isdir(tp_dir):
        os.mkdir(tp_dir)
    tp_model_dir = os.path.join(tp_dir, inputs["models"][1])
    if not os.path.isdir(tp_model_dir):
        os.mkdir(tp_model_dir)
    statistics_dir = os.path.join(tp_model_dir, "statistics")
    if not os.path.isdir(statistics_dir):
        os.mkdir(statistics_dir)

    # If APOE has been added as a group delete it will be added twice with
    # the covariates. Delete the first APOE from the covariates.
    if "keep_apoe_cov_stats" in group_extraction_info[analysis_name].keys():
        if not group_extraction_info[analysis_name]["keep_apoe_cov_stats"]:
            COVARIATES_RSFMRI_ANALYSES = [
                x for x in COVARIATES_RSFMRI_ANALYSES if x != "APOE4"]

    if tp == "M0":
        covariates_tp = COVARIATES_RSFMRI_ANALYSES + ["Delay_IntIRMM0_V1_days"]
    else:
        covariates_tp = COVARIATES_RSFMRI_ANALYSES + ["intIRM36-V1"]
    if conn_file_additional_covariates is not None:
        covariates_extract_conn = covariates_tp + conn_file_additional_covariates
    else:
        covariates_extract_conn = covariates_tp

    # Extract connectivities
    conn_file = extract_connectivities(inputs["group_name"], tp=tp,
                                       center_name=inputs["center"],
                                       covariates=covariates_extract_conn,
                                       outdir=tp_model_dir)

    # > Extract subgroup if necessary
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

    # > Get number of subjects by group
    conn_data = pd.read_csv(conn_file).dropna()
    sub_col = conn_data[group_colname]
    sub_groups = np.unique(sub_col)
    subjects_info["with_covariates"][tp] = {}
    for gpe in sub_groups:
        cpt_gpe = 0
        for elt in sub_col:
            if elt == gpe:
                cpt_gpe += 1
        subjects_info["with_covariates"][tp][str(gpe)] = str(cpt_gpe)

    # > Run statistical test

    # >> Replace intIRM36-V1 by intIRM36_V1 for RScript
    covariates_tp_rscript = []
    for elt in covariates_tp:
        if elt != "intIRM36-V1":

            # Keep or not APOE E4 covariate
            if "keep_apoe_cov_stats" in group_extraction_info[
                    analysis_name].keys():
                if (not group_extraction_info[analysis_name][
                        "keep_apoe_cov_stats"] and elt == "APOE4"):
                    continue
            covariates_tp_rscript.append(elt)
        else:
            covariates_tp_rscript.append("intIRM36_V1")

    # > If defined check that the number of expected connections is good
    if inputs["expected_nb_conn"] is not None:
        conn_data_cols = conn_data.columns
        nb_conn = len(conn_data_cols) - 1 - len(covariates_tp_rscript)
        if nb_conn != inputs["expected_nb_conn"]:
            raise ValueError(
                "Not expected number of connections in file {0}".format(
                    conn_file))

    # >> Run script
    cmd = ["Rscript", "--vanilla", SCRIPTS_STATS[inputs["models"][1]],
           "-d", conn_file,
           "-o", statistics_dir,
           "-g", group_colname,
           "-c", " ".join(covariates_tp_rscript)]

    if inputs["models"][1] == "ancova":

        # Covariate type for "age", "sexe", "NIVSCOL", "APOE4" and delay
        # Check that APOE has to be kept in the covariates
        type_covariates = []
        for elt in covariates_tp_rscript:
            if elt in ["sexe", "NIVSCOL", "APOE4"]:
                type_covariates.append("Categorical")
            else:
                type_covariates.append("Continuous")
        cmd += ["-t", " ".join(type_covariates)]
    cmd += ["-v"]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))
    msg_out = stdout.decode("utf-8").strip("\n").split("[1] \"Ouput : \"")
    msg = msg_out[0]
    output = msg_out[-1]
    output_test = output.replace("[1] ", "").replace("\"", "")
    outputs["with_covariates"][tp] = output_test
    commands_stats["with_covariates"][tp] = " ".join(cmd)
    print(output_test)


"""
Write outputs
"""
outfile = os.path.join(inputs["outdir"], inputs["analysis_name"] + ".json")
with open(outfile, "wt") as open_file:
    outdata = {"output_files": outputs, "subjects_info": subjects_info,
               "commands_RScripts": commands_stats}
    json.dump(outdata, open_file,
              sort_keys=True, check_circular=True, indent=4)
if verbose > 0:
    print("[Outputs]:")
    print(outfile)
