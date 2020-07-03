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
from pyphd.constants import (COVARIATES_RSFMRI_ANALYSES,
                             MAPT_RSFMRI_CONTRAST_FILES, ANALYSES_DETAILS_JSON)


# Script documentation
DOC = """
Compute MAPT fc analysis
---------------------------------

Steps
1) Extract all connectivities
2) Extract subgroups if needed
3) Save test results + nb of subjects in each groups

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_run_prepare_permutations.py \
    -a Placebo_only_CDR0_vs_CDR05 \
    -m ttest glm \
    -g scoreCDR1 \
    -o /tmp \
    -G \
    -V 2

#     -J ~/PHD/SCRIPTS/fMRI_ANALYSIS/trash_scripts/rerun_analysis_inter_intra_networks.json \
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
        prog="python mapt_rsfmri_run_prepare_permutations",
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
        "-F", "--fsl-sh", type=is_file, default=os.path.join(
            os.getenv("HOME"), "fsl_init.sh"), help="Path to FSL sh file")
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        default=["M0", "M36", "M36-M0"], help="Timepoints.")
    parser.add_argument(
        "-E", "--expected-nb-conn", type=int,
        help="Expected number of connections. If set, check that there is the"
              "right number of connections.")
    parser.add_argument(
        "-N", "--center", type=str, default="0399", help="Center.")
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
    "tool": "mapt_rsfmri_run_prepare_permutations"
}
outputs = {}
if verbose > 0:
    pprint("[info] Run fc analysis...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


# Get prepare permutation script
prepare_permutation_script = os.path.join(
    os.path.dirname(__file__), "pyphd_prepare_permutation_connectivities.py")

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
if ("additional_conn_file_covariates"
        in group_extraction_info[analysis_name].keys()):
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]

"""
Process without covariates
"""
subjects_info = {"without_covariates": {},
                 "with_covariates": {}}
commands_stats = {"without_covariates": {},
                  "with_covariates": {}}
outdirs = {}
outdirs["without_covariates"] = {}
outdirs["with_covariates"] = {}

# Create analysis outdir
analysis_dir = os.path.join(inputs["outdir"], analysis_name)
if not os.path.isdir(analysis_dir):
    os.mkdir(analysis_dir)

for tp in inputs["timepoints"]:

    # TODO: extract_connectivities and extract_group overwrite connectivity
    # files generated by mapt_rsfmri_run_parametric_tests.py
    # Add options to both function not to overwrite and/or change outdir value
    # > Extract connectivity file
    conn_file = extract_connectivities(
        inputs["group_name"], tp=tp, center_name=inputs["center"],
        covariates=conn_file_additional_covariates, network=None,
        conn_file_pattern=None, datafile=None, outdir=None, conn_datapath=None,
        tp_name=None)

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
    # Note: 05/02/2020
    # changed conn_data = pd.read_csv(conn_file).dropna()
    # to conn_data = pd.read_csv(conn_file, dtype=str).dropna()
    # because there was a bug with Omega3_low and 0/0.0
    conn_data = pd.read_csv(conn_file, dtype=str).dropna()
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

    # > Prepare permutation input
    # >> Prepare for each group a column
    gpe_separate_cols = []
    for idx, gpe in enumerate(sorted(sub_groups)):
        gpe_separate_cols.append("{0}/{1}".format(str(gpe), str(idx)))

    # >> Create model outdir
    model_dir = os.path.join(analysis_dir, tp, inputs["models"][0])
    if not os.path.isdir(model_dir):
        os.makedirs(model_dir)

    # >> Find contrast
    if inputs["models"][0] == "anova":
        contrast_file = MAPT_RSFMRI_CONTRAST_FILES["anova_4gpes"][0]
        contrast_f_file = MAPT_RSFMRI_CONTRAST_FILES["anova_4gpes"][1]
    elif inputs["models"][0] == "ttest":
        contrast_file = MAPT_RSFMRI_CONTRAST_FILES["ttest"]
        contrast_f_file = None
    else:
        raise ValueError("Unknown test : {0}".format(inputs["models"][0]))
    cmd = ["python3", prepare_permutation_script,
           "-i", conn_file,
           "-g", group_colname,
           "-s"]
    cmd += gpe_separate_cols
    cmd += ["-c", contrast_file, "-f", inputs["fsl_sh"]]
    cmd += ["-o", model_dir]
    if contrast_f_file is not None:
        cmd += ["-F", contrast_f_file]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))
    commands_stats["without_covariates"][tp] = " ".join(cmd)
    outdirs["without_covariates"][tp] = model_dir

"""
Process with covariates
"""
for tp in inputs["timepoints"]:

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
    conn_file = extract_connectivities(
        inputs["group_name"], tp=tp, center_name=inputs["center"],
        covariates=covariates_extract_conn)

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
    # Note: 05/02/2020
    # changed conn_data = pd.read_csv(conn_file).dropna()
    # to conn_data = pd.read_csv(conn_file, dtype=str).dropna()
    # because there was a bug with Omega3_low and 0/0.0
    conn_data = pd.read_csv(conn_file, dtype=str).dropna()
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
    covariates_tp_script = []
    for elt in covariates_tp:
        if elt != "intIRM36-V1":

            # Keep or not APOE E4 covariate
            if ("keep_apoe_cov_stats"
                    in group_extraction_info[analysis_name].keys()):
                if (not group_extraction_info[
                        analysis_name]["keep_apoe_cov_stats"] and
                        elt == "APOE4"):
                    continue
            covariates_tp_script.append(elt)
        else:
            covariates_tp_script.append("intIRM36_V1")

    # > If defined check that the number of expected connections is good
    if inputs["expected_nb_conn"] is not None:
        conn_data_cols = conn_data.columns
        nb_cols = len(conn_data_cols) - 1 - len(covariates_tp_script)
        if nb_cols != inputs["expected_nb_conn"]:
            raise ValueError(
                "Not expected number of connections in file {0}".format(
                    conn_file))

    # > Prepare permutation input
    # >> Prepare for each group a column
    gpe_separate_cols = []
    for idx, gpe in enumerate(sorted(sub_groups)):
        gpe_separate_cols.append("{0}/{1}".format(str(gpe), str(idx)))

    # >> Create model outdir
    model_dir = os.path.join(analysis_dir, tp, inputs["models"][1])
    if not os.path.isdir(model_dir):
        os.makedirs(model_dir)

    # >> Find contrast
    if inputs["models"][1] == "anova":
        contrast_file = MAPT_RSFMRI_CONTRAST_FILES["anova_4gpes"][0]
        contrast_f_file = MAPT_RSFMRI_CONTRAST_FILES["anova_4gpes"][1]
    elif inputs["models"][1] == "ancova":
        if len(covariates_tp_script) == 4:
            contrast_file = MAPT_RSFMRI_CONTRAST_FILES["ancova_4gpes_4covs"][0]
            contrast_f_file = MAPT_RSFMRI_CONTRAST_FILES[
                "ancova_4gpes_4covs"][1]
        elif len(covariates_tp_script) == 5:
            contrast_file = MAPT_RSFMRI_CONTRAST_FILES["ancova_4gpes_5covs"][0]
            contrast_f_file = MAPT_RSFMRI_CONTRAST_FILES[
                "ancova_4gpes_5covs"][1]
        else:
            msg = "No contrast file available for ancova with "
            msg += "{0} covariates".format(str(len(covariates_tp_script)))
            raise ValueError(msg)
    elif inputs["models"][1] == "glm":
        if len(covariates_tp_script) == 4:
            contrast_file = MAPT_RSFMRI_CONTRAST_FILES["glm_4covs"]
            contrast_f_file = None
        elif len(covariates_tp_script) == 5:
            contrast_file = MAPT_RSFMRI_CONTRAST_FILES["glm_5covs"]
            contrast_f_file = None
        else:
            raise ValueError(
                "No contrast file available for glm {0} covariates".format(
                    str(len(covariates_tp_script))))
    else:
        raise ValueError("Unknown model {0}".format(inputs["models"][1]))
    cmd = ["python3", prepare_permutation_script,
           "-i", conn_file,
           "-g", group_colname,
           "-s"]
    cmd += gpe_separate_cols
    cmd += ["-c", contrast_file, "-f", inputs["fsl_sh"]]
    cmd += ["-o", model_dir, "-D"]
    cmd += ["-C"] + covariates_tp_script
    if contrast_f_file is not None:
        cmd += ["-F", contrast_f_file]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))
    commands_stats["with_covariates"][tp] = " ".join(cmd)
    outdirs["with_covariates"][tp] = model_dir


"""
Write outputs
"""
outfile = os.path.join(inputs["outdir"], inputs["analysis_name"] + ".json")
with open(outfile, "wt") as open_file:
    outdata = {"outdirs": outdirs, "subjects_info": subjects_info,
               "commands_scripts": commands_stats}
    print(outdata)
    json.dump(outdata, open_file,
              sort_keys=True, check_circular=True, indent=4)
print(outfile)
