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
from pyphd.constants import (
    MAPT_RSFMRI_CONTRAST_FILES, ANALYSES_DETAILS_JSON,
    MAPT_RSFMRI_CONTRAST_FILES_INTERACTIONS, MODELS_COVARIATES,
    RSFMRI_TEMPORAL_COVARIATES, SCRIPT_DIR)


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
    SCRIPT_DIR, "GIT_REPOS", "pyphd", "pyphd", "scripts",
                "pyphd_prepare_permutation_connectivities.py")

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
if "covariates" in group_extraction_info[analysis_name].keys():
    covariates_info = group_extraction_info[analysis_name]["covariates"]
else:
    covariates_info = {}
add_cov_spe_timepoints = group_extraction_info[
    analysis_name]["add_cov_spe_timepoints"]

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

    # For each model, extract connectivities, group data and run test
    for model in inputs["models"]:

        # > Create tp/model dir and statistics dir
        # Note 06/08/2020 : previous code version could go wrong
        # if multiple processes with the same group_name and
        # subgroup extraction, with access to the SAME files were run in
        # parallel.
        # To avoid this issue analysis are run in specific analysis dir.
        tp_dir = os.path.join(analysis_dir, tp)
        if not os.path.isdir(tp_dir):
            os.mkdir(tp_dir)
        tp_model_dir = os.path.join(tp_dir, model)
        if not os.path.isdir(tp_model_dir):
            os.mkdir(tp_model_dir)

        # A model can be tested with different covariates, information about
        # all this is store in all_model_information
        # > outdir
        # > list of covs
        all_model_information = {}

        # Get if model needs covariates
        if model in covariates_info.keys():

            # Model can be tested with different sets of covariates
            for model_cov in covariates_info[model]:

                # create model cov dir
                model_cov_dir = os.path.join(tp_model_dir, model_cov)
                if not os.path.isdir(model_cov_dir):
                    os.mkdir(model_cov_dir)

                # get set of covariates
                list_covariates = MODELS_COVARIATES[model_cov]

                # some covariates may be added depending on the timepoint for
                # the data (e.g : cov -> delay between mri and beginning of
                # the treatment)
                if add_cov_spe_timepoints:
                    list_covariates = list_covariates + [
                        RSFMRI_TEMPORAL_COVARIATES[tp]]

                # append subgroup (named here additional covariates in the
                # covariates list) : used to create subgroups of subjects
                if conn_file_additional_covariates is not None:
                    covariates_extract_conn = list_covariates + \
                                              conn_file_additional_covariates
                else:
                    covariates_extract_conn = list_covariates

                all_model_information[model_cov] = {}
                all_model_information[model_cov][
                    "covariates"] = list_covariates
                all_model_information[model_cov][
                    "covariates_extract"] = covariates_extract_conn
                all_model_information[model_cov][
                    "has_covariates"] = True
                all_model_information[model_cov][
                    "outdir"] = model_cov_dir
        else:
            all_model_information[model] = {}
            all_model_information[model]["covariates"] = None
            all_model_information[model][
                "covariates_extract"] = conn_file_additional_covariates
            all_model_information[model]["has_covariates"] = False
            all_model_information[model]["outdir"] = tp_model_dir

        # For each model specification extract connectivities
        # and subgroup of subjects
        for model_spe, model_spe_data in all_model_information.items():
            conn_file = extract_connectivities(
                inputs["group_name"], tp=tp, center_name=inputs["center"],
                covariates=model_spe_data["covariates_extract"], network=None,
                conn_file_pattern=None, datafile=None,
                outdir=model_spe_data["outdir"], conn_datapath=None,
                tp_name=None)

            group_colname = inputs["group_name"]

            # If two-way anova create interaction column in extract_group
            # procedure
            if "subgroup_two_way_anova" in group_extraction_info[
                    analysis_name].keys():
                create_interaction_col = True
                interaction_cols = [group_extraction_info[analysis_name][
                                    "subgroup_two_way_anova"],
                                    inputs["group_name"]]
            else:
                create_interaction_col = False
                interaction_cols = []
            if inputs["extract_group"]:
                groups_info = {}
                for elt, elt_data in inputs["groups_info"].items():
                    groups_info[int(elt)] = elt_data
                conn_file = extract_group(
                    conn_file, groups_info=groups_info,
                    rename_cols=inputs["rename_cols"],
                    erase_cols=inputs["erase_cols"],
                    rename_file=inputs["rename_file"],
                    create_interaction_col=create_interaction_col,
                    interaction_cols=interaction_cols)
                if group_colname in inputs["rename_cols"]:
                    group_colname = inputs["rename_cols"][group_colname]

            # > Get number of subjects by group
            conn_data = pd.read_csv(conn_file).dropna()

            # > if interaction is tested change group colname
            if create_interaction_col:
                group_colname = "_:x:_".join(interaction_cols)

            # > get group subjects unique values
            sub_col = conn_data[group_colname]
            sub_groups = np.unique(sub_col)
            if not model_spe_data["has_covariates"]:
                subjects_info["without_covariates"][tp] = {}
                for gpe in sub_groups:
                    cpt_gpe = 0
                    for elt in sub_col:
                        if elt == gpe:
                            cpt_gpe += 1
                    subjects_info["without_covariates"][tp][
                        str(gpe)] = str(cpt_gpe)
            else:
                if model_spe not in subjects_info["with_covariates"].keys():
                    subjects_info["with_covariates"][model_spe] = {}
                subjects_info["with_covariates"][model_spe][tp] = {}
                for gpe in sub_groups:
                    cpt_gpe = 0
                    for elt in sub_col:
                        if elt == gpe:
                            cpt_gpe += 1
                    subjects_info["with_covariates"][model_spe][tp][
                        str(gpe)] = str(cpt_gpe)

            # > If defined check that number of expected connections is good
            if inputs["expected_nb_conn"] is not None:
                conn_data_cols = conn_data.columns
                if len(conn_data_cols) - 1 != inputs["expected_nb_conn"]:
                    raise ValueError(
                        "Not expected nb of connections in file {0}".format(
                            conn_file))

            # > Prepare permutation input
            # >> Prepare for each group a column
            gpe_separate_cols = []
            for idx, gpe in enumerate(sorted(sub_groups)):
                gpe_separate_cols.append("{0}/{1}".format(str(gpe), str(idx)))

            # Find contrast file
            # > in MAPT_RSFMRI_CONTRAST_FILES if it is a simple group
            # > comparison
            if model_spe_data["has_covariates"]:
                len_covs = str(len(model_spe_data["covariates"]))
            else:
                len_covs = "0"

            if not create_interaction_col:
                info_contrast = MAPT_RSFMRI_CONTRAST_FILES
                nb_groups = str(len(sub_groups))
            else:
                info_contrast = MAPT_RSFMRI_CONTRAST_FILES_INTERACTIONS

                # >> Find nb of subjects in each gpe
                interaction_cols_levels = {}
                for idx, var in enumerate(sorted(sub_groups)):
                    var = var.split("_:x:_")
                    for idx_level, gpe_level in enumerate(var):
                        if idx_level not in interaction_cols_levels.keys():
                            interaction_cols_levels[idx_level] = []
                        interaction_cols_levels[idx_level].append(gpe_level)
                for idx, var_data in interaction_cols_levels.items():
                    interaction_cols_levels[idx] = np.unique(var_data)
                length_levels_groups = []
                for col in sorted(interaction_cols_levels.keys()):
                    length_levels_groups.append(
                        len(interaction_cols_levels[col]))
                nb_groups = "_x_".join([str(x) for x in length_levels_groups])

            # Get t contrast file
            contrast_file = info_contrast[
                "nb_groups"][nb_groups]["nb_covs"][len_covs]["tcontrast"]

            # Get f contrast file if specified
            if "fcontrast" in info_contrast["nb_groups"][nb_groups][
                    "nb_covs"][len_covs].keys():
                contrast_f_file = info_contrast["nb_groups"][nb_groups][
                    "nb_covs"][len_covs]["fcontrast"]
            else:
                contrast_f_file = None

            # Create cmd
            cmd = ["python3.5", prepare_permutation_script,
                   "-i", conn_file,
                   "-g", group_colname,
                   "-s"]
            cmd += gpe_separate_cols
            cmd += ["-c", contrast_file, "-f", inputs["fsl_sh"]]
            cmd += ["-o", model_spe_data["outdir"]]

            if model_spe_data["has_covariates"]:

                # > Covariate names is transfromed from - to _
                # > in extract_group function, do the same for the script
                if model_spe_data["covariates"] is not None:
                    covariates_model = []
                    for cov in model_spe_data["covariates"]:
                        covariates_model.append(cov.replace("-", "_"))
                cmd += ["-C"] + covariates_model
                cmd += ["-D"]

            if contrast_f_file is not None:
                cmd += ["-F", contrast_f_file]
            proc = subprocess.Popen(cmd,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                raise ValueError(
                    "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))

            if model_spe_data["has_covariates"]:
                if model_spe not in outdirs["with_covariates"].keys():
                    outdirs["with_covariates"][model_spe] = {}
                if model_spe not in commands_stats["with_covariates"].keys():
                    commands_stats["with_covariates"][model_spe] = {}
                outdirs["with_covariates"][model_spe][
                    tp] = model_spe_data["outdir"]
                commands_stats["with_covariates"][model_spe][
                    tp] = " ".join(cmd)
            else:
                outdirs["without_covariates"][tp] = model_spe_data[
                    "outdir"]
                commands_stats[
                    "without_covariates"][tp] = " ".join(cmd)


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
