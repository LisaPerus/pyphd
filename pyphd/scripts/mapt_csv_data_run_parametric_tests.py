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
from pyphd.fc_analyses.conn import extract_group, extract_clinical_data_for_mapt_metrics
from pyphd.constants import (SCRIPTS_STATS, MODELS_COVARIATES,
                             RSFMRI_TEMPORAL_COVARIATES, ANALYSES_DETAILS_JSON,
                             MAPT_RSFMRI_COV_TYPES)

# Script documentation
DOC = """
Run parametric tests on MAPT csv data.
--------------------------------------


Steps
1) Extract subgroups if needed
2) Run test
3) Save test results + nb of subjects in each groups


cp $HOME/PHD/NOTES/MAPT_fc_analyses/avancee_article_MAPT_v9/SUB_ANAT_ANALYSES_RESULTS/Cortical_thickness_ROI_cambridge_Freesurfer_6/Regions_mean_cortical_thickness_M36-M0.csv /tmp/Regions_mean_cortical_thickness_M36-M0.csv
python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_csv_data_run_parametric_tests.py \
    -a CDR0_stable_only_O3_vs_noO3 \
    -c $HOME/PHD/NOTES/Conn_input_Montpellier_Bordeaux_Toulouse_MAPT.csv \
    -v /tmp/Regions_mean_cortical_thickness_M36-M0.csv \
    -t "M36-M0" \
    -m ttest glm \
    -g gpeMapt4c \
    -o /tmp/test \
    -A "Data not found" \
    -G \
    -N 0399 \
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
        "-c", "--clinical-datafile", type=str, required=True,
        help="Datafile with clinical data.")
    required.add_argument(
        "-v", "--variable-datafile", type=str, required=True, nargs="+",
        help="Datafile with continuous variable to test. This file must "
             "contain a 'Sid' col with subjects identifier. There is one "
             "file by timepoint.")
    required.add_argument(
        "-t", "--timepoints", type=str, required=True, nargs="+",
        help="Timepoints.")
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
        "-A", "--na-variable-datafile-values", type=str,
        help="NA value to delete subjects in variable datafile. If value is "
             "found for one subject (one line) in one col, this subject is "
             "discarded.")
    parser.add_argument(
        "-G", "--extract-group", action="store_true",
        help="Option to extract subgroup.")
    parser.add_argument(
        "-D", "--do-not-add-tp-covs", action="store_true",
        help="Do not add timepoint covariates.")
    parser.add_argument(
        "-J", "--group-json", type=is_file,
        help="Json with group extraction info", default=ANALYSES_DETAILS_JSON)
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
    "tool": "mapt_csv_data_run_parametric_tests.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Run parametric test on MAPT csv data...")
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
covariates_info = group_extraction_info[analysis_name]["covariates"]

if not inputs["do_not_add_tp_covs"]:
    add_cov_spe_timepoints = group_extraction_info[
        analysis_name]["add_cov_spe_timepoints"]
else:
    add_cov_spe_timepoints = False

conn_file_additional_covariates = None
if "additional_conn_file_covariates" in group_extraction_info[
        analysis_name].keys():
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]
datafile = inputs["clinical_datafile"]
varfiles = inputs["variable_datafile"]


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

for idx_tp, tp in enumerate(inputs["timepoints"]):

    # For each model, extract group data and run test
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

            # > Extract clinical data
            # >> To double-check if file was correctly extracted delete_sid_col
            # can be put to False.
            # >> Extraction was checked for one example.
            group_colname = inputs["group_name"]
            covariate_append = model_spe_data["covariates_extract"] + [
                group_colname]
            replace_cov_names = {}
            for cov in covariate_append:
                if "-" in cov:
                    replace_cov_names[cov] = cov.replace("-", "_")

            # > Save one file with subjects ID to double check data extraction
            tp_varfile = extract_clinical_data_for_mapt_metrics(
                    datafile=datafile,
                    mri_metrics_csv_file=varfiles[idx_tp],
                    outdir=model_spe_data["outdir"],
                    subjects_na_var_delete=[inputs[
                        "na_variable_datafile_values"]],
                    columns_na_var_delete=["NaN", "NA"],
                    tp=tp,
                    center_name=inputs["center"],
                    covariates=covariate_append,
                    cov_add_file_beginning=group_colname,
                    replace_cov_names=replace_cov_names,
                    delete_sid_col=False)
            ext = os.path.splitext(tp_varfile)[1]
            os.rename(tp_varfile, tp_varfile.replace(ext, "_sids" + ext))

            # > Save without subjects IDs
            tp_varfile = extract_clinical_data_for_mapt_metrics(
                    datafile=datafile,
                    mri_metrics_csv_file=varfiles[idx_tp],
                    outdir=model_spe_data["outdir"],
                    subjects_na_var_delete=[inputs[
                        "na_variable_datafile_values"]],
                    columns_na_var_delete=["NaN", "NA"],
                    tp=tp,
                    center_name=inputs["center"],
                    covariates=covariate_append,
                    cov_add_file_beginning=group_colname,
                    replace_cov_names=replace_cov_names,
                    delete_sid_col=True)

            # > Extract by group
            if inputs["extract_group"]:
                groups_info = {}
                for elt, elt_data in inputs["groups_info"].items():
                    groups_info[int(elt)] = elt_data
                tp_varfile, _ = extract_group(
                    tp_varfile, groups_info=groups_info,
                    rename_cols=inputs["rename_cols"],
                    erase_cols=inputs["erase_cols"],
                    rename_file=inputs["rename_file"])
                if group_colname in inputs["rename_cols"]:
                    group_colname = inputs["rename_cols"][group_colname]

            # > Get number of subjects by group
            conn_data = pd.read_csv(tp_varfile).dropna()
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

            # > Create statistics_dir
            statistics_dir = os.path.join(
                model_spe_data["outdir"], "statistics")
            if not os.path.isdir(statistics_dir):
                os.mkdir(statistics_dir)

            # > If covariates get their type
            covariates_types = []
            if model_spe_data["covariates"] is not None:
                for cov in model_spe_data["covariates"]:
                    covariates_types.append(MAPT_RSFMRI_COV_TYPES[cov])

            # > Transform - to _ in covariate names for Rscript
            if model_spe_data["covariates"] is not None:
                covariates_model = []
                for cov in model_spe_data["covariates"]:
                    covariates_model.append(cov.replace("-", "_"))
            else:
                covariates_model = None

            # > If filename is too long R cannot read it, move to directory
            # and reduce namefile
            if len(tp_varfile) > 300:
                os.chdir(os.path.dirname(tp_varfile))
                tp_varfile_copy = os.path.join(
                    os.path.dirname(tp_varfile), "varfile.txt")
                shutil.copyfile(tp_varfile, tp_varfile_copy)
                tp_varfile = os.path.basename(tp_varfile_copy)
                statistics_dir = os.path.basename(statistics_dir)

            # > Run statistical test for ttest/anova
            if model in ["ttest", "anova"]:
                cmd = ["Rscript", "--vanilla", SCRIPTS_STATS[model],
                       "-d", tp_varfile,
                       "-o", statistics_dir,
                       "-g", group_colname,
                       "-v"]

            # > glm/ancova
            elif model in ["glm", "ancova"]:
                cmd = ["Rscript", "--vanilla", SCRIPTS_STATS[model],
                       "-d", tp_varfile,
                       "-o", statistics_dir,
                       "-g", group_colname,
                       "-c", " ".join(covariates_model),
                       "-t", " ".join(covariates_types)]

            elif model == "two_way_anova":

                # >> get subgroup col name
                subgroup_colname = group_extraction_info[
                    analysis_name]["subgroup_two_way_anova"]

                # >> two way anova
                cmd = ["Rscript", "--vanilla", SCRIPTS_STATS[model],
                       "-d", tp_varfile,
                       "-o", statistics_dir,
                       "-g", group_colname,
                       "-s", subgroup_colname]
                if covariates_model is not None:
                    cmd += ["-c", " ".join(covariates_model)]
                cmd += ["-t", " ".join(covariates_types)]
            else:
                raise ValueError("Unknown model : {0}".format(model))

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
            output = msg_out[-1]
            output_test = output.replace("[1] ", "").replace("\"", "")

            if model_spe_data["has_covariates"]:
                if model_spe not in outputs["with_covariates"].keys():
                    outputs["with_covariates"][model_spe] = {}
                if model_spe not in commands_stats["with_covariates"].keys():
                    commands_stats["with_covariates"][model_spe] = {}
                outputs["with_covariates"][model_spe][tp] = output_test
                commands_stats[
                    "with_covariates"][model_spe][tp] = " ".join(cmd)
            else:
                outputs["without_covariates"][tp] = output_test
                commands_stats[
                    "without_covariates"][tp] = " ".join(cmd)
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
