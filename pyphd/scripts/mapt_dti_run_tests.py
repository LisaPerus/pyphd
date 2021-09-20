#!/usr/bin/env python
# Script used to run tests on MAPT DTI data.
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
from pprint import pprint
from collections import OrderedDict
import pandas as pd
from pyconnectome.wrapper import FSLWrapper
from pyphd.constants import (COVARIATES_RSFMRI_ANALYSES, ANALYSES_DETAILS_JSON,
                             MODELS_COVARIATES, RSFMRI_TEMPORAL_COVARIATES)
from pyphd.fc_analyses.conn import extract_group
from pyphd.dti.utils import extract_dti_group_analysis_data
from pyphd.fsl.utils import (
    create_group_design_matrix, randomise, find_contrast, text2vest)


# Script documentation
DOC = """

Script used to run permutation tests on MAPT DTI data.
------------------------------------------------------

This script will ...
/!\ Dont forget common mask for longitudinal analyses /!\

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_dti_run_tests.py \
    -i /home/lp259104/PHD/NOTES/MAPT_DTI_Lisa/DTI_DATAFILE/dti_data_0399_M0_M36_M36-M0_only_subjects_common_M0_M36_clinical_data.csv \
    -a CDR_x_gpeMapt4c \
    -g gpeMapt4c \
    -f ~/fsl_init.sh \
    -o /tmp/test \
    -m two_way_anova \
    -T M0 M36 \
    -G \
    -A FA MD \
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
        prog="python mapt_dti_run_tests.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--dti-input-data", type=is_file, required=True,
        help="CSV file listing MNI warped scalar image for each subject at "
             "each timepoint and also covariates.")
    required.add_argument(
        "-a", "--analysis-name", type=str, required=True,
        help="Analysis name.")
    required.add_argument(
        "-g", "--group-name", type=str, required=True,
        help="Group name.")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-m", "--models", type=str, required=True, nargs="+",
        help="Models to utilize.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-G", "--extract-group", action="store_true",
        help="Option to extract subgroup.")
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        help="Timepoint for each input scalar list.")
    parser.add_argument(
        "-A", "--scalars", type=str, default=["FA"], nargs="+",
        help="Scalar names.")
    parser.add_argument(
        "-E", "--demean", action="store_true",
        help="Demean covariates in randomise.")
    parser.add_argument(
        "-N", "--nb-perms", type=int,
        help="Number of permutations.")
    parser.add_argument(
        "-F", "--manually-demean", action="store_true",
        help="Manually demean covariates when creating design file.")
    parser.add_argument(
        "-J", "--group-json", type=is_file,
        help="Json with group extraction info", default=ANALYSES_DETAILS_JSON)
    parser.add_argument(
        "-M", "--mask", type=is_file,
        help="File to mask skeleton with. Can be used if we want to compare "
             "data between different timepoints. E.g : do a stack at baseline "
             "that can be compared at another stack at three years."
             "The mask is the common mask between timepoints.")
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
fsl_version = FSLWrapper([], shfile=inputs["fsl_config"]).version
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "mapt_dti_run_tests.py",
    "fsl_version": fsl_version
}

outputs = {}

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
add_cov_spe_timepoints = group_extraction_info[
    analysis_name]["add_cov_spe_timepoints"]

conn_file_additional_covariates = None
if "additional_conn_file_covariates" in group_extraction_info[
        analysis_name].keys():
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]


"""
Process without covariates
"""
subjects_info = {}
commands_stats = {}
outputs = {}
outputs = {}

# Create analysis dir
analysis_dir = os.path.join(inputs["outdir"], inputs["analysis_name"])
if not os.path.isdir(analysis_dir):
    os.mkdir(analysis_dir)

for scalar in inputs["scalars"]:

    outputs[scalar] = {}
    commands_stats[scalar] = {}
    subjects_info[scalar] = {}
    outputs[scalar]["with_covariates"] = {}
    outputs[scalar]["without_covariates"] = {}
    commands_stats[scalar]["with_covariates"] = {}
    commands_stats[scalar]["without_covariates"] = {}
    subjects_info[scalar]["with_covariates"] = {}
    subjects_info[scalar]["without_covariates"] = {}

    # Create scalar dir
    scalar_dir = os.path.join(analysis_dir, scalar)
    if not os.path.isdir(scalar_dir):
        os.mkdir(scalar_dir)

    for tp in inputs["timepoints"]:

        # For each model, extract covariate data and characteristics
        # NB : code is the same as for the script
        # mapt_rsfmri_run_parametric_tests.py
        for model in inputs["models"]:

            # > Create tp/model dir and statistics dir
            # Note 06/08/2020 : previous code version could go wrong
            # if multiple processes with the same group_name and
            # subgroup extraction, with access to the SAME files were run in
            # parallel.
            # To avoid this issue analysis are run in specific analysis dir.
            tp_dir = os.path.join(scalar_dir, tp)
            if not os.path.isdir(tp_dir):
                os.mkdir(tp_dir)
            tp_model_dir = os.path.join(tp_dir, model)
            if not os.path.isdir(tp_model_dir):
                os.mkdir(tp_model_dir)

            # A model can be tested with different covariates,information about
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

                    # some covariates may be added depending on the timepoint
                    # for the data (e.g : cov -> delay between mri and
                    # beginning of the treatment)
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

            # For each model specification extract covariates of interest,
            # subgroup of subjects and scalar column
            for model_spe, model_spe_data in all_model_information.items():

                # Extract scalar col and covariates of interest
                # NB : order of columns is important if extract_group function
                # is used
                # group_name should always be first column and
                # additional_conn_file_covariates the last
                extract_cols = [inputs["group_name"]] + [
                    scalar + "_" + tp] + model_spe_data[
                        "covariates_extract"]
                outfile = extract_dti_group_analysis_data(
                    dti_datafile=inputs["dti_input_data"],
                    cols_to_keep=extract_cols, outdir=model_spe_data["outdir"],
                    additional_outfile_name=None)

                # > Extract subgroup if necessary
                if inputs["extract_group"]:
                    groups_info = {}
                    for elt, elt_data in inputs["groups_info"].items():
                        groups_info[int(elt)] = elt_data

                    # Extract group
                    outfile = extract_group(
                        outfile, groups_info=groups_info,
                        rename_cols=inputs["rename_cols"],
                        erase_cols=inputs["erase_cols"],
                        rename_file=inputs["rename_file"])

                # > If type of analysis with covariates, dropna() and NA
                lines_with_na = []
                lines_deleted = []
                outdata = pd.read_csv(outfile)
                for index, row in outdata.iterrows():
                    if ("NA" in list(row) or "NaN" in list(row) or
                            "nan" in list(row) or "NAN" in list(row)):
                        lines_with_na.append(index)
                        lines_deleted.append(
                            ",".join([str(x) for x in list(row)]))
                for idx in lines_with_na:
                    outdata = outdata.drop(idx)
                nan_values_df = outdata[outdata.isnull().values]
                for index, row in nan_values_df.iterrows():
                    lines_deleted.append(",".join([str(x) for x in list(row)]))
                nan_subjects_csv = os.path.join(
                    model_spe_data["outdir"],
                    "input_deleted_lines_nan_subjects.csv")
                with open(nan_subjects_csv, "wt") as open_file:
                    open_file.write(",".join(list(outdata.columns)))
                    open_file.write("\n")
                    for line in lines_deleted:
                        open_file.write(line)
                        open_file.write("\n")
                outdata = outdata.dropna()
                outdata.to_csv(outfile, index=False)

                # > Stack 4D image
                outdata = pd.read_csv(outfile)
                scalar_col = scalar + "_" + tp
                scalar_files = list(outdata[scalar_col])
                outname = os.path.basename(os.path.splitext(outfile)[0])
                outnii = os.path.join(
                    model_spe_data["outdir"], outname + ".nii.gz")
                cmd = ["fslmerge", "-t", outnii]
                for fid in scalar_files:
                    cmd += [fid]
                fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
                fslprocess()

                # > Mask if needed
                # Note 03/07/2020 : since the mask is used in randomise we DO
                # NOT need to mask the 4D stack image. However, to be able to
                # determine quickly which mask was used just by looking at the
                # stack image we still do it.
                if inputs["mask"] is not None:
                    cmd = ["fslmaths", outnii, "-mas", inputs["mask"], outnii]
                    fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
                    fslprocess()

                # > Create design matrix file
                gpe_col = inputs["group_name"]
                if gpe_col in inputs["rename_cols"].keys():
                    gpe_col = inputs["rename_cols"][gpe_col]

                # If test is a two way anova, add a merge option in
                # create_group_design_matrix to merge columns with main group
                # and second group
                # NOTE : Always order columns by growing number of elements
                # e.g : ["scoreCDR1", "gpeMapt4c"]
                # scoreCDR1 : CDR0 CDR0.5 -> 2 elts
                # gpeMapt4c : 1:omega3+IM 2:omega3 3:IM 4:ctrl -> 4 elts
                # Like this contrast will be found more easily
                if model == "two_way_anova":
                    merge_cols_non_ordered = [
                        gpe_col] + conn_file_additional_covariates
                    merge_cols = []
                    len_merge_cols = []
                    nb_elts_cols = {}
                    for col in merge_cols_non_ordered:
                        nb_col = len(set(outdata[col]))
                        if nb_col not in nb_elts_cols.keys():
                            nb_elts_cols[nb_col] = [col]
                        else:
                            nb_elts_cols[nb_col] += [col]
                    for nb_col in sorted(nb_elts_cols.keys()):
                        merge_cols += nb_elts_cols[nb_col]
                        len_merge_cols += [nb_col for i in range(
                            len(nb_elts_cols[nb_col]))]
                else:
                    merge_cols = None
                    len_merge_cols = None

                # Note 03/07/2020 : since we are looking at group difference
                # and not at mean within a group itself there is no need to
                # demean the covariates
                # (see http://mumford.fmripower.org/mean_centering/
                # and https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;b0aa0023.1412)
                # However as it does not change anything and can be needed for
                # other future design we can leave the demeaning option at
                # True.
                (design_file, _, nb_gpe_cols,
                 nb_subjects) = create_group_design_matrix(
                    datafile=outfile,
                    gpe_col=gpe_col,
                    outdir=model_spe_data["outdir"],
                    covariates=model_spe_data["covariates"],
                    merge_cols=merge_cols,
                    demean=inputs["manually_demean"])

                # > Find contrast file
                if model_spe_data["covariates"] is None:
                    nb_covariates = None
                else:
                    nb_covariates = len(model_spe_data["covariates"])

                # Find contrast file
                contrast_file, fts_txt_file = find_contrast(
                    nb_gpe_cols=nb_gpe_cols,
                    nb_covariates=nb_covariates,
                    model=model,
                    multiple_levels_nb=len_merge_cols)

                # > Transform contrast and design files
                design_mat_file = design_file.replace(".txt", ".mat")
                text2vest(
                    indata=design_file,
                    outdata=design_mat_file,
                    fsl_sh=inputs["fsl_config"])

                # If contrast is in text file, convert in .con file
                constrat_con_file = os.path.join(
                    model_spe_data["outdir"],
                    os.path.basename(contrast_file).replace(".txt", ".con"))
                text2vest(
                    indata=contrast_file,
                    outdata=constrat_con_file,
                    fsl_sh=inputs["fsl_config"])

                # If f contrast text file is specified, convert it to fts file
                if fts_txt_file is not None:
                    f_contrast_file = os.path.join(
                        model_spe_data["outdir"],
                        os.path.basename(fts_txt_file).replace(".txt", ".fts"))
                    text2vest(
                        indata=fts_txt_file,
                        outdata=f_contrast_file,
                        fsl_sh=inputs["fsl_config"])
                else:
                    f_contrast_file = None

                # > Run randomise
                stat_dir = os.path.join(model_spe_data["outdir"], "stats")
                if not os.path.isdir(stat_dir):
                    os.mkdir(stat_dir)
                out_test_file = os.path.join(stat_dir, model + "_" + model_spe)
                stat_files, tfce_files, randomise_cmd = randomise(
                    input_file=outnii,
                    output_file=out_test_file,
                    design_file=design_mat_file,
                    contrast_file=constrat_con_file,
                    f_contrast=f_contrast_file,
                    fsl_sh=inputs["fsl_config"],
                    mask=inputs["mask"],
                    tfce_2d_opt=True,
                    nb_perms=inputs["nb_perms"],
                    demean=inputs["demean"])

                # Save outputs
                if model_spe_data["has_covariates"]:
                    if model_spe not in outputs[
                            scalar]["with_covariates"].keys():
                        outputs[scalar]["with_covariates"][model_spe] = {}
                    if model_spe not in commands_stats[
                            scalar]["with_covariates"].keys():
                        commands_stats[
                            scalar]["with_covariates"][model_spe] = {}
                    if model_spe not in subjects_info[
                            scalar]["with_covariates"].keys():
                        subjects_info[
                            scalar]["with_covariates"][model_spe] = {}
                    outputs[scalar]["with_covariates"][model_spe][tp] = {
                        "stat_files": stat_files,
                        "tfce_fwe_corr_files": tfce_files}
                    commands_stats[scalar][
                        "with_covariates"][model_spe][tp] = randomise_cmd
                    subjects_info[scalar][
                        "with_covariates"][model_spe][tp] = nb_subjects
                else:

                    outputs[scalar]["without_covariates"][tp] = {
                        "stat_files": stat_files,
                        "tfce_fwe_corr_files": tfce_files}
                    commands_stats[scalar][
                        "without_covariates"][tp] = randomise_cmd
                    subjects_info[scalar][
                        "without_covariates"][tp] = nb_subjects


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
outfile = os.path.join(inputs["outdir"], inputs["analysis_name"] + ".json")
with open(outfile, "wt") as open_file:
    outdata = {"output_files": outputs, "subjects_info": subjects_info,
               "commands_randomise": commands_stats, "runtime": runtime,
               "inputs": inputs}
    json.dump(outdata, open_file,
              sort_keys=True, check_circular=True, indent=4)
if verbose > 0:
    print("[Outputs]:")
    print(outfile)
    pprint(outputs)
