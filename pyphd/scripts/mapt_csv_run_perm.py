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
                             MODELS_COVARIATES, RSFMRI_TEMPORAL_COVARIATES,
                             MAPT_RSFMRI_COV_TYPES)
from pyphd.fc_analyses.conn import (
    extract_group, extract_clinical_data_for_mapt_metrics)
from pyphd.dti.utils import extract_dti_group_analysis_data
from pyphd.fsl.utils import (
    create_group_design_matrix, randomise, find_contrast, text2vest, palm)
from statsmodels.stats.multitest import multipletests


# Script documentation
DOC = """

Script used to run permutation tests on MAPT csv data.
------------------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_csv_run_perm.py \
    -v Regions_mean_cortical_thickness_M36-M0_qced_mapt_rsfmri_subs.csv \
    -c /home/lp259104/PHD/NOTES/Conn_input_Montpellier_Bordeaux_Toulouse_MAPT.csv \
    -a CDR_x_gpeMapt4c \
    -g gpeMapt4c \
    -f ~/fsl_init.sh \
    -o /tmp/test \
    -m two_way_ancova \
    -T M0 \
    -G \
    -D \
    -E \
    -O \
    -P \
    -C two_way_ancova/m1 \
    -M /home/lp259104/PHD/TOOLS/palm-alpha119/palm \
    -X lh.ct_cambridge_basc_multiscale_asym_scale036_cluster_6.0.txt \
       lh.ct_cambridge_basc_multiscale_asym_scale036_cluster_5.0.txt \
       rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_6.0.txt \
       lh.ct_cambridge_basc_multiscale_asym_scale036_cluster_20.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_20.0.txt \
       rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_1.0.txt \
       rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_5.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_1.0.txt \
       lh.ct_cambridge_basc_multiscale_asym_scale036_cluster_13.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_6.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_23.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_5.0.txt \
       rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_13.0.txt \
       lh_plus_rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_13.0.txt \
       lh.ct_cambridge_basc_multiscale_asym_scale036_cluster_23.0.txt \
       rh.ct_cambridge_basc_multiscale_asym_scale036_cluster_20.0.txt \
    -V 1 \
    -N 20


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
        prog="python mapt_csv_run_perm.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-c", "--clinical-datafile", type=str, required=True,
        help="Datafile with clinical data.")
    required.add_argument(
        "-v", "--variable-datafile", type=str, required=True, nargs="+",
        help="Datafile with continuous variable to test. This file must "
             "contain a 'Sid' col with subjects identifier and may have "
             "multiple column variables to test. There is one "
             "file by timepoint.")
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
        "-A", "--na-variable-datafile-values", type=str,
        help="NA value to delete subjects in variable datafile. If value is "
             "found for one subject (one line) in one col, this subject is "
             "discarded.")
    parser.add_argument(
        "-K", "--center", type=str, help="Center.")
    parser.add_argument(
        "-D", "--do-not-add-tp-covs", action="store_true",
        help="Do not add timepoint covariates.")
    parser.add_argument(
        "-G", "--extract-group", action="store_true",
        help="Option to extract subgroup.")
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        help="Timepoint for each input scalar list.")
    parser.add_argument(
        "-C", "--covariates-models", type=str, nargs="+",
        help="Specify models with specific covariates seperated by a slash "
             "(e.g : two_way_ancova/m1, two_way_ancova/m5)."
             "Equivalent of covariates field in json analyses file.")
    parser.add_argument(
        "-S", "--add-sid-col", type=str, default="Sid",
        help="Add Sid column to intermediate files generated by the script."
             "Can be useful to control if the file are correctly generated.")
    parser.add_argument(
        "-E", "--demean", action="store_true",
        help="Demean covariates in randomise.")
    parser.add_argument(
        "-O", "--two-tail", action="store_true",
        help="Do two tail tests.")
    parser.add_argument(
        "-P", "--save-parametric", action="store_true",
        help="Save test parametric outputs.")
    parser.add_argument(
        "-N", "--nb-perms", type=int, default=5000,
        help="Number of permutations.")
    parser.add_argument(
        "-M", "--palm-path", type=str,
        help="Path to palm bin file.")
    parser.add_argument(
        "-F", "--manually-demean", action="store_true",
        help="Manually demean covariates when creating design file.")
    parser.add_argument(
        "-J", "--group-json", type=is_file,
        help="Json with group extraction info", default=ANALYSES_DETAILS_JSON)
    parser.add_argument(
        "-X", "--exclude-variables", type=str, nargs="+",
        help="Exclude some variables from variable datafile if needed.")
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
    "tool": "mapt_csv_run_perm.py",
    "fsl_version": fsl_version
}

outputs = {}
outputs["palm_tstat_val_files"] = {}
outputs["palm_fstat_val_files"] = {}
outputs["palm_tstat_pval_unc_files"] = {}
outputs["palm_fstat_pval_unc_files"] = {}
outputs["palm_p_tstat_fwe_files"] = {}
outputs["palm_p_fstat_fwe_files"] = {}
outputs["palm_p_tstat_uncparap_files"] = {}
outputs["palm_p_fstat_uncparap_files"] = {}
outputs["palm_time_elasped_file"] = {}
outputs["palm_config_file"] = {}


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

if inputs["covariates_models"] is not None:
    covariates_info = {}
    for elt in inputs["covariates_models"]:
        model_list, model_list_covs = elt.split("/")
        if model_list not in covariates_info.keys():
            covariates_info[model_list] = []
        covariates_info[model_list].append(model_list_covs)
else:
    covariates_info = group_extraction_info[analysis_name]["covariates"]

if inputs["do_not_add_tp_covs"]:
    add_cov_spe_timepoints = False
else:
    add_cov_spe_timepoints = group_extraction_info[
        analysis_name]["add_cov_spe_timepoints"]

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
subjects_info = {}
commands_stats = {}

# Create analysis dir
analysis_dir = os.path.join(inputs["outdir"], inputs["analysis_name"])
if not os.path.isdir(analysis_dir):
    os.mkdir(analysis_dir)

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

for idx_tp, tp in enumerate(inputs["timepoints"]):

    # Get all variables to test
    varfile_data = pd.read_csv(varfiles[idx_tp])
    var_to_test = []
    for col in varfile_data.columns:
        if col != "Sid":
            if inputs["exclude_variables"] is not None:
                if col not in inputs["exclude_variables"]:
                    var_to_test.append(col)
            else:
                var_to_test.append(col)

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
        tp_dir = os.path.join(analysis_dir, tp)
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

            # > Extract clinical data
            outfile = extract_clinical_data_for_mapt_metrics(
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

            # > Extract subgroup if necessary
            extract_group_del_subjects = None
            if inputs["extract_group"]:
                groups_info = {}
                for elt, elt_data in inputs["groups_info"].items():
                    groups_info[int(elt)] = elt_data

                # Extract group
                outfile, extract_group_del_subjects = extract_group(
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

            # >> Add lines that were deleted fron extract group
            if extract_group_del_subjects is not None:
                for idx, row in extract_group_del_subjects.iterrows():
                    lines_deleted.append(
                        ",".join([str(x) for x in list(row)]))

            # > Write nan data
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
            if model == "two_way_anova" or model == "two_way_ancova":
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

            # > Get categorical variables
            if model_spe_data["covariates"] is not None:
                categorical_covs = []
                for cov in model_spe_data["covariates"]:
                    if MAPT_RSFMRI_COV_TYPES[cov] == "Categorical":
                        categorical_covs.append(cov)
            else:
                categorical_covs = None

            # Create design file
            if model_spe_data["covariates"] is not None:
                if inputs["add_sid_col"] is not None:
                    design_covs = model_spe_data["covariates"] + [
                        inputs["add_sid_col"]]
                else:
                    design_covs = model_spe_data["covariates"]
            else:
                if inputs["add_sid_col"] is not None:
                    design_covs = [inputs["add_sid_col"]]
                else:
                    design_covs = None
            (design_file, _, nb_gpe_cols,
             nb_subjects) = create_group_design_matrix(
                datafile=outfile,
                gpe_col=gpe_col,
                outdir=model_spe_data["outdir"],
                covariates=design_covs,
                categorical_covariates_to_levels=categorical_covs,
                delete_sid_col=inputs["add_sid_col"],
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

            # For each variable extract values and write them as palm input
            # files
            # > Create palm inputs dir
            palm_input_dir = os.path.join(
                    model_spe_data["outdir"], "palm_inputs")
            if not os.path.isdir(palm_input_dir):
                os.mkdir(palm_input_dir)

            # > Create palm outputs dir
            palm_output_dir = os.path.join(
                model_spe_data["outdir"], "palm_outputs")
            if not os.path.isdir(palm_output_dir):
                os.mkdir(palm_output_dir)

            # > Extract values and run palm
            with open(outfile, "rt") as open_file:
                lines = open_file.readlines()
            var_idx_corr = {}
            header = lines[0].split(",")
            for idx, var in enumerate(header):
                if inputs["exclude_variables"] is not None:
                    if var in inputs["exclude_variables"]:
                        continue
                var_idx_corr[var] = idx

            outputs["palm_tstat_val_files"][tp] = {}
            outputs["palm_fstat_val_files"][tp] = {}
            outputs["palm_tstat_pval_unc_files"][tp] = {}
            outputs["palm_fstat_pval_unc_files"][tp] = {}
            outputs["palm_p_tstat_fwe_files"][tp] = {}
            outputs["palm_p_fstat_fwe_files"][tp] = {}
            outputs["palm_p_tstat_uncparap_files"][tp] = {}
            outputs["palm_p_fstat_uncparap_files"][tp] = {}
            outputs["palm_time_elasped_file"][tp] = {}
            outputs["palm_config_file"][tp] = {}

            merge_results = {}
            for var in var_to_test:

                # >> Extract and write values
                val_col_values = []

                for line in lines[1:]:
                    val_col_values.append(line.split(",")[var_idx_corr[var]])
                infile_palm = os.path.join(
                    palm_input_dir, "{0}_palm_input.csv".format(var))
                with open(infile_palm, "wt") as open_file:
                    for val in val_col_values:
                        open_file.write(val)
                        open_file.write("\n")

                # > Run palm
                palm_output_basename = os.path.join(
                    palm_output_dir, "palm_{0}".format(var))

                (palm_tstat_val_files, palm_fstat_val_files,
                 palm_tstat_pval_unc_files, palm_fstat_pval_unc_files,
                 palm_p_tstat_fwe_files, palm_p_fstat_fwe_files,
                 palm_p_tstat_uncparap_files, palm_p_fstat_uncparap_files,
                 palm_time_elasped_file, palm_config_file) = palm(
                    indata=infile_palm,
                    design_file=design_mat_file,
                    contrast_file=constrat_con_file,
                    f_contrast=f_contrast_file,
                    output_basename=palm_output_basename,
                    nb_permutations=inputs["nb_perms"],
                    twotail=inputs["two_tail"],
                    saveparametric=inputs["save_parametric"],
                    alternate_palm_bin=inputs["palm_path"])

                # > Save outputs
                merge_results[var] = {}
                local_variables = locals().copy()
                for output_var in [palm_tstat_val_files, palm_fstat_val_files,
                                   palm_tstat_pval_unc_files,
                                   palm_fstat_pval_unc_files,
                                   palm_p_tstat_fwe_files,
                                   palm_p_fstat_fwe_files,
                                   palm_p_tstat_uncparap_files,
                                   palm_p_fstat_uncparap_files,
                                   palm_time_elasped_file, palm_config_file]:

                    # >> Find output_varname
                    output_varname = None
                    for varname_data in local_variables.items():
                        if type(varname_data[1]) == type(output_var):
                            if varname_data[1] == output_var:
                                output_varname = varname_data[0]

                    output_var_list = []
                    if output_var is None:
                        outputs[output_varname][tp][var] = None
                    else:
                        for fid in output_var:
                            with open(fid, "rt") as open_file:
                                val_fid = open_file.readlines()[0].strip("\n")
                                output_var_list.append(val_fid)
                        outputs[output_varname][tp][var] = output_var
                    merge_results[var][output_varname] = output_var_list

            # Merge all variables results
            # > Check if we take pvalue for t or f stat
            choose_fstat_pval = True
            if f_contrast_file is None:
                choose_fstat_pval = False

            # > Get pvalues and correct them with BH
            pvalues = []
            pvalues_fdr_corr = []
            for var in var_to_test:
                if choose_fstat_pval:
                    pval = merge_results[var]["palm_fstat_pval_unc_files"]
                else:
                    pval = merge_results[var]["palm_tstat_pval_unc_files"]

                # > Check that there is no confusion between multiple contrasts
                # to choose results from
                if len(pval) > 1:
                    raise ValueError("Multiple contrast to choose pval from")
                pval = pval[0]
                pvalues.append(pval)
            pvalues = [float(x) for x in pvalues]
            _, pvalues_fdr_corr, _, _ = multipletests(
                pvalues, alpha=0.05, method='fdr_bh')

            stat_outfile = os.path.join(
                model_spe_data["outdir"], "statistics_summary_file.csv")
            with open(stat_outfile, "wt") as open_file:
                header = ["Variables", "Test", "T-value", "F-value", "P-value",
                          "Survive individual test (p < 0.05)",
                          "Corrected p-val (BH)",
                          "Survive FDR correction at q-FDR < 0.05"]
                open_file.write(",".join(header))
                open_file.write("\n")
                for idx_var, var in enumerate(var_to_test):
                    line = [var, "perm"]
                    line += ["/".join(merge_results[var][
                             "palm_tstat_val_files"])]
                    line += ["/".join(merge_results[var][
                             "palm_fstat_val_files"])]
                    line += [str(pvalues[idx_var])]
                    if pvalues[idx_var] < 0.05:
                        line += ["True"]
                    else:
                        line += ["False"]
                    line += [str(pvalues_fdr_corr[idx_var])]
                    if pvalues_fdr_corr[idx_var] < 0.05:
                        line += ["True"]
                    else:
                        line += ["False"]
                    open_file.write(",".join(line))
                    open_file.write("\n")
            outputs["statistics_summary_file"] = stat_outfile


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
outfile = os.path.join(inputs["outdir"], inputs["analysis_name"] + ".json")
with open(outfile, "wt") as open_file:
    outdata = {"output_files": outputs, "runtime": runtime, "inputs": inputs}
    json.dump(outdata, open_file,
              sort_keys=True, check_circular=True, indent=4)
if verbose > 0:
    print("[Outputs]:")
    print(outfile)
    pprint(outputs)
