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
from pyphd.constants import (COVARIATES_RSFMRI_ANALYSES, ANALYSES_DETAILS_JSON)
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

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_dti_run_tests.py \
    -i ${MEDIA_SCRIPT}/MAPT_DTI_LISA/FA_list_0399_common_subjects_prepost_M0.csv ${MEDIA_SCRIPT}/MAPT_DTI_LISA/FA_list_0399_common_subjects_prepost_M36.csv \
    -c ${MEDIA_SCRIPT}/MAPT_DTI_LISA/0399_common_subjects_prepost_clinical_data/subjects_info_0399_common_subjects_prepost.csv \
    -a IM_vs_Placebo \
    -f ~/fsl_init.sh \
    -o ${MEDIA_SCRIPT}/MAPT_DTI_LISA/STATISTICS
    -o /tmp/test \
    -T M0 M36 \
    -S Sid \
    -C Sid \
    -M ${MEDIA_SCRIPT}/MAPT_DTI/MAPT_RERUN/tbss/common_mask_M0_M36_M60_mean_FA_skel.nii.gz \
    -D FA_skel \
    -V 1


#############################################################################
# Quick command test
python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_dti_run_tests.py \
    -i ${MEDIA_SCRIPT}/MAPT_DTI_LISA_TEST/test_subjects_sample.txt \
    -c ${MEDIA_SCRIPT}/MAPT_DTI_LISA/0399_common_subjects_prepost_clinical_data/subjects_info_0399_common_subjects_prepost.csv \
    -a CDR05_vs_CDR0 \
    -f ~/fsl_init.sh \
    -o ${MEDIA_SCRIPT}/MAPT_DTI_LISA_TEST/test_without_masking_stack \
    -T M36 \
    -S Sid \
    -C Sid \
    -M ${MEDIA_SCRIPT}/MAPT_DTI/MAPT_RERUN/tbss/common_mask_M0_M36_M60_mean_FA_skel.nii.gz \
    -D FA_skel \
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
        "-i", "--input-scalar-list", type=is_file, required=True, nargs="+",
        help="CSV file listing for MNI warped scalar image. "
             "Header with Sid File, comma-separated."
             "Each csv file corresponds to a timepoint.")
    required.add_argument(
        "-c", "--clinical-datafile", type=is_file, required=True,
        help="File with clinical data information. "
             "Must be semicolon-separated.")
    required.add_argument(
        "-a", "--analysis-name", type=str, required=True,
        help="Analysis name.")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        help="Timepoint for each input scalar list.")
    parser.add_argument(
        "-G", "--groups-info", type=str,
        help="Alternative group info to the one in ANALYSES_DETAILS_JSON "
             " for subgroup extraction. Dict is directly passed, surrounded "
             "by simple quotes (double quotes for dict keys and values).")
    parser.add_argument(
        "-S", "--scalar-datafile-subcol", type=str, default="Sid",
        help="Colname of subjects IDs in scalar list datafile.")
    parser.add_argument(
        "-D", "--scalar-datafile-dticol", type=str, default="FA_skel",
        help="Colname of scalar files in scalar list datafile.")
    parser.add_argument(
        "-A", "--scalar", type=str, default="FA",
        help="Scalar name.")
    parser.add_argument(
        "-C", "--clinical-datafile-subcol", type=str, default="Sid",
        help="Colname of subjects IDs in clinical datafile.")
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
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "mapt_dti_run_tests.py"
}
outputs = {}

"""
Charge group extraction info
"""
with open(inputs["group_json"], "rt") as open_file:
    group_extraction_info = json.load(open_file)
analysis_name = inputs["analysis_name"]
if group_extraction_info[analysis_name]["extract_group"]:
    inputs["rename_cols"] = group_extraction_info[analysis_name]["rename_cols"]
    inputs["rename_file"] = group_extraction_info[analysis_name]["rename_file"]
    inputs["erase_cols"] = group_extraction_info[analysis_name]["erase_cols"]
    inputs["extract_group"] = group_extraction_info[
        analysis_name]["extract_group"]
else:
    inputs["rename_cols"] = {}
    inputs["rename_file"] = {}
    inputs["erase_cols"] = []

# If groups info is not specified take group_json groups_info data
if inputs["groups_info"] is None:
    if group_extraction_info[analysis_name]["extract_group"]:
        inputs["groups_info"] = group_extraction_info[
            analysis_name]["groups_info"]
        else:
            inputs["groups_info"] = {}

group_name = group_extraction_info[analysis_name]["group_name"]
conn_file_additional_covariates = None
if "additional_conn_file_covariates" in group_extraction_info[
        analysis_name].keys():
    conn_file_additional_covariates = group_extraction_info[
        analysis_name]["additional_conn_file_covariates"]

# Get timepoints and timepoints data
if inputs["timepoints"] is None:
    print(
        "Only processing scalar file {0}".format(
            inputs["input_scalar_list"][0]))
    timepoints_data = {"NA": inputs["input_scalar_list"][0]}
else:
    timepoints_data = OrderedDict()
    for idx, tp in enumerate(inputs["timepoints"]):
        timepoints_data[tp] = inputs["input_scalar_list"][idx]

"""
Process without covariates
"""
subjects_info = {"without_covariates": {},
                 "with_covariates": {}}
commands_stats = {"without_covariates": {},
                  "with_covariates": {}}
outputs["without_covariates"] = {}
outputs["with_covariates"] = {}
for tp, dti_file in timepoints_data.items():

    for type_analysis in ["without_covariates", "with_covariates"]:

        test = None
        if type_analysis == "without_covariates":
            test = group_extraction_info[analysis_name]["models"][0]
        elif type_analysis == "with_covariates":
            test = group_extraction_info[analysis_name]["models"][1]
        else:
            raise ValueError("Unknown analysis {0}".format(type_analysis))

        print(analysis_name, tp, test)

        outputs[type_analysis][tp] = {}
        outdir = os.path.join(
            inputs["outdir"], analysis_name, inputs["scalar"], tp, test)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        if type_analysis == "without_covariates":
            covariates = []
            cov_plus_add_cols = []
        else:

            # If APOE has been added as a group delete it will be added twice
            # with the covariates. Delete the first APOE from the covariates.
            if ("keep_apoe_cov_stats"
                    in group_extraction_info[analysis_name].keys()):
                if not group_extraction_info[analysis_name][
                        "keep_apoe_cov_stats"]:
                    covariates = [
                        x for x in COVARIATES_RSFMRI_ANALYSES if x != "APOE4"]
            else:
                covariates = COVARIATES_RSFMRI_ANALYSES

            if tp == "M0":
                covariates = covariates + ["Delay_IntIRMM0_V1_days"]
            else:
                covariates = covariates + ["intIRM36-V1"]

        if conn_file_additional_covariates is not None:
            cov_plus_add_cols = covariates + conn_file_additional_covariates
        else:
            cov_plus_add_cols = covariates

        if len(cov_plus_add_cols) == 0:
            cov_plus_add_cols = None
        if len(covariates) == 0:
            covariates = None

        # Extract file with scalar path and clinical
        outfile = extract_dti_group_analysis_data(
            dti_datafile=dti_file,
            clinical_datafile=inputs["clinical_datafile"],
            dti_data_sid_col=inputs["scalar_datafile_subcol"],
            clinical_data_sid_col=inputs["clinical_datafile_subcol"],
            group_name=group_name,
            outdir=outdir,
            covariates=cov_plus_add_cols,
            additional_outfile_name=None)

        # Change group_name to be the first column and
        # subgroup col to be the last
        new_cols_order = [group_name]
        outdata = pd.read_csv(outfile)
        cols = outdata.columns
        for col in cols:
            if col in [group_name]:
                continue
            if conn_file_additional_covariates is not None:
                if col in conn_file_additional_covariates:
                    continue
            new_cols_order.append(col)
        if conn_file_additional_covariates is not None:
            for col in conn_file_additional_covariates:
                new_cols_order.append(col)
        outdata = outdata[new_cols_order]
        outdata.to_csv(outfile, index=False)

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
            if "NA" in list(row) or "NaN" in list(row):
                lines_with_na.append(index)
                lines_deleted.append(",".join([str(x) for x in list(row)]))
        for idx in lines_with_na:
            outdata = outdata.drop(idx)
        nan_values_df = outdata[outdata.isnull().values]
        for index, row in nan_values_df.iterrows():
            lines_deleted.append(",".join([str(x) for x in list(row)]))
        nan_subjects_csv = os.path.join(
            outdir, "input_deleted_lines_nan_subjects.csv")
        with open(nan_subjects_csv, "wt") as open_file:
            open_file.write(",".join(list(outdata.columns)))
            open_file.write("\n")
            for line in lines_deleted:
                open_file.write(line)
                open_file.write("\n")
        outputs[type_analysis][tp]["Deleted NaN subjects"] = nan_subjects_csv
        outdata = outdata.dropna()
        outdata.to_csv(outfile, index=False)

        # > Stack 4D image
        outdata = pd.read_csv(outfile)
        scalar_col = inputs["scalar_datafile_dticol"]
        scalar_files = list(outdata[scalar_col])
        outname = os.path.basename(os.path.splitext(outfile)[0])
        outname += "_" + scalar_col
        outnii = os.path.join(outdir, outname + ".nii.gz")
        cmd = ["fslmerge", "-t", outnii]
        for fid in scalar_files:
            cmd += [fid]
        fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
        fslprocess()

        # > Mask if needed
        # Note 03/07/2020 : since the mask is used in randomise we DO NOT
        # need to mask the 4D stack image. However, to be able to determine
        # quickly which mask was used just by looking at the stack image
        # we still do it.
        if inputs["mask"] is not None:
            cmd = ["fslmaths", outnii, "-mas", inputs["mask"], outnii]
            fslprocess = FSLWrapper(cmd, shfile=inputs["fsl_config"])
            fslprocess()
        outputs[type_analysis][tp]["Stack image"] = outnii

        # > Create design matrix file
        gpe_col = group_name
        if gpe_col in inputs["rename_cols"].keys():
            gpe_col = inputs["rename_cols"][gpe_col]

        # Note 03/07/2020 : since we are looking at group difference
        # and not at mean within a group itself there is no need to demean
        # the covariates (see http://mumford.fmripower.org/mean_centering/
        # and https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;b0aa0023.1412)
        # However as it does not change anything and can be needed for other
        # future design we leave the demeaning option at True.
        design_file, _, nb_gpe_cols = create_group_design_matrix(
            datafile=outfile,
            gpe_col=gpe_col,
            outdir=outdir,
            covariates=covariates,
            demean=True)

        # > Find contrast file
        if covariates is None:
            nb_covariates = None
        else:
            nb_covariates = len(covariates)
        contrast_file, fts_txt_file = find_contrast(
            nb_gpe_cols=nb_gpe_cols,
            nb_covariates=nb_covariates,
            model=test)

        # > Transform contrast and design files
        design_mat_file = design_file.replace(".txt", ".mat")
        text2vest(
            indata=design_file,
            outdata=design_mat_file,
            fsl_sh=inputs["fsl_config"])
        outputs[type_analysis][tp]["Design mat file"] = design_mat_file

        # If contrast is in text file, convert in .con file
        constrat_con_file = os.path.join(
            outdir,
            os.path.basename(contrast_file).replace(".txt", ".con"))
        text2vest(
            indata=contrast_file,
            outdata=constrat_con_file,
            fsl_sh=inputs["fsl_config"])
        outputs[type_analysis][tp]["Contrast con file"] = constrat_con_file

        # If f contrast text file is specified, convert it to fts file
        if fts_txt_file is not None:
            f_contrast_file = os.path.join(
                outdir,
                os.path.basename(fts_txt_file).replace(".txt", ".fts"))
            text2vest(
                indata=fts_txt_file,
                outdata=f_contrast_file,
                fsl_sh=inputs["fsl_config"])
        else:
            f_contrast_file = None
        outputs[type_analysis][tp]["F Contrast fts file"] = f_contrast_file

        # > Run randomise
        out_test_file = os.path.join(outdir, test)
        stat_files, tfce_files, randomise_cmd = randomise(
            input_file=outnii,
            output_file=out_test_file,
            design_file=design_mat_file,
            contrast_file=constrat_con_file,
            f_contrast=f_contrast_file,
            fsl_sh=inputs["fsl_config"],
            mask=inputs["mask"],
            tfce_2d_opt=True)
        outputs[type_analysis][tp]["Stat files"] = stat_files
        outputs[type_analysis][tp]["TFCE FWE corrected files"] = tfce_files
        outputs[type_analysis][tp]["Randomise cmd"] = randomise_cmd


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], analysis_name, "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)

output_basename = "mapt_dti_run_tests_{0}".format(inputs["scalar"])
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
