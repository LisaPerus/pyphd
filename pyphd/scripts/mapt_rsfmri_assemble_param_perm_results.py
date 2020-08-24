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
from decimal import Decimal

# Third party imports
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

# Script documentation
DOC = """
Assemble MAPT resting state analysis parametric and permutation test results
----------------------------------------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_assemble_param_perm_results.py \
    -a Ge2RF_DecMMSE0399MRI_MMSEinf30_APOE4_HRD_only_IM_1_Pl_0 \
    -p $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_TESTS/Ge2RF_DecMMSE0399MRI_MMSEinf30_APOE4_HRD_only_IM_1_Pl_0.json \
    -e $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PERMUTATION_TESTS/Ge2RF_DecMMSE0399MRI_MMSEinf30_APOE4_HRD_only_IM_1_Pl_0.json \
    -o /tmp/assemble

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_assemble_param_perm_results.py \
    -a CDR05_only_gpeMapt4c \
    -p $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_TESTS/CDR05_only_gpeMapt4c.json \
    -e $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PERMUTATION_TESTS/CDR05_only_gpeMapt4c.json \
    -o /tmp/assemble \
    -P M0:$MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PERMUTATION_TESTS/CDR05_only_gpeMapt4c/M0 \
       M36:$MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PERMUTATION_TESTS/CDR05_only_gpeMapt4c/M36 \
       M36-M0:$MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_ALL_SUBJECTS_AT_TRANSVERSAL/PERMUTATION_TESTS/CDR05_only_gpeMapt4c/M36-M0 \
    -S "M36-M0:$MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_ALL_SUBJECTS_AT_TRANSVERSAL/PERMUTATION_TESTS/statistics/M36-M0/<model>_conn_connectivities_CDR05_only_gpeMapt4c_0399_diff_rscores_POST-PRE*.csv"


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
        prog="python mapt_rsfmri_assemble_param_perm_results.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--analysis-name", type=str, required=True,
        help="Analysis name.")
    required.add_argument(
        "-p", "--parametric-json", type=is_file, required=True,
        help="Paramatric tests json summary file.")
    required.add_argument(
        "-e", "--permutation-json", type=is_file, required=True,
        help="Permutation tests json summary file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")
    parser.add_argument(
        "-T", "--timepoints", type=str, nargs="+",
        default=["M0", "M36", "M36-M0"], help="Timepoints")
    parser.add_argument(
        "-P", "--alt-perm-outdir", type=str, nargs="+",
        help="If previous permutation tests have been run for specific "
             "timepoints, specify alternative output "
             "directory with permutation results, written in the form "
             "timepoint:outdir. Eg : M0:/tmp")
    parser.add_argument(
        "-S", "--alt-perm-statfile-pattern", type=str, nargs="+",
        help="If previous permutation tests have been run for specific "
             "timepoints, specify alternative output "
             "file pattern containing assembled permutation results."
             "This file is the same as separate connections results stored "
             "in alt-perm-outdir, but it may happen that to free space disk "
             "only assembled perm file was kept and not perm outdir."
             "Argument is given in the form file:<model>_file, where <model> "
             "is a pattern that will be replaced by each model tested."
             "Eg : M0:<model>_conn_connectivities_CDR05_only_gpeMapt4c_0399_"
             "diff_rscores_POST-PRE*s.csv")
    parser.add_argument(
        "-V", "--verbose",
        type=int, choices=[0, 1, 2], default=1,
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
    "tool": "mapt_rsfmri_assemble_param_perm_results.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Assemble rsfmri param+perm data...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Functions
"""


def get_perm_conn_statvals(conn, type_of_test, perm_outdir=None,
                           alt_perm_dir=None, alt_perm_file=None,
                           expect_to_find_perm_data=True):
    """
    Return stat and pval values for a connection
    --------------------------------------------

    Steps:
    1) Determine stat filename patterns corresponding to the type of test run
    2) Check if permutation data exist for the connection in perm_outdir,
       alt_perm_outdir or alt_perm_file (in this order)
    3) Check if it conn permutation value was expected to be computed
    4) Return perm stat and p values

    Parameters:
    -----------
    conn: str
        connection name.
    type_of_test: str
        type of test. Can be ttest, glm, anova or ancova.
    perm_outdir: str
        path to permutation outdir. Can be None if it is not specified in json
        permutation file.
    alt_perm_dir: str
        if perm_outdir does not contain data for connection check if data
        is there.
    alt_perm_file: str
        if conn permutation result is not found in perm_outdir, alt_perm_dir
        check if it exist in alt_perm_file.
    expect_to_find_perm_data: bool
        if connection is not found in perm_outdir, alt_perm_dir or
        alt_perm_file it may not have been computed if its distribution was
        normal. If connection is not expected to be found and is not found in
        perm_outdir, alt_perm_dir or alt_perm_file returns "Not Computed"
        for stat and p vals.

    Returns:
    --------
    perm_stat_val: float
        perm stat val
    perm_pval : float
        perm pval
    """

    # Define stat and pval file pattern
    stat_file_pattern = None
    pval_file_pattern = None

    if type_of_test in ["ttest", "glm"]:
        stat_file_pattern = "{0}_*dat_tstat.csv".format(conn)
        pval_file_pattern = "{0}_*dat_tstat_uncp.csv".format(conn)
    elif type_of_test in ["anova", "ancova"]:
        stat_file_pattern = "{0}_*dat_fstat_c4.csv".format(conn)
        pval_file_pattern = "{0}_*dat_fstat_uncp_c4.csv".format(conn)
    else:
        raise ValueError("Unknown type of test {0}".format(type_of_test))

    # Check for conn permutation files in perm outdir
    files_in_perm_outdir = True
    if perm_outdir is not None:
        permutation_stat_file = glob.glob(
            os.path.join(perm_outdir, "*palm_output", stat_file_pattern))
        permutation_pval_file = glob.glob(
            os.path.join(perm_outdir, "*palm_output", pval_file_pattern))
        if len(permutation_stat_file) > 1:
            raise ValueError(
                "Multiple permutations file {0}".format(permutation_stat_file))
        if len(permutation_stat_file) == 0:
            files_in_perm_outdir = False
    else:
        files_in_perm_outdir = False

    # Check for conn permutation files in alt perm outdir
    files_in_alt_perm_outdir = True
    if not files_in_perm_outdir:
        if alt_perm_dir is None:
            files_in_alt_perm_outdir = False
        else:
            permutation_stat_file = glob.glob(
                os.path.join(
                    alt_perm_dir, type_of_test, "*palm_output",
                    stat_file_pattern))
            permutation_pval_file = glob.glob(
                os.path.join(
                    alt_perm_dir, type_of_test, "*palm_output",
                    pval_file_pattern))
            if len(permutation_stat_file) > 1:
                raise ValueError(
                    "Multiple permutations file {0}".format(
                        permutation_stat_file))
            if len(permutation_stat_file) == 0:
                files_in_alt_perm_outdir = False

    # Check for conn permutation files in alt perm statfile
    if (not files_in_perm_outdir) and (not files_in_alt_perm_outdir):
        if alt_perm_file is not None:

            # > Define column names for stat and pval
            stat_col = None
            pval_col = "P-value"
            if type_of_test in ["ttest", "glm"]:
                stat_col = "T-value"
            elif type_of_test in ["anova", "ancova"]:
                stat_col = "F-value"
            else:
                raise ValueError(
                    "Unknown type of test {0}".format(type_of_test))

            # > Read file
            stat_data = pd.read_csv(alt_perm_file, index_col=0)
            perm_stat_val = stat_data.loc[conn, stat_col]
            perm_pval = stat_data.loc[conn, pval_col]
        else:

            # > If not permutation data is not found it is either expected
            # because permutation may not have been run if connectivity
            # distribution was normal (often the case with the latest test,
            # to save time permutation were run only on connectiions with
            # non normal distribution)
            if expect_to_find_perm_data:
                raise ValueError(
                    "Could not find permutation data for conn {0}".format(
                        conn))
            else:
                perm_stat_val = "Not computed"
                perm_pval = "Not computed"
    else:
        # Read stat and pval and return values
        permutation_stat_file = permutation_stat_file[0]
        permutation_pval_file = permutation_pval_file[0]
        with open(permutation_stat_file, "rt") as open_file:
            perm_stat_val = open_file.readline().strip("\n")
        with open(permutation_pval_file, "rt") as open_file:
            perm_pval = open_file.readline().strip("\n")
        perm_stat_val = float(perm_stat_val)
        perm_pval = float(perm_pval)

    return perm_stat_val, perm_pval


def get_stats_values(conn, type_of_test, normal_distribution,
                     analysis_parametric_data,
                     perm_dir=None, alt_perm_dir=None, alt_perm_file=None):
    """
    Function to get parametric and permutation test values and pvalues
    ------------------------------------------------------------------

    Steps:
    1) Get parametric tests stat and pvals
    2) Get permutation tests stat and pvals
    3) Round stat values at 4 digits
    4) Check if parametric and permutation stat value are the same (not for all
       connections)
    5) Return stat and pvalues and if stat values have been checked between
       parametric and permutation tests.

    Parameters:
    -----------
    conn: str
        connection name.
    type_of_test: str
        type of test. Can be ttest, glm, anova or ancova.
    normal_distribution: bool
        if connection has normal distribution.
    analysis_parametric_data: pandas Dataframe
        dataframe with connections and results from R scripts.
    perm_dir: str
        path to permutation outdir. Can be None if it is not specified in json
        permutation file.
    alt_perm_dir: str
        if perm_outdir does not contain data for connection check if data
        is there.
    alt_perm_file: str
        if conn permutation result is not found in perm_outdir, alt_perm_dir
        check if it exist in alt_perm_file.

    Returns:
    --------
    parametric_stat: str/float
        parametric test stat value
    parametric_pval: str/float
        parametric test pvalue
    round_permutation_stat: str/float
        permutation test stat value, rounded.
    permutation_pval:
        permutation test p value
    check_adequacy_between_param_perm: bool
        whether stat values have been compared between parametric and
        permutation tests.
    """

    check_adequacy_between_param_perm = False

    # Get parametric stat and pval
    conn_name = conn.replace(".", "_to_")
    parametric_pval = analysis_parametric_data.loc[conn, "P-value"]
    if type_of_test == "ttest":
        parametric_stat = analysis_parametric_data.loc[conn, "T-value/Rank"]
    elif type_of_test == "glm":
        parametric_stat = analysis_parametric_data.loc[conn, "t-value"]
    else:
        parametric_stat = analysis_parametric_data.loc[conn, "F-value"]

    # Get permutation stat and pvals
    expect_to_find_perm_data = None
    if normal_distribution:
        expect_to_find_perm_data = False
    else:
        expect_to_find_perm_data = True
    perm_stat_val, permutation_pval = get_perm_conn_statvals(
        conn=conn_name,
        type_of_test=type_of_test,
        perm_outdir=perm_dir,
        alt_perm_dir=alt_perm_dir,
        alt_perm_file=alt_perm_file,
        expect_to_find_perm_data=expect_to_find_perm_data)

    # Round values and check that stats values are the same
    if perm_stat_val != "Not computed":
        perm_stat_val = Decimal(perm_stat_val)
        round_permutation_stat = round(perm_stat_val, 4)
        abs_permutation_stat = abs(round_permutation_stat)
    else:
        round_permutation_stat = "Not computed"

    parametric_stat = Decimal(parametric_stat)
    round_parametric_stat = round(parametric_stat, 4)
    abs_parametric_stat = abs(round_parametric_stat)

    # Check if parametric and permutation stat value are the same
    if perm_stat_val != "Not computed":
        if type_of_test in ["ttest", "anova"]:

            # > If non parametric was used in parametric tests scripts
            # (by default R scripts perform non parametric tests if connection
            # distribution is not normal) or if Welch correction was used for
            # ttest do not check if permutation stat value and stat value
            # from R 'parametric test' scripts are the same (they will not be).
            # There should be at least one connection with normal distribution
            # where permutation test was also run to check that parametric
            # and permutation scripts to the same test (same absolute stat
            # value)
            if analysis_parametric_data.loc[conn, "Test"] not in [
                "Wilcoxon rank sum test",
                "Wilcoxon rank sum test with continuity correction",
                "Welch Two Sample t-test",
                    "Kruskal-Wallis"]:
                if abs_parametric_stat != abs_permutation_stat:
                    print(conn_name, type_of_test)
                    raise ValueError(
                        "Different stat values {0} and {1}".format(
                            str(abs_parametric_stat),
                            str(abs_permutation_stat)))
                check_adequacy_between_param_perm = True
        else:
            if abs_parametric_stat != abs_permutation_stat:
                print(conn_name, type_of_test)
                raise ValueError(
                    "Different stat values {0} and {1}".format(
                        str(abs_parametric_stat), str(abs_permutation_stat)))
            check_adequacy_between_param_perm = True

    return (parametric_stat, parametric_pval, round_permutation_stat,
            permutation_pval, check_adequacy_between_param_perm)


"""
End Functions
"""
# Read data
with open(inputs["parametric_json"], "rt") as open_file:
    parametric_json = json.load(open_file)
with open(inputs["permutation_json"], "rt") as open_file:
    permutation_json = json.load(open_file)

# Check same number of subjects in parametric & permutation tests
if parametric_json["subjects_info"] != permutation_json["subjects_info"]:
    print(parametric_json["subjects_info"])
    print(permutation_json["subjects_info"])
    raise ValueError(
        "Different number of subjects in parametric and permutation json.")

# Prepare json output for parametric + permutations tests
json_output = {"subjects_info": permutation_json["subjects_info"],
               "commands_scripts": permutation_json["commands_scripts"],
               "output_files": {
               "with_covariates": {}, "without_covariates": {}}}

# Get alternative perm dir stat files
alternative_permdirs = {}
alternative_statfiles_patterns = {}
if inputs["alt_perm_outdir"] is not None:
    for elt in inputs["alt_perm_outdir"]:
        tp, permdir = elt.split(":")
        alternative_permdirs[tp] = permdir
if inputs["alt_perm_statfile_pattern"] is not None:
    for elt in inputs["alt_perm_statfile_pattern"]:
        tp, perm_statfile_pattern = elt.split(":")
        alternative_statfiles_patterns[tp] = perm_statfile_pattern

# Assemble data
for type_analysis in ["without_covariates", "with_covariates"]:
    for tp in inputs["timepoints"]:

        # Get for an analyses files with/without covariates
        analysis_parametric_file = parametric_json[
            "output_files"][type_analysis][tp].strip("\n")
        type_of_test = os.path.basename(analysis_parametric_file).split("_")[0]

        # Check if permutation dir exists in json file
        analysis_permutation_dir = None
        if tp in permutation_json["outdirs"][type_analysis].keys():
            analysis_permutation_dir = permutation_json["outdirs"][
                type_analysis][tp]

        print(inputs["analysis_name"], tp, type_of_test)

        # Read parametric data
        analysis_parametric_data = pd.read_csv(
            analysis_parametric_file, index_col=0)

        # Stores results
        connections = analysis_parametric_data.index
        perm_stat_values = []
        perm_pvalues = []

        # Check that a normally distributed connection has been confronted
        # to permutation test
        normally_distribution_connectivity_check = False
        nb_normally_distribution_connectivity_check = 0
        normal_distribution_col = "Normal distribution within the groups (Shapiro p-val > 0.05)"
        if normal_distribution_col not in analysis_parametric_data.columns:
            normal_distribution_col = "Normality of residuals (Shapiro p-val > 0.05)"

        # Check if alt_perm_dir is given for this timepoint
        alt_perm_dir = None
        if tp in alternative_permdirs.keys():
            alt_perm_dir = alternative_permdirs[tp]

        # Check if alt_perm_pattern is given for this timepoint and if
        # an alternative perm stat file can be found
        alt_perm_file = None
        if tp in alternative_statfiles_patterns.keys():
            alt_perm_pattern = alternative_statfiles_patterns[tp]
            alt_perm_pattern = alt_perm_pattern.replace(
                "<model>", type_of_test)
            alt_perm_file = glob.glob(alt_perm_pattern)
            if len(alt_perm_file) == 0:
                raise ValueError(
                    "Cannot find alt stat file from pattern {0}".format(
                        alternative_statfiles_patterns[tp]))
            if len(alt_perm_file) > 1:
                raise ValueError(
                    "Multiple alt files found from pattern {0} : {1}".format(
                        alternative_statfiles_patterns[tp], alt_perm_file))
            alt_perm_file = alt_perm_file[0]

        for conn in connections:
            conn_distribution = analysis_parametric_data.loc[
                conn, normal_distribution_col]
            (parametric_stat, parametric_pval, round_permutation_stat,
             permutation_pval,
             check_adequacy_between_param_perm) = get_stats_values(
                conn=conn,
                type_of_test=type_of_test,
                normal_distribution=conn_distribution,
                analysis_parametric_data=analysis_parametric_data,
                perm_dir=analysis_permutation_dir,
                alt_perm_dir=alt_perm_dir,
                alt_perm_file=alt_perm_file)

            # > Add perm stat value
            perm_stat_values.append(round_permutation_stat)

            # > If normal distribution, add to permutation pvalues parametric
            # pvalue
            if conn_distribution:
                perm_pvalues.append(parametric_pval)
            else:
                perm_pvalues.append(permutation_pval)

            # > If stat val difference between parametric and permutation was
            # tested, report it
            if check_adequacy_between_param_perm:
                normally_distribution_connectivity_check = True
                nb_normally_distribution_connectivity_check += 1

        # Check that stat val difference between parametric and permutation was
        # tested, for at least one connection
        # If it is not the case, fail the assembling of all results.
        if not normally_distribution_connectivity_check:
            raise ValueError("Difference between parametric and permutation "
                             "tests should be checked for at least "
                             "one connection.")

        # Create new output data
        new_analysis_data = analysis_parametric_data.copy()
        new_analysis_data["Permutation stat val"] = perm_stat_values
        new_analysis_data["Permutation pval (only for test with non normal data)"] = perm_pvalues

        # Correct using fdr
        rejected, fdr_corr_pvalues = fdrcorrection(
            pvals=new_analysis_data["Permutation pval (only for test with non normal data)"],
            alpha=0.05,
            method='indep')
        new_analysis_data["Corrected p-val (BH) (include permutation results)"] = fdr_corr_pvalues
        new_analysis_data["Survive FDR correction at q-FDR < 0.05 (include permutation results)"] = rejected

        # Create csv output for analysis and timepoint and model
        tp_dir = os.path.join(inputs["outdir"], tp)
        if not os.path.isdir(tp_dir):
            os.mkdir(tp_dir)
        out_csv = os.path.join(
            tp_dir, os.path.basename(analysis_parametric_file).replace(
                ".csv", "_plus_perm_results.csv"))
        new_analysis_data.to_csv(out_csv, index=True)

        # Add output file to json output log
        json_output["output_files"][type_analysis][tp] = out_csv


"""
Write logs
"""

# Write json output log
json_log_file = os.path.join(
    inputs["outdir"], inputs["analysis_name"] + ".json")
with open(json_log_file, "wt") as open_file:
    json.dump(json_output, open_file, sort_keys=True, check_circular=True,
              indent=4)
outputs["Json file"] = json_log_file

# Write inputs/outputs logs
logdir = os.path.join(inputs["outdir"], "logs", inputs["analysis_name"])
if not os.path.isdir(logdir):
    os.makedirs(logdir)
output_basename = "mapt_rsfmri_assemble_param_perm_results"
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
