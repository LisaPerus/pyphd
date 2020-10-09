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
from collections import OrderedDict

# Third party imports
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

# Script documentation
DOC = """
Assemble MAPT resting state analysis parametric and permutation test results
----------------------------------------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_assemble_param_perm_results.py \
    -r $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_TESTS/CDR_x_gpeMapt4c/M0/two_way_anova/m1/statistics/two_way_anova__conn_connectivities_CDR_x_gpeMapt4c_0399_rscores_common_subjects_timesteps_PREPOST_at_PREage_sexe_NIVSCOL_UNI_Delay_IntIRMM0_V1_days_.csv \
    -l $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PERMUTATION_TESTS/CDR_x_gpeMapt4c/M0/two_way_anova/m1/CDR_x_gpeMapt4c_palm_output/perm_results_assembled.csv \
    -k palm_dat_fstat_c10 palm_dat_fstat_c9 palm_dat_fstat_uncp_c10 palm_dat_fstat_uncp_c9 palm_dat_fstat_uncparap_c10 palm_dat_fstat_uncparap_c9 palmtwo_tail_dat_tstat_c8 \
       palmtwo_tail_dat_tstat_uncp_c8 palmtwo_tail_dat_tstat_uncparap_c8 \
    -o /tmp/test \
    -a mainGroup:secondGroup_P \
    -m palm_dat_fstat_uncp_c10 \
    -t mainGroup:secondGroup_F_splithere_palm_dat_fstat_c10 mainGroup_gpeMapt4c_F_splithere_palm_dat_fstat_c9 \
    -C mainGroup:secondGroup_P_splithere_palm_dat_fstat_uncparap_c10 mainGroup_gpeMapt4c_P_splithere_palm_dat_fstat_uncparap_c9 subGroup_scoreCDR1_P_splithere_palmtwo_tail_dat_tstat_uncparap_c8 \
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
        prog="python mapt_rsfmri_assemble_param_perm_results.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-r", "--parametric-file", type=is_file, required=True,
        help="Parametric file.")
    required.add_argument(
        "-l", "--permutation-file", type=is_file, required=True,
        help="Permutation file.")
    required.add_argument(
        "-k", "--perm-keep-cols", type=str, required=True, nargs="+",
        help="Columns to keep from permutation file.")
    required.add_argument(
        "-a", "--r-pval-col", required=True, type=str,
        help="R pval column of interest in parametric file. This parametric "
             "pvalues will be added to the final pvalues if connectivity "
             "distribution is normal.")
    required.add_argument(
        "-m", "--perm-pval-col", required=True, type=str,
        help="Pval column of interest in permutation file. This permutation "
             "pvalues will be added to the parametric pvalues if connectivity "
             "distribution is not normal.")
    required.add_argument(
        "-t", "--stat-check-col", required=True, type=str, nargs="+",
        help="Columns with stat values to check. First value is R stat column "
             "second column is Palm stat column. Each pair of column is "
             "separated by || . E.g : T-value/Rank||palm_dat_tstat_c1. "
             "Check will be done between R and palm stat values for each "
             "connection.")
    required.add_argument(
        "-o", "--outdir", type=is_directory, required=True, metavar="<path>",
        help="Output directory.")

    # Optional arguments
    parser.add_argument(
        "-C", "--check-parametric-pval-col", type=str, nargs="+",
        help="If parametric test pvalue has been computed through palm, "
             "compare this pvalue to R parametric pval. "
             "First value is R pval column, second column is Palm pval column."
             "Each pair of column is separated by || ."
             "E.g : P-value/Rank||palm_dat_tstat_uncparap_c1. "
             "Check will be done between R and palm pvalues for each "
             "connection.")
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


def compare_palm_and_r_values(data, filename, palm_stat_col, r_stat_col,
                              r_pcol, palm_paramp_col=[]):
    """
    Function to compare results between R and palm
    ----------------------------------------------

    Parameters:
    -----------
    data: pandas dataframe
        Dataframe with assembled R and palm results
    filename: str
        Parametric R results file name.
    palm_stat_col: list of str
        Palm stat value column names. Each column in the list will be
        compared to colum in r_stat_col at the same index.
    r_stat_col: list of str
        R stat value column names.
    palm_paramp_col: list of str
        Palm data eventual column with parametric pvalues. Each column in the
        list will be compared to colum in r_pcol at the same index.
    r_pcol: list of str
        R data pvalues colname.

    Returns:
    --------
    connections_checked_stats: list of str
        list of connections for which stat value were checked between R
        and Palm
    connections_checked_pvals: list of str
        list of connections for which pvalues were checked between R
        and Palm
    """

    connections_checked_stats = []
    connections_checked_pvals = []
    for conn in data.index:

        # Check stat values similarities
        connection_checked = 0
        for idx_col, rcol in enumerate(r_stat_col):

            # If R test is wilcoxon, welch or kruskall wallis test it is
            # useless to compare to palm test, skipping.
            # If permutation was not performed also skip.
            if "Test" in data.columns:
                if (data.loc[conn, "Test"] in [
                    "Wilcoxon rank sum test",
                    "Wilcoxon rank sum test with continuity correction",
                        "Welch Two Sample t-test", "Kruskal-Wallis"]):
                    continue
            if data.loc[conn, palm_stat_col[idx_col]] == "Not computed":
                continue

            r_stat_val = data.loc[conn, rcol]
            palm_stat_val = data.loc[conn, palm_stat_col[idx_col]]
            r_stat_val = Decimal(r_stat_val)
            r_stat_val = float(round(r_stat_val, 4))

            if r_stat_val != palm_stat_val:
                msg = "Different stat values file {0} conn {1}".format(
                    filename, conn)
                msg += " contrast {0} : ".format(palm_stat_col[idx_col])
                msg += "{0} vs {1}".format(str(r_stat_val), str(palm_stat_val))
                raise ValueError(msg)
            connection_checked += 1
        if connection_checked == len(r_stat_col):
            connections_checked_stats.append(conn)

        # Eventually check parametric pvalues similarities
        connection_checked = 0
        for idx_col, palmcol in enumerate(palm_paramp_col):

            # If R test is wilcoxon, welch or kruskall wallis test it is
            # useless to compare to palm test, skipping.
            # If permutation was not performed also skip.
            if "Test" in data.columns:
                if (data.loc[conn, "Test"] in [
                    "Wilcoxon rank sum test",
                    "Wilcoxon rank sum test with continuity correction",
                        "Welch Two Sample t-test", "Kruskal-Wallis"]):
                    continue
            if data.loc[conn, palmcol] == "Not computed":
                continue

            r_param_pval = data.loc[conn, r_pcol[idx_col]]
            palm_param_pval = data.loc[conn, palmcol]
            r_param_pval = Decimal(r_param_pval)
            r_param_pval = float(round(r_param_pval, 4))
            if r_param_pval != palm_param_pval:
                print(data.loc[conn, palm_stat_col[idx_col]])
                msg = "Different pvalues file {0} conn {1}".format(
                    filename, conn)
                msg += " contrast {0} : ".format(palmcol)
                msg += "{0} vs {1}".format(
                    str(r_param_pval), str(palm_param_pval))
                raise ValueError(msg)
            connection_checked += 1
        if connection_checked == len(palm_paramp_col):
            connections_checked_pvals.append(conn)

    return connections_checked_stats, connections_checked_pvals


def add_pval_col(
    param_data, perm_data, perm_pval_col, r_pval_col,
        distribution_col="Normality of residuals (Shapiro p-val > 0.05)"):
    """
    Function to compute add new pvalues column.
    -------------------------------------------

    Parameters:
    -----------
    param_data: pandas dataframe
        parametric data
    perm_data: pandas dataframe
        permutation data
    perm_pval_col: str
        column name for permutation pvalues.
    r_pval_col: str
        column name for parametric test (r) pvalues.
    distribution_col: str
        column name for normality of the test

    Returns:
    -------
    final_pvals: list of float
        list of assembled pvalues
    """
    final_pvals = []
    for conn in param_data.index:
        normal_distrib_conn = param_data.loc[conn, distribution_col]
        if normal_distrib_conn:
            final_pvals.append(param_data.loc[conn, r_pval_col])
        else:
            final_pvals.append(param_data.loc[conn, r_pval_col])
    return final_pvals


def fdrcorrection_intra_inter_networks(
    conns, pvalues,
        inter=[["Sal_", "DMN_"], ["Sal_", "ECN_"], ["DMN_", "ECN_"]]):

    """
    FDR correction by inter-intra networks
    --------------------------------------

    Parameters:
    -----------
    conns: list of str
        list of connections.
    pvalues: list of float
        list of pvalues for each connection.
    inter: list of str
        pattern in names of inter networks connections

    Returns:
    --------
    ordered_pvalues_list:
        ordered pvalues corrected by intra/internetwork fdr correction.
    """

    # Get all inter and intra networks connections
    inter_networks_conns = {}
    intra_networks_conns = {}
    for idx_conn, conn in enumerate(conns):
        inter_elts_idx = None
        for idx, inter_elts in enumerate(inter):
            elts_in_conn = 0
            for elt in inter_elts:
                if elt in conn:
                    elts_in_conn += 1
            if elts_in_conn == len(inter_elts):
                if inter_elts_idx is None:
                    inter_elts_idx = idx
                else:
                    msg = "Cant determine internetwork connection for conn {0}"
                    raise ValueError(msg.format(conn))
                inter_networks_conns[conn] = pvalues[idx_conn]
        if conn not in inter_networks_conns.keys():
            intra_networks_conns[conn] = pvalues[idx_conn]

    # Check expected number of connection
    # TODO : change these harcoded values
    if len(intra_networks_conns.keys()) != 408:
        raise ValueError(
            "Unexpected number of intranetworks connections : {0}".format(
                str(len(intra_networks_conns.keys()))))

    if len(inter_networks_conns.keys()) != 817:
        raise ValueError(
            "Unexpected number of intranetworks connections : {0}".format(
                str(len(inter_networks_conns.keys()))))

    # Correct pvalues for inter networks
    internetwork_conns = list(inter_networks_conns.keys())
    internetwork_conns_pvalues = []
    for conn in internetwork_conns:
        internetwork_conns_pvalues.append(inter_networks_conns[conn])
    rejected_inter_net, fdr_corr_pvalues_inter_net = fdrcorrection(
        pvals=internetwork_conns_pvalues,
        alpha=0.05,
        method='indep')

    # Same for intranetworks
    intranetwork_conns = list(intra_networks_conns.keys())
    intranetwork_conns_pvalues = []
    for conn in intranetwork_conns:
        intranetwork_conns_pvalues.append(intra_networks_conns[conn])
    rejected_intra_net, fdr_corr_pvalues_intra_net = fdrcorrection(
        pvals=intranetwork_conns_pvalues,
        alpha=0.05,
        method='indep')

    # Concatenante all pvalues together
    final_pvalues = {}
    for idx, conn in enumerate(internetwork_conns):
        final_pvalues[conn] = {
            "p": fdr_corr_pvalues_inter_net[idx],
            "rejected": rejected_inter_net[idx]}
    for idx, conn in enumerate(intranetwork_conns):
        final_pvalues[conn] = {
            "p": fdr_corr_pvalues_intra_net[idx],
            "rejected": rejected_intra_net[idx]}
    if len(final_pvalues.keys()) != 1225:
        raise ValueError(
            "Unexpected number of concatenated pvals : {0}".format(
                str(len(final_pvalues.keys()))))
    ordered_pvalues_list = []
    rejected_pvalues_list = []
    for conn in conns:
        ordered_pvalues_list.append(final_pvalues[conn]["p"])
        rejected_pvalues_list.append(final_pvalues[conn]["rejected"])

    return ordered_pvalues_list, rejected_pvalues_list


"""
End Functions
"""
# Get param + perm data
param_file = inputs["parametric_file"]
perm_file = inputs["permutation_file"]

# Read data
param_data = pd.read_csv(param_file, index_col=0)
perm_data = pd.read_csv(perm_file)

# Rename connection in perm data
perm_data.index = [x.replace("_to_", ".") for x in perm_data["Conn"]]

# Add all cols to each param file
perm_keep_cols = inputs["perm_keep_cols"]
perm_keep_cols_info = OrderedDict()
for conn in param_data.index:
    for col in perm_keep_cols:
        if col not in perm_keep_cols_info.keys():
            perm_keep_cols_info[col] = []

        # > Add not computed value
        if conn not in perm_data.index:
            perm_keep_cols_info[col].append("Not computed")
        else:
            perm_keep_cols_info[col].append(perm_data.loc[conn, col])
for col in sorted(perm_keep_cols):
    param_data[col] = perm_keep_cols_info[col]

# Check stat val for normal distrib + eventually pvals uncparap
stats_check_r = []
stats_check_palm = []
for both_cols in inputs["stat_check_col"]:
    both_cols = both_cols.split("_splithere_")
    r_col = both_cols[0]
    palm_col = both_cols[1]
    stats_check_r.append(r_col)
    stats_check_palm.append(palm_col)
pvals_check_r = []
pvals_check_palm = []

for both_cols in inputs["check_parametric_pval_col"]:
    both_cols = both_cols.split("_splithere_")
    r_col = both_cols[0]
    palm_col = both_cols[1]
    pvals_check_r.append(r_col)
    pvals_check_palm.append(palm_col)
conn_checked_stats, conn_checked_pvals = compare_palm_and_r_values(
    data=param_data,
    filename=param_file,
    palm_stat_col=stats_check_palm,
    r_stat_col=stats_check_r,
    palm_paramp_col=pvals_check_palm,
    r_pcol=pvals_check_r)
if len(conn_checked_stats) == 0:
    raise ValueError(
        "{0} vs {1} : difference between R and palm should be checked for "
        "at least 1 conn.".format(param_file, perm_file))

# Create col with pvalues from parametric and perm tests
assembled_pvals = add_pval_col(
    param_data=param_data,
    perm_data=perm_data,
    perm_pval_col=inputs["perm_pval_col"],
    r_pval_col=inputs["r_pval_col"],
    distribution_col="Normality of residuals (Shapiro p-val > 0.05)")
out_data = param_data
out_data[
    "Permutation pval (only for test with non normal data)"] = assembled_pvals

# FDR correction
rejected, fdr_corr_pvalues = fdrcorrection(
    pvals=assembled_pvals,
    alpha=0.05,
    method='indep')
out_data[
    "Corrected p-val (BH) (include permutation results)"] = fdr_corr_pvalues
col = "Survive FDR correction at q-FDR < 0.05 (include permutation results)"
out_data[col] = rejected

# Bonus : FDR corr by intra-inter network
(ordered_pvalues_list,
    rejected_pvalues_list) = fdrcorrection_intra_inter_networks(
        conns=list(out_data.index),
        pvalues=assembled_pvals,
        inter=[["Sal_", "DMN_"], ["Sal_", "ECN_"], ["DMN_", "ECN_"]])
colname = "Corrected p-val (BH-inter + BH-intra) (include permutation results)"
out_data[colname] = ordered_pvalues_list
colname = "Survive FDR correction at q-FDR < 0.05 (BH-inter + BH-intra) "
colname += "(include permutation results)"
out_data[colname] = rejected_pvalues_list

# Create csv output for analysis and timepoint and model
out_csv = os.path.join(
    inputs["outdir"], os.path.basename(param_file).replace(
        ".csv", "_plus_perm_results.csv"))
out_data.to_csv(out_csv, index=True)

# Add output file to json output log
outputs["output_file"] = out_csv
outputs["nb_conns_with_stat_doublechecked_with_perm"] = len(
    conn_checked_stats)
outputs["nb_conns_with_pvals_doublechecked_with_perm"] = len(
    conn_checked_pvals)


"""
Write logs
"""
# Write inputs/outputs logs
logdir = os.path.join(inputs["outdir"], "logs")
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
