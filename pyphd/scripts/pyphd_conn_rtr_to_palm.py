#!/usr/bin/env python
# Script used to extract ROI-to-ROI Conn 1rst level analysis connectivity
# scores and input them to FSL PALM for permutation testing.
# Copyright (C) 2019  Lisa Perus
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
from datetime import datetime
from pprint import pprint
import progressbar

# Third party imports
import pandas as pd
import numpy as np
from numpy import tanh
from scipy.io import loadmat
from statsmodels.stats.multitest import fdrcorrection

# Pyphd imports
from pyphd.fsl.utils import palm, text2vest

# Script documentation
DOC = """
Script used to extract ROI-to-ROI Conn 1rst level analysis connectivity
scores and input them to FSL PALM for permutation testing.
-----------------------------------------------------------------------

1) Extract Z-score matrix
2) Compute r-scores from Z-score matrix
3) Use permutation testing to get some difference between groups, etc.
4) Perform multiple testing correction (BH)

Example on MAPT data:
cat input.csv
Sid,Conn_file,Design_group1,Design_group2
Subject001,resultsROI_Subject001_Condition001.mat,1,0
Subject013,resultsROI_Subject013_Condition001.mat,0,1
...
Subject00n,resultsROI_Subject0n_Condition001.mat,0,1

python3 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/pyphd_conn_rtr_to_palm.py \
    -i /tmp/test_palm/input.csv \
    -c $HOME/PHD/DATA/RESULTS/MAPT/CONNECTIVITY_ANALYSIS/SBC_05_Cambridge_R7/stats/one_sided_ttest.txt \
    -N 10000 \
    -f $HOME/fsl_init.sh \
    -o /tmp/test_palm \
    -V 1

NB: for an ANOVA give t-contrast file (4gpes):
1 0 0 -1
0 1 0 -1
0 0 1 -1
and f contrast file:
1 1 1
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
        prog="python pyphd_conn_rtr_to_palm.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-file", type=is_file, required=True,
        help="Csv file listing for each subject its conn first level "
             "condition file and design matrix. Has a header.")
    required.add_argument(
        "-c", "--contrast", type=is_file, help="Contrast file (text file).")
    required.add_argument(
        "-f", "--fsl-config", metavar="<path>", type=is_file,
        help="Path to fsl sh config file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-O", "--output-name", type=str, help="Output name.")
    parser.add_argument(
        "-D", "--demean", type=int, nargs="+",
        help="Demean covariates in col 0, 1, ..., n IN DESIGN FILE as "
             "suggested by FSL see "
             "(https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM/Faq) "
             "/!\ First column has a 0 index and corresponds to column "
             "index is design file /!\ ")
    parser.add_argument(
        "-M", "--multiple-testing-constrasts", type=str, nargs="+",
        help="Correct for multiple testing on these contrasts. E.g : c1, c2.")
    parser.add_argument(
        "-C", "--connections", type=str, nargs="+",
        help="Connections that are to be kept only."
             "Each element of a connection must be separated by ':' "
             "(e.g DMN:Visual)")
    parser.add_argument(
        "-F", "--f-contrast", type=is_file,
        help="F contrast txt file. If set :  test the t contrasts in a single "
             "F contrast.")
    parser.add_argument(
        "-T", "--two-tail", action="store_true",
        help="If set, perform two-tailed test.")
    parser.add_argument(
        "-N", "--nb-permutations", type=int, default=10000,
        help="Number of permutation for permutation testing.")
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
    "tool": "pyphd_conn_rtr_to_palm.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Conn connectivity scores testing...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


"""
Step 0 - Load data
"""
conn_mat_files = []
subjects = []
designs = []
nan_subjects = []  # Record subjects with no information about the design
nan_subjects_sids = []
with open(inputs["input_file"], "rt") as open_file:
    lines = open_file.readlines()
input_data_path, input_data_ext = os.path.splitext(inputs["input_file"])
header = lines[0].strip("\n").split(",")
for idx, line in enumerate(lines[1:]):
    line = line.strip("\n").split(",")
    design = ",".join(line[2:])
    if ("NaN" in design or "NA" in design or "na" in design
        or "" in [x.strip("\n").strip(" ") for x in design.split(",")]):
        nan_subjects.append(idx)
        nan_subjects_sids.append(line[0])
    subjects.append(line[0])
    conn_mat_files.append(line[1])
    designs.append(design)

"""
Step 1 : Get connectivity scores
"""
r_scores_connectivity_matrices = []
z_scores_connectivity_matrices = []
for conn_mat_file in conn_mat_files:
    conn_data = loadmat(conn_mat_file)

    # Create z-score dataframe
    # > Get row and column names
    source_names = [x[0] for x in conn_data["names"][0]]
    target_names = [x[0] for x in conn_data["names2"][0]]
    zscores_df = pd.DataFrame(
        conn_data["Z"], index=source_names, columns=target_names)

    # Compute r-scores
    rscores_df = zscores_df.apply(tanh)

    r_scores_connectivity_matrices.append(rscores_df)
    z_scores_connectivity_matrices.append(zscores_df)


"""
Step 2 : Do permutation testing for each connection
"""
# Get connections
if inputs["connections"] is None:
    connections = []
    for matrix in z_scores_connectivity_matrices:
        for row in matrix.index:
            for col in matrix.columns:
                connections.append(row + ":" + col)

    # Get unique connections
    connections = np.unique(np.array(connections))
else:
    connections = inputs["connections"]

# Keep only connections one way
delete_connections = []
for idx, conn in enumerate(connections):
    if conn in delete_connections:
        continue
    conn = conn.split(":")
    if len(conn) > 2:
        print("Multiple ':' in connection {0} -> ambiguous splitting.".format(
              ":".join(conn)))
    opp_conn = conn[1] + ":" + conn[0]
    if opp_conn in connections:
        delete_connections.append(opp_conn)
connections = [x for x in connections if x not in delete_connections]

# Write design mat file
design_file = os.path.join(
    inputs["outdir"], "_".join(header[2:]) + "_design.txt")
with open(design_file, "wt") as open_file:
    for idx, design in enumerate(designs):
        if idx in nan_subjects:
            continue
        design = design.strip("\n").split(",")
        open_file.write("\t".join(design))
        open_file.write("\n")

# If needed mean center columns
if inputs["demean"] is not None:
    design_data = pd.read_csv(design_file, sep="\t", header=None)
    for col in inputs["demean"]:
        col_values = design_data.iloc[:, col].astype(float)
        mean_val_col = np.mean(col_values)
        mean_centered_col_values = col_values - mean_val_col
        design_data.iloc[:, col] = mean_centered_col_values
    design_data.to_csv(design_file, sep='\t', header=False, index=False)

design_mat_file = design_file.replace(".txt", ".mat")
text2vest(
    indata=design_file,
    outdata=design_mat_file,
    fsl_sh=inputs["fsl_config"])

# If contrast is in text file, convert in .con file
if inputs["contrast"].endswith(".txt"):
    constrat_file = os.path.join(
        inputs["outdir"],
        os.path.basename(inputs["contrast"]).replace(".txt", ".con"))
    text2vest(
        indata=inputs["contrast"],
        outdata=constrat_file,
        fsl_sh=inputs["fsl_config"])
else:
    constrat_file = inputs["contrast"]

# If f contrast text file is specified, convert it to fts file
if inputs["f_contrast"] is not None:
    f_contrast_file = os.path.join(
        inputs["outdir"],
        os.path.basename(inputs["f_contrast"]).replace(".txt", ".fts"))
    text2vest(
        indata=inputs["f_contrast"],
        outdata=f_contrast_file,
        fsl_sh=inputs["fsl_config"])
else:
    f_contrast_file = None

# Write connections values
connectivity_dir = os.path.join(
    inputs["outdir"], os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "") + "_connections")
if not os.path.isdir(connectivity_dir):
    os.mkdir(connectivity_dir)

palm_outdir = os.path.join(
    inputs["outdir"], os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "") + "_palm_output")
if not os.path.isdir(palm_outdir):
    os.mkdir(palm_outdir)
outputs["palm_outdir"] = palm_outdir
with progressbar.ProgressBar(max_value=len(connections),
                             redirect_stdout=True) as bar:
    for idx_conn, connection in enumerate(connections):
        print("Processing {0}...".format(connection))
        connections_values = []
        source, target = connection.split(":")
        source = source
        target = target
        for idx, matrix in enumerate(z_scores_connectivity_matrices):
            sub_connection = matrix.loc[source, target]
            if idx not in nan_subjects:
                connections_values.append(sub_connection)

        # Write connections file
        connection_file = os.path.join(
            connectivity_dir, "{0}_to_{1}_connectivity_values.csv".format(
                source.replace(" ", ""), target.replace(" ", "")))
        with open(connection_file, "wt") as open_file:
            for val in connections_values:
                open_file.write(str(val))
                open_file.write("\n")
        outputs["{0}_to_{1}_connectivity_values".format(
                source, target)] = connection_file

        # Run Palm
        palm_output_basename = os.path.join(
            palm_outdir, "{0}_to_{1}_palm".format(
                source.replace(" ", ""), target.replace(" ", "")))
        if inputs["two_tail"]:
            palm_output_basename += "two_tail"
        stat_val, pval_unc, p_fwe = palm(
            indata=connection_file,
            design_file=design_mat_file,
            contrast_file=constrat_file,
            f_contrast=f_contrast_file,
            output_basename=palm_output_basename,
            nb_permutations=inputs["nb_permutations"],
            twotail=inputs["two_tail"])
        outputs["{0}_to_{1}_palm_stat_val".format(source, target)] = stat_val
        outputs["{0}_to_{1}_palm_pval_unc".format(source, target)] = pval_unc
        outputs["{0}_to_{1}_palm_p_fwe".format(source, target)] = p_fwe
        bar.update(idx_conn)


"""
Step 4 : Perform multiple testing correction
"""
if inputs["multiple_testing_constrasts"] is not None:
    contrasts_results = {}
    for contrast in inputs["multiple_testing_constrasts"]:
        contrasts_results[contrast] = {"tvalues": [], "unc_pvalues": [],
                                       "fdr_corrected_pvalues": [],
                                       "rejected_fdr_corr_pvals_p005": []}

        # Get all contrast tvals or fstats
        stats_files = []
        for conn in connections:
            stat_file = glob.glob(
                os.path.join(palm_outdir, "{0}*stat_{1}.csv".format(
                    conn.replace(":", "_to_"), contrast)))
            if len(stat_file) == 0:
                print("Could not find stat file for {0}".format(conn))
            else:
                stat_file = stat_file[0]
                stats_files.append(stat_file)

        uncp_files = []
        t_stat = True
        for stat_file in stat_files:
            if "tstat" in stat_file:
                uncp_file = stat_file.replace("tstat", "tstat_uncp")
            if "fstat" in stat_file:
                uncp_file = stat_file.replace("fstat", "fstat_uncp")
                t_val = False
            uncp_files.append(uncp_file)

        stat_values = []
        for fid in stat_files:
            with open(fid, "rt") as open_file:
                lines = open_file.readlines()
            val = float(lines[0].strip("\n"))
            stat_values.append(val)

        unc_pvalues = []
        for fid in uncp_files:
            with open(fid, "rt") as open_file:
                lines = open_file.readlines()
            unc_pval = float(lines[0].strip("\n"))
            unc_pvalues.append(unc_pval)

        if t_stat:
            contrasts_results[contrast]["tvalues"] = stat_values
        else:
            contrasts_results[contrast]["fvalues"] = stat_values
        contrasts_results[contrast]["unc_pvalues"] = unc_pvalues

        # Apply Benjamini-Hochberg FDR correction
        # FDR correction
        rejected, fdr_corr_pvalues = fdrcorrection(
            pvals=unc_pvalues,
            alpha=0.05,
            method='indep')
        contrasts_results[contrast]["fdr_corrected_pvalues"] = fdr_corr_pvalues
        contrasts_results[contrast]["rejected_fdr_corr_pvals_p005"] = rejected

        contrast_summary_file = os.path.join(
            inputs["outdir"], os.path.basename(inputs["input_file"]).replace(
                input_data_ext, "") + "{0}_palm_results.csv".format(contrast))
        with open(contrast_summary_file, "wt") as open_file:
            if t_stat:
                open_file.write("Connection,Tstat,Pval_unc,Pval_FDR_corr,")
            else:
                open_file.write("Connection,Fstat,Pval_unc,Pval_FDR_corr,")
            open_file.write("Reject_H0_FDR\n")
            for idx, conn in enumerate(connections_names):
                line = ",".join([conn, str(stat_values[idx]),
                                 str(unc_pvalues[idx]),
                                 str(fdr_corr_pvalues[idx]),
                                 str(rejected[idx])])
                line += "\n"
                open_file.write(line)
        outputs["Summary_contrast_{0}".format(
                contrast)] = contrast_summary_file


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
if inputs["two_tail"]:
    output_basename = os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "_two_tail_palm")
else:
    output_basename = os.path.basename(inputs["input_file"]).replace(
        input_data_ext, "_one_tail_palm")
if inputs["output_name"] is not None:
    output_basename += "_" + inputs["output_name"]
nan_subjects_file = os.path.join(
    inputs["outdir"], output_basename + "_nan_subjects.csv")
with open(nan_subjects_file, "wt") as open_file:
    open_file.write("Discarded subjects ")
    open_file.write("(with NaN value for at least one variable\n")
    for sid in nan_subjects_sids:
        open_file.write(sid)
        open_file.write("\n")
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
