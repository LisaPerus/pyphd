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
from pyphd.constants import SCRIPT_DIR
from statsmodels.stats.multitest import fdrcorrection

# Script documentation
DOC = """
Correct for multiple intervention to placebo testing mapt rsfmri analysis
-------------------------------------------------------------------------

/!\ Script adapted only to MAPT rsfmri output files /!\

Test on dummy data:

echo "Conn,Permutation pval (only for test with non normal data)">> /tmp/test1.txt
echo "VDMN7.VDMN5,0.005" >> /tmp/test1.txt
echo "VDMN7.ASal3,0.1" >> /tmp/test1.txt
echo "RECN3.RECN1,0.03" >> /tmp/test1.txt

echo "Conn,Permutation pval (only for test with non normal data)">> /tmp/test2.txt
echo "VDMN7.VDMN5,0.0001" >> /tmp/test2.txt
echo "VDMN7.ASal3,0.18" >> /tmp/test2.txt
echo "RECN3.RECN1,0.3" >> /tmp/test2.txt

echo "Conn,Permutation pval (only for test with non normal data)">> /tmp/test3.txt
echo "VDMN7.VDMN5,0.2" >> /tmp/test3.txt
echo "VDMN7.ASal3,0.001" >> /tmp/test3.txt
echo "RECN3.RECN1,0.007" >> /tmp/test3.txt


python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_mtesting_corr_interventions.py \
    -i /tmp/test1.txt /tmp/test2.txt /tmp/test3.txt \
    -t test1 test2 test3 \
    -w test \
    -o /tmp

Test on MAPT data:

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_rsfmri_mtesting_corr_interventions.py \
    -i $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_PLUS_PERMUTATION_TESTS/M0/ttest_conn_connectivities_CDR0_only_IM_1_Placebo_0_0399_rscores_common_subjects_timesteps_PREPOST_at_PRE_plus_perm_results.csv \
    $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_PLUS_PERMUTATION_TESTS/M0/ttest_conn_connectivities_CDR0_only_Omega3_1_Placebo_0_0399_rscores_common_subjects_timesteps_PREPOST_at_PRE_plus_perm_results.csv \
    $MEDIA_SCRIPT/MAPT_rsfmri/MONTPELLIER_ONLY/TESTS_WITH_SUBJECTS_COMMON_AT_M0_M36/PARAMETRIC_PLUS_PERMUTATION_TESTS/M0/ttest_conn_connectivities_CDR0_only_Omega3_plus_IM_1_Placebo_0_0399_rscores_common_subjects_timesteps_PREPOST_at_PRE_plus_perm_results.csv \
    -t IM O3 O3pIM \
    -w CDR0_IntCorr_vs_Pl \
    -o /tmp/test_dir

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
        prog="python mapt_rsfmri_mtesting_corr_interventions.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-files", type=is_file, required=True, nargs="+",
        help="Input files for each test that has to be corrected.")
    required.add_argument(
        "-t", "--test-names", type=str, required=True, nargs="+",
        help="Each test name corresponding to each input file.")
    required.add_argument(
        "-w", "--whole-testname", type=str, required=True,
        help=".")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

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
    "tool": "mapt_rsfmri_mtesting_corr_interventions.py"
}
outputs = {}
if verbose > 0:
    pprint(
        "[info] Correct MAPT interventions with Hochberg (1988) correction...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Get all connections
all_conns = []
for fid in inputs["input_files"]:
    data = pd.read_csv(fid, index_col=0)
    all_conns += list(data.index)
all_conns = list(set(all_conns))

pvalues_data = {}
for idx, fid in enumerate(inputs["input_files"]):

    data = pd.read_csv(fid, index_col=0)

    # Get pvalues
    for conn in all_conns:

        if conn not in pvalues_data.keys():
            pvalues_data[conn] = {}

        pvalues_data[conn][fid] = data.loc[
            conn, "Permutation pval (only for test with non normal data)"]

# Get corrected pvalues
corrected_pvalues_data = {}
test_dir = os.path.join(inputs["outdir"], inputs["whole_testname"])
conn_outdir = os.path.join(test_dir, "conns_corrected_pvalues")
if not os.path.isdir(conn_outdir):
    os.makedirs(conn_outdir)
for conn in all_conns:

    # Correct pvalues
    cmd = ["Rscript", "--vanilla",
                      os.path.join(SCRIPT_DIR, "GIT_REPOS", "RPhd", "scripts",
                                   "hochberg_correction.R")]
    conn_pvals = []
    conn_tests = []
    for idx, fid in enumerate(inputs["input_files"]):
        test_name = inputs["test_names"][idx]
        conn_pvals.append(str(pvalues_data[conn][fid]))
        conn_tests.append(test_name)
    outfile = os.path.join(conn_outdir, conn.replace(".", "_to_") + ".csv")
    cmd = cmd + ["-p"] + [" ".join(conn_pvals)]
    cmd += ["-n"] + [" ".join(conn_tests)]
    cmd += ["-o"] + [outfile]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1}".format(" ".join(cmd), stderr))

    # Get corrected pvalues
    conn_results = pd.read_csv(outfile)
    if conn_results.shape[0] != 1:
        raise ValueError("Not one line for {0}".format(outfile))

    for idx, fid in enumerate(inputs["input_files"]):
        test_name = inputs["test_names"][idx]
        if fid not in corrected_pvalues_data.keys():
            corrected_pvalues_data[fid] = {}
        corrected_pvalues_data[fid][conn] = conn_results.loc[0, test_name]

# Write results
for idx, fid in enumerate(inputs["input_files"]):
    outdata = pd.read_csv(fid, index_col=0)
    test_name = inputs["test_names"][idx]
    corrected_hochberg_pvals_col = []
    for conn in outdata.index:
        corrected_hochberg_pvals_col.append(corrected_pvalues_data[fid][conn])
    outdata["Permutation pval (only for test with non normal data) corrected for intervention (Hochberg)"] = corrected_hochberg_pvals_col

    # Apply FDR correction
    rejected, fdr_corr_pvalues = fdrcorrection(
        pvals=outdata["Permutation pval (only for test with non normal data) corrected for intervention (Hochberg)"],
        alpha=0.05,
        method='indep')
    outdata["Corrected p-val (BH) (include permutation results + intervention corr)"] = fdr_corr_pvalues
    outdata["Survive FDR correction at q-FDR < 0.05 (include permutation results + intervention corr)"] = rejected

    # Save file
    _, outfile_ext = os.path.splitext(fid)
    outfile = os.path.join(
        inputs["outdir"],
        os.path.basename(fid).replace(
            outfile_ext, "_corr_intervention" + outfile_ext))
    outdata.to_csv(outfile, index=True)
    outputs["{0} stat file with corrected pvals".format(test_name)] = outfile


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(os.path.join(inputs["outdir"], "logs"))
if not os.path.isdir(logdir):
    os.mkdir(logdir)
script_name = "mapt_rsfmri_mtesting_corr_interventions"
for name, final_struct in [
    ("{0}_{1}_inputs".format(script_name, inputs["whole_testname"]), inputs),
    ("{0}_{1}_outputs".format(script_name, inputs["whole_testname"]), outputs),
    ("{0}_{1}_runtime".format(script_name, inputs["whole_testname"]), runtime)
        ]:
    log_file = os.path.join(logdir, "{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[info] Outputs:")
    pprint(outputs)
