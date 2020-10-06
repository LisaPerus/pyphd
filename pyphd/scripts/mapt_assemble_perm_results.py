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

# Script documentation
DOC = """

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/mapt_assemble_perm_results.py \
    -p CDR_x_gpeMapt4c_palm_output \
    -o /tmp

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
        prog="python mapt_assemble_perm_results.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-p", "--permutation-dir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the permutation directory.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional arguments
    parser.add_argument(
        "-S", "--split-pattern", type=str, default="_connectivity_values_",
        help="Pattern to split files 'type of file' name and beginning of "
             "filename. 'Type of file' name is always in the second part of "
             "the whole filename.")
    parser.add_argument(
        "-I", "--ignore-file-pattern", type=str, nargs="+",
        default=["palm_elapsed.csv", "palmtwo_tail_elapsed.csv"],
        help="Ignore files with specific patterns.")
    parser.add_argument(
        "-A", "--stack-file-pattern", type=str, nargs="+",
        default=["palm_palmconfig.txt", "palmtwo_tail_palmconfig.txt"],
        help="Stack content of files with specific patterns.")
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
    pprint("[info] Assemble rsfmri perm data...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)


# List all files in permutation outdir
all_files = glob.glob(os.path.join(inputs["permutation_dir"], "*"))
all_files = [x for x in all_files if os.path.isfile(x)]
files_pattern = [x.split(inputs["split_pattern"])[1] for x in all_files]
files_pattern = set(files_pattern)

# Get files that have to be assembled separetely
files_pattern = [
    x for x in files_pattern if x not in inputs["ignore_file_pattern"] and
    x not in inputs["stack_file_pattern"]]

# For each file, read file value
results = {}
all_connections = []
for fid in all_files:
    fid_pat = fid.split(inputs["split_pattern"])
    conn = os.path.basename(fid_pat[0])
    fid_pat = fid_pat[1]
    if fid_pat not in files_pattern:
        continue

    # > Read value
    with open(fid, "rt") as open_file:
        lines = open_file.readlines()
    if len(lines) > 1:
        raise ValueError(
            "Permutation file {0} is expected to have only one line".format(
                fid))

    # > Store value
    value = lines[0].strip("\n")
    if fid_pat not in results.keys():
        results[fid_pat] = {}
    results[fid_pat][conn] = value
    all_connections.append(conn)
all_connections = set(all_connections)

# Write stack result
outputs["data_stacked"] = {}
for pat in inputs["stack_file_pattern"]:
    stat_result_output = os.path.join(
        inputs["outdir"], "all_perm_stacked_" + pat)
    files_pat = glob.glob(
        os.path.join(inputs["permutation_dir"], "*" + pat + "*"))
    with open(stat_result_output, "wt") as open_file:
        for fid in files_pat:
            with open(fid, "rt") as palm_file:
                lines = palm_file.readlines()
            open_file.write(fid + "\n")
            for line in lines:
                open_file.write(line)
    outputs["data_stacked"][pat] = stat_result_output

# Write output file
output_file = os.path.join(inputs["outdir"], "perm_results_assembled.csv")
if os.path.isfile(output_file):
    raise ValueError(
        "File {0} already exist, cannot overwrite.".format(output_file))
sorted_patterns = sorted(results.keys())

with open(output_file, "wt") as open_file:
    header = "Conn"
    header += "," + ",".join([x.replace(".csv", "") for x in sorted_patterns])
    open_file.write(header)
    open_file.write("\n")
    for conn in all_connections:
        line = [conn]
        for pat in sorted_patterns:
            line.append(results[pat][conn])
        open_file.write(",".join(line))
        open_file.write("\n")
outputs["Assembled perm results"] = output_file


"""
Update the outputs and save them and the inputs in a 'logs' directory.
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
output_basename = "mapt_assemble_perm_results"
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(logdir, output_basename + "_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)
if verbose > 0:
    print("[final]")
    pprint(outputs)
