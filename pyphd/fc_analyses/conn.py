# Functions to use Conn outputs
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

# System imports
import os
import glob
import argparse

# Third party imports
from pyphd.constants import CONN_INPUTS
import pandas as pd
from scipy.io import loadmat


def extract_connectivities(group_name, tp=None, center_name=None,
                           covariates=None, network=None,
                           conn_file_pattern=None, datafile=None, outdir=None,
                           conn_datapath=None, tp_name=None):
    """Extract connectivities values from Conn mat files.

    Parameters
    ----------
    group_name: str
        Add information from column group_name for each subject.
        Must be the same column name as in datafile.
    tp: str
        Timepoint
    center_name : str
        Center name
    covariates : list of str
        list of covariates. Names must be the same as those in datafile
        columns.
    network: str
        Network name. If set to None, select all inter/intra networks
        connections, else select only intra-network connections.
    conn_file_pattern: str
        Conn mat files pattern. If set to None with tp not set to None,
        select CONN_INPUTS default pattern.
    datafile : str
        path to csv datafile listing subjects clinical characteristics.
        Subjects must be in the same order as the Conn mat files.
        E.g : subject1 -> conn_mat_file_1.mat
        If set to None with tp not set to None, select CONN_INPUTS default
        datafile.
    outdir : str
        path to output directory. If set to None with tp not set to None,
        select CONN_INPUTS default outdir.
    conn_datapath : str
        path to conn data. If set to None with tp not set to None,
        select CONN_INPUTS default conn datapath.
    tp_name: str
        Name added to output file. If set to None with tp not set to None,
        select CONN_INPUTS default tp_name.

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # If conn_file_pattern, datafile, outdir, conn_datapath or tp_name are set
    # to None, select default values from CONN_INPUTS (tp must be specified!)
    if conn_file_pattern is None:
        if tp is not None:
            conn_file_pattern = CONN_INPUTS[tp]["conn_file_pattern"]
        else:
            raise ValueError(
                "Please specify a pattern for conn_file_pattern, or set a"
                " timepoint for default values.")
    if datafile is None:
        if tp is not None:
            datafile = CONN_INPUTS[tp]["datafile"]
        else:
            raise ValueError(
                "Please specify a path for datafile, or set a timepoint"
                " for default values.")
    if outdir is None:
        if tp is not None:
            outdir = CONN_INPUTS[tp]["outdir"]
        else:
            raise ValueError(
                "Please specify a path for outdir, or set a timepoint"
                " for default values.")
    if conn_datapath is None:
        if tp is not None:
            conn_datapath = CONN_INPUTS[tp]["conn_datapath"]
        else:
            raise ValueError(
                "Please specify a path for conn_datapath, or set a timepoint"
                " for default values.")
    if tp_name is None:
        if tp is not None:
            tp_name = CONN_INPUTS[tp]["timepoint_name"]
        else:
            raise ValueError(
                "Please specify a path for tp_name, or set a timepoint"
                " for default values.")

    # Read data
    with open(datafile, "rt") as open_file:
        lines = open_file.readlines()
    header = lines[0].split(",")

    # Get group col id
    if group_name is not None:

        # Find group id
        gpe_idx = None
        for idx, elt in enumerate(header):
            if elt.strip("\n") == group_name:
                gpe_idx = idx
        if gpe_idx is None:
            raise ValueError("Could not find group {0}".format(group_name))

    # Keep only center subjects
    keep_subjects = []
    if center_name is not None:
        for idx, line in enumerate(lines[1:]):
            sid = line.split(",")[0]
            center = sid.replace("sub-", "")[0:4]
            if center == center_name:
                keep_subjects.append(idx)  # Begins at 0
    else:
        keep_subjects = [x for x in range(0, len(lines) - 1)]

    # Keep only files of interest
    # It is assumed that Conn file when sorted are in the same order that the
    # subjects in datafile. Code will not work if this is not the case.
    conn_files = sorted(
        glob.glob(os.path.join(conn_datapath, conn_file_pattern)))

    # Get all connections associated with a network
    conn_data = loadmat(conn_files[0])
    all_networks = [x[0] for x in conn_data["names"][0]]
    target_networks = [x[0] for x in conn_data["names2"][0]]
    connections = []
    if network is not None:
        source_networks = [x for x in all_networks if network in x]
    else:
        source_networks = all_networks

    for source in source_networks:
        for target in source_networks:
            if source == target:
                continue
            if [target, source] in connections:
                continue
            else:
                connections.append([source, target])
    print("{0} connections to evaluate...".format(len(connections)))

    # Add covariates if necessary
    covariates_results = {}
    covariates_col_ids = {}
    if covariates is not None:
        for cov in covariates:
            covariates_results[cov] = []
        for cov in covariates:
            idx_col = None
            for idx, elt in enumerate(header):
                if elt.strip("\n").strip(" ") == cov:
                    covariates_col_ids[cov] = idx
                    idx_col = idx
            if idx_col is None:
                print("Could not find {0} in header".format(cov))

    # Get for each connection the connectivity scores of all subjects
    conn_results = {}
    groups_order = []
    for idx, mat_file in enumerate(conn_files):
        if idx not in keep_subjects:
            continue
        conn_data = loadmat(mat_file)
        zscores_df = pd.DataFrame(
            conn_data["Z"], index=all_networks,
            columns=target_networks)
        if group_name is None:
            gpe = "1"
        else:
            gpe = lines[1:][idx].split(",")[gpe_idx].strip("\n")
        groups_order.append(gpe)
        for conn in connections:
            source = conn[0]
            tg = conn[1]
            zscore = zscores_df.loc[source, tg]
            if ":".join(conn) not in conn_results.keys():
                conn_results[":".join(conn)] = [str(zscore)]
            else:
                conn_results[":".join(conn)].append(str(zscore))

        # Get covariate if necessary
        if covariates is not None:
            for cov in covariates:
                sid_cov_data = lines[1:][idx].strip("\n").split(",")[
                                    covariates_col_ids[cov]]
                covariates_results[cov].append(sid_cov_data)

    # Write results
    outfile = os.path.join(outdir, "conn_connectivities.txt")
    if group_name is not None:
        outfile = outfile.replace(".txt", "_" + group_name + ".txt")
    if center_name is not None:
        outfile = outfile.replace(".txt", "_" + center_name + ".txt")
    if network is not None:
        outfile = outfile.replace(".txt", "_" + network + ".txt")
    if tp is not None:
        outfile = outfile.replace(".txt", "_" + tp_name + ".txt")
    if covariates is not None:
        outfile = outfile.replace(".txt", "_".join(covariates) + ".txt")
    all_conns = sorted(conn_results.keys())
    with open(outfile, "wt") as open_file:
        if group_name is None:
            open_file.write("GROUP")
        else:
            open_file.write(group_name)

        # Write header
        for conn in all_conns:
            open_file.write("," + conn)

        # Add covariates if necessary
        if covariates is not None:
            for cov in covariates:

                # Write covariate and convert - to _ char (for R script)
                open_file.write("," + cov.replace("-", "_"))
        open_file.write("\n")

        # Write connection for each subject
        for idx, gpe in enumerate(groups_order):
            open_file.write(gpe)
            for conn in all_conns:
                open_file.write("," + conn_results[conn][idx])

            # Add covariates if necessary
            if covariates is not None:
                for cov in covariates:
                    open_file.write("," + covariates_results[cov][idx])
            open_file.write("\n")
    return outfile


def extract_group(conn_file, groups_info={}, rename_cols={}, erase_cols=[],
                  rename_file={}):
    """Extract from connectivities values file group and subgroups of subjects.

    /!\ conn_file MUST HAVE A ONE LINE HEADER /!\

    Parameters
    ----------
    conn_file: str
        File with connectivities values for each subject : one subject per
        line and the connections on the columns.
                    dmn1 <-> dmn2, dmn3 <-> dmn4
        E.g : Sid1        0.9           0.2
        Usually outputted by extract_connectivities function.
    groups_info : dict
        Dictionnary with information to extract groups. Keys are group columns,
        and values are dictionnary with values to keep and new values to
        replace them with.
        E.g : {0 : {"4:ctrl" : "0", "3:IM" : "1"}, -1 : {"1" : "1"}}
        will keep subjects with values 4:ctrl and 3:IM in the first column (
        and replace these values by 0 and 1) but also have value 1 in the last
        column.
    rename_cols : dict
        Keys are col ids in file, and values are new columns name in output
        file.
    erase_cols : list of str
        List of cols to delete at the end (by column name).
    rename_file: dict
        dict which rename parts of infile to create outfile.
        E.g : infile = "connectitivies_gpeMapt4c.txt"
              rename_file = {"gpeMapt4c" : "IM_1_Placebo_0"}
              outfile = "connectitivies_IM_1_Placebo_0"

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # Read conn file
    conn_data = pd.read_csv(conn_file, header=0, dtype=str)
    outdata = conn_data.copy()

    # Select subjects in subgroups
    for col, col_data in groups_info.items():
        outdata = outdata[outdata.iloc[:, col].isin(list(col_data.keys()))]

    # Rename subgroup subjects with new values
    for col, col_data in groups_info.items():
        for old_name, new_name in col_data.items():
            outdata.iloc[:, col][outdata.iloc[:, col] == old_name] = new_name

    # Erase columns
    for col_to_erase in erase_cols:
        del outdata[col_to_erase]

    # Rename columns
    new_colnames = []
    for old_name, new_name in rename_cols.items():
        for idx, col in enumerate(outdata.columns):
            if col == old_name:
                new_colnames.append(new_name)
            else:
                new_colnames.append(col)
    outdata.columns = new_colnames

    # Create outfile
    outfile = conn_file
    for elt_old, elt_new in rename_file.items():
        outfile = outfile.replace(elt_old, elt_new)

    # Save output file
    outdata.to_csv(outfile, index=False)
    return outfile


def sort_connectivities_by_pvals(conn_result_file, p_val_thresh,
                                 p_val_colname, additional_infos_columns):
    """Extract from a statistic test result file connectitivies under
       a specific pval threshold.

    Parameters
    ----------
    conn_result_file: str
        Path to file with statistics results for each connection.
        Must have at least one column for pvalues.
        Connectivities are supposed to be in the first column.
    p_val_thresh : float
        pval threshold. Can be higher than 0.05 if we want to observe trends
        (especially if a correction has been applied).
    p_val_colname : str
        colname containing pvalues on which the pvalues have to be thresholded.
    additional_infos_columns: dict
        Indicates supplementary info to save for significant connections.
        Key is column number, value is the name of var of interest.

    Returns
    -------
    results: List of dict
        List of dictionnaries with, for each dictionnary
    """

    # Read infile
    data = pd.read_csv(conn_result_file)

    # Select connections under threshold
    significant_data = data[data[p_val_colname] < p_val_thresh]
    significant_data = significant_data.sort_values(by=[p_val_colname])

    # Store results
    results = []
    for idx, row in significant_data.iterrows():
        results_conn = {}
        for col, info_new_name in additional_infos_columns.items():
            results_conn[info_new_name] = row[col]
        results[conn] = results_conn
    return results
