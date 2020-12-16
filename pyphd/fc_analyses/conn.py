# Functions to use Conn outputs and general connectivity data
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
import numpy as np
from pyphd.constants import CONN_INPUTS
from pyphd.fsl.utils import text2vest
import pandas as pd
from scipy.io import loadmat
from statsmodels.stats.multitest import fdrcorrection


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
                  rename_file={}, create_interaction_col=False,
                  interaction_cols=[]):
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
        and replace these values by 0 and 1) but and subjects with value 1 in
        the last column.
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
    create_interaction_col: bool
        if True, create a column with caracteristics from multiple cols
        specified in interaction_cols. For example if
        interaction_cols = ["subGpe", "mainGpe"] with
        mainGpe    subGpe
          gpe1      s1
          gpe2      s2
          gpe3      s2
        a new column subGpe_mainGpe will be created:
        subGpe_mainGpe
            s1_gpe1
            s2_gpe2
            s2_gpe3
    /!\ WARNING : this step is done before renaming columns /!\
    interaction_cols : list of str
        list of columns used to create interaction columns.

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # Read conn file
    conn_data = pd.read_csv(conn_file, header=0, dtype=str)
    outdata = conn_data.copy()

    # Select subjects in subgroups
    # Note 30/09/2020 : if there are nan values in subgroups it may delete
    # some lines
    for col, col_data in groups_info.items():
        outdata = outdata[outdata.iloc[:, col].isin(list(col_data.keys()))]

    # Rename subgroup subjects with new values
    for col, col_data in groups_info.items():
        for old_name, new_name in col_data.items():
            outdata.iloc[:, col][outdata.iloc[:, col] == old_name] = new_name

    # Eventually create interaction columns
    if create_interaction_col:
        interaction_col = []
        for index, row in outdata.iterrows():
            row_int_value = row[interaction_cols]
            row_int_value = [str(x) for x in row_int_value]
            row_int_value = "_:x:_".join(row_int_value)
            interaction_col.append(row_int_value)
        int_colname = "_:x:_".join(interaction_cols)
        outdata[int_colname] = interaction_col

    # Erase columns
    for col_to_erase in erase_cols:
        del outdata[col_to_erase]

    # Rename columns
    if len(rename_cols) != 0:
        new_colnames = []
        if len(rename_cols) != 1:
            raise NotImplementedError(
                "Following code works only for renaming one col. "
                "Lazy developper will update it later.")
        for old_name, new_name in rename_cols.items():
            for idx, col in enumerate(outdata.columns):
                if col == old_name:
                    new_colnames.append(new_name)
                else:
                    new_colnames.append(col)
        outdata.columns = new_colnames

    # Create outfile
    outdir = os.path.dirname(conn_file)
    outfile = os.path.basename(conn_file)
    for elt_old, elt_new in rename_file.items():
        outfile = outfile.replace(elt_old, elt_new)
    outfile = os.path.join(outdir, outfile)

    # Save output file
    outdata.to_csv(outfile, index=False)
    return outfile


def sort_connectivities_by_pvals(conn_result_file, p_val_thresh,
                                 p_val_colname, save_info_columns):
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
    save_info_columns: list
        Indicates which information to save (by column name).

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
        for col in save_info_columns:
            results_conn[col] = row[col]
        results.append(results_conn)
    return results


def prepare_permutation_on_connectivities(input_file, gpe_col,
                                          gpe_separate_cols, contrast_file,
                                          fsl_config, outdir, covariates=None,
                                          demean=True, f_contrast=None):
    """From a file extracted with extract_connectivities, create necessary
    inputs to run permutation on a connection with FSL PALM.

    Parameters
    ----------
    input_file: str
        Csv file listing for each subject its group, values for the
        different dependant variables, and eventual covariates.
        Must be comma-separated.
    gpe_col: str
        Name of column with subject's group appartenance.
    gpe_separate_cols: list of str
        Indicates how to divide one col with groups of subject in 0 and 1.
        3:IM/0 4:ctrl/1 indicates that subjects 3:IM will be in the first
        column, and 4:ctrl in the second col.
        It is very important because in order to known which contrast is
        specified.
    contrast_file: str
        contrast file.
    fsl_config: str
        path to fsl init sh file.
    outdir : str
        path to output directory.
    covariates: list of str
        Covariates names in input csv file header.
    demean: bool
        if true, demean covariates.
    f_contrast: str
        Path to F contrast txt file. If set :  test the t contrasts in a single
        F contrast.

    Returns
    -------
    design_mat_file: str
        path to design mat file
    constrat_file: str
        path to contrast file
    out_design_csv: str
        path to design mat file with header
    f_contrast_file: str
        path to f contrast file
    """

    input_data_path, input_data_ext = os.path.splitext(input_file)
    indata = pd.read_csv(input_file, dtype=str)

    # Get rid of NaN data and save NaN subjects
    nan_values_df = indata[indata.isnull().values]
    nan_subjects_csv = os.path.join(
        outdir, os.path.basename(input_file).replace(
            input_data_ext, "_nan_subjects" + input_data_ext))
    nan_values_df.to_csv(nan_subjects_csv, index=False)
    indata = indata.dropna()

    """
    Create design matrix
    """
    gpe_order = {}
    for gpe_info in gpe_separate_cols:
        gpe_info = gpe_info.split("/")
        gpe = gpe_info[0]
        gpe_col = int(gpe_info[1])
        gpe_order[gpe_col] = gpe

    # Get dependant variables values
    not_dependant_variables = [gpe_col]
    if covariates is not None:
        not_dependant_variables += covariates
    dependant_variables = [
        x for x in indata.columns if x not in not_dependant_variables]

    # Create empty dataframe for design matrix
    # > Create new cols for design matrix df
    cols = []
    for gpe_col in sorted(gpe_order.keys()):
        cols.append(gpe_order[gpe_col].replace(":", "_") + "_" + str(gpe_col))
    if covariates is not None:
        cols += covariates

    # > Create new index for design matrix df
    index_design = [x for x in range(indata.shape[0])]

    # > Create empty design matrix df
    design_mat_df = pd.DataFrame(index=index_design, columns=cols)

    # > Fill design matrix by columns
    # >> with group columns
    for gpe_col in sorted(gpe_order.keys()):
        gpe_col_name = gpe_order[gpe_col].replace(":", "_") + "_" + str(
            gpe_col)
        col_values = []
        for elt in indata[gpe_col]:
            if elt == gpe_order[gpe_col]:
                col_values.append("1")
            else:
                col_values.append("0")
        design_mat_df[gpe_col_name] = col_values

    # >> with covariates
    if covariates is not None:
        for cov in covariates:
            col_values = indata[cov]
            if demean:
                mean_val_col = np.mean([float(x) for x in col_values])
                col_values = [float(x) - mean_val_col for x in col_values]
            design_mat_df[cov] = col_values

    # > Save intermediate desin matrix output
    out_design_csv = os.path.join(outdir, "design.csv")
    design_mat_df.to_csv(out_design_csv, index=False)

    # > Save design mat file
    design_txt_file = os.path.join(outdir, "design.txt")
    design_mat_df.to_csv(design_txt_file, index=False, header=False, sep=" ")

    """
    Tranform design matrix/contrast files to valid PALM inputs
    """
    design_mat_file = design_txt_file.replace(".txt", ".mat")
    text2vest(
        indata=design_txt_file,
        outdata=design_mat_file,
        fsl_sh=fsl_config)

    # If contrast is in text file, convert in .con file
    if contrast_file.endswith(".txt"):
        constrat_file = os.path.join(
            outdir,
            os.path.basename(contrast_file).replace(".txt", ".con"))
        text2vest(
            indata=contrast_file,
            outdata=constrat_file,
            fsl_sh=fsl_config)
    else:
        constrat_file = contrast_file

    # If f contrast text file is specified, convert it to fts file
    if f_contrast is not None:
        f_contrast_file = os.path.join(
            outdir,
            os.path.basename(f_contrast).replace(".txt", ".fts"))
        text2vest(
            indata=f_contrast,
            outdata=f_contrast_file,
            fsl_sh=fsl_config)
    else:
        f_contrast_file = None

    return design_mat_file, constrat_file, out_design_csv, f_contrast_file


def create_adjacency_matrix(connections, connections_values,
                            fill_diagonal_with_zeros=False):
    """
    Create adjancency matrix for connections

    Parameters
    ----------
    connections: list of list of str
        List of list containing connected elements.
    connections_values: list of float
        List containing connection strength.
    fill_diagonal_with_zeros: bool
        For specific purpose fill matrix diagonal with zeroes (for example
        if you use the output matrix with mne plot connectome function)

    Returns
    -------
    adjacency_matrix: pandas Dataframe
        Adjacency matrix of all connections.
    node_list: list of str
        list of nodes for adjacency_matrix header
    """

    # Create a n x n matrix with n being the number of elements connected
    # Get connections names without blank space because if may be a problem
    # to use pandas loc function if index or columns names have blank names
    all_elements = []
    conn_wo_blanks = []
    for conn in connections:
        all_elements.append(conn[0])
        all_elements.append(conn[1])
        conn_wo_blanks.append(
            [conn[0].replace(" ", "_"), conn[1].replace(" ", "_")])
    all_elements = np.unique(all_elements)
    all_elements_wo_blanks = [x.replace(" ", "_") for x in all_elements]
    adjacency_matrix = np.zeros((len(all_elements), len(all_elements)))

    # Fill diagonal with 1
    if not fill_diagonal_with_zeros:
        np.fill_diagonal(adjacency_matrix, 1)

    # Transform adjacency matrix to pandas dataframe
    adjacency_matrix = pd.DataFrame(
        adjacency_matrix, index=all_elements_wo_blanks,
        columns=all_elements_wo_blanks)

    # Fill matrix with connectivities values
    for idx, conn in enumerate(conn_wo_blanks):
        adjacency_matrix.loc[
            conn[0], conn[1]] = connections_values[idx]
        adjacency_matrix.loc[
            conn[1], conn[0]] = connections_values[idx]

    return adjacency_matrix, all_elements


def parse_conn_roi_to_roi_output_conn_line(line, statistic):
    """
    Parse conn roi to roi output textfile exported from the GUI connection
    line.

    Parameters
    ----------
    line: str
        line to parse.
    statistic: str
        statistic used (T or F).

    Returns
    -------
    conn: list of str
        name of ROIs in the connection
    conn_name: str
        connection name.
    stats_info: list of str
        list of pvalues
    df: str
        degree of freedom
    residuals: str
        residuals for connection
    """
    if statistic == "T":
        line = line.split("T(")
    else:
        line = line.split("F(")

    # Get conn name
    conn = line[0]
    conn = conn.replace("Connection ", "")
    conn = conn.strip(" ")
    conn = conn.split(" -")
    conn = [x.strip(" ") for x in conn]
    conn_name = ".".join(conn)
    if len(conn) != 2:
        raise ValueError(
            "Could not split correctly connection : " + str(conn_name))
    conn_details = line[1]
    conn_details = conn_details.split(") = ")
    df_and_residuals = conn_details[0]
    stats_info = conn_details[1]
    df_and_residuals = df_and_residuals.split(",")
    df = df_and_residuals[0]
    residuals = df_and_residuals[1]
    stats_info = stats_info.split(" ")
    stats_info = [x for x in stats_info if len(x) != 0]

    return conn, conn_name, stats_info, df, residuals


def parse_conn_roi_to_roi_output_textfile(
        conn_textfile, method, statistic, conn_version, outfile):
    """
    Parse conn roi to roi output textfile exported from the GUI, to
    reorganize it into a csv file with information for each connection.

    Parameters
    ----------
    conn_textfile: str
        path to Conn ROI-to-ROI textfile.
    method: str
        name of the correction method, can be 2_SPC, 3_TFCE, 6_NBS, 4_PUS or
        0_ALL_CONNECTIONS.
    statistic : str
        name of statistic used, T or F.
    conn_version: str
        conn version
    outfile: str
        path to output file

    Returns
    -------
    outfile: str
        path to parsed conn file
    """
    if conn_version != "Conn19c":
        raise ValueError(
            "Parsing is not available for versions other than Conn19c")

    # Read data
    with open(conn_textfile, "rt") as open_file:
        lines = open_file.readlines()

    # SPC or NBS
    if method == "2_SPC" or method == "6_NBS":
        if method == "2_SPC":
            unit = "Cluster"
        else:
            unit = "Network"
        conn_roi_header = ["ROI1", "ROI2", unit]
        header = [unit + "_score_stat", unit + "_score_punc",
                  unit + "_score_pFDR",
                  unit + "_score_pFWE", unit + "_Mass_stat",
                  unit + "_Mass_punc", unit + "_Mass_pFDR",
                  unit + "_Mass_pFWE", unit + "_size_stat",
                  unit + "_size_punc", unit + "_size_pFDR",
                  unit + "_size_pFWE"]
        connections_header_cols = ["Conn_stat", "Conn_punc", "Conn_FDR"]
        results = {}
        unit_nb = 0
        nb_line = 0
        for line in lines[1:]:
            nb_line += 1
            if len(line.strip("\n").strip(" ")) == 0:
                continue
            if ("Cluster" in line) or ("Network" in line):
                unit_nb += 1
                if unit_nb in results.keys():
                    raise ValueError(
                        "Cluster/Network doublon for : {0}".format(
                            conn_textfile))
                results[unit_nb] = {
                    "Connections": {}}
                line = line.strip("\n").split("Score = ")
                line_scores = line[1].split(" ")
                line_scores = [x for x in line_scores if len(x) != 0]

                # > 4 elements are expected (score, punc, pFDR, pFWE)
                if len(line_scores) != 4:
                    msg = "Unexpected nb of elements at line " + str(nb_line)
                    msg += " file " + conn_textfile
                    raise ValueError(msg)

                results[unit_nb][unit + "_score_stat"] = line_scores[0]
                results[unit_nb][unit + "_score_punc"] = line_scores[1]
                results[unit_nb][unit + "_score_pFDR"] = line_scores[2]
                results[unit_nb][unit + "_score_pFWE"] = line_scores[3]
            elif "Connection " not in line:
                if "Mass = " in line:
                    line_split = "Cluster_Mass"
                    type_score = "Mass"
                    line = line.strip("\n").split("Mass = ")
                elif "Size = ":
                    line_split = "Cluster_size"
                    line = line.strip("\n").split("Size = ")
                    type_score = "size"
                else:
                    raise ValueError(
                        "{0} cannot determine how to process line {1}".format(
                            conn_textfile, str(nb_line)))
                line_scores = line[1].split(" ")
                line_scores = [x for x in line_scores if len(x) != 0]

                # > 4 elements are expected (score, punc, pFDR, pFWE)
                if len(line_scores) != 4:
                    msg = "Unexpected nb of elements at line " + str(nb_line)
                    msg += " file " + conn_textfile
                    raise ValueError(msg)

                results[unit_nb][
                    unit + "_" + type_score + "_stat"] = line_scores[0]
                results[unit_nb][
                    unit + "_" + type_score + "_punc"] = line_scores[1]
                results[unit_nb][
                    unit + "_" + type_score + "_pFDR"] = line_scores[2]
                results[unit_nb][
                    unit + "_" + type_score + "_pFWE"] = line_scores[3]
            else:
                (conn, conn_name, pvals_info, df,
                    residuals) = parse_conn_roi_to_roi_output_conn_line(
                        line, statistic=statistic)
                if conn_name in results[unit_nb]["Connections"].keys():
                    raise ValueError(
                        "Doublon connections : {0}".format(conn_name))
                results[unit_nb]["Connections"][conn_name] = {}
                results[unit_nb]["Connections"][conn_name]["ROI1"] = conn[0]
                results[unit_nb]["Connections"][conn_name]["ROI2"] = conn[1]
                results[unit_nb]["Connections"][conn_name][
                    "Conn_stat"] = pvals_info[0]
                results[unit_nb]["Connections"][conn_name][
                    "Conn_punc"] = pvals_info[1]
                results[unit_nb]["Connections"][conn_name][
                    "Conn_FDR"] = pvals_info[2]

    # TFCE
    elif method == "3_TFCE":
        unit = "Cluster"
        conn_roi_header = ["ROI1", "ROI2", unit]
        header = ["Cluster_Peak_TFCE_stat", "Cluster_Peak_TFCE_punc",
                  "Cluster_Peak_TFCE_pFDR", "Cluster_Peak_TFCE_pFWE"]
        connections_header_cols = ["Conn_stat", "Conn_punc", "Conn_FDR"]
        results = {}
        cluster_nb = 0
        nb_line = 0
        for line in lines[1:]:
            nb_line += 1
            if len(line.strip("\n").strip(" ")) == 0:
                continue
            if unit in line:
                cluster_nb += 1
                if cluster_nb in results.keys():
                    raise ValueError(
                        "Cluster doublon for 2_SPC : {0} cluster {1}".format(
                            conn_textfile, str(cluster_nb)))
                results[cluster_nb] = {
                    "Cluster_Peak_TFCE_stat": None,
                    "Cluster_Peak_TFCE_punc": None,
                    "Cluster_Peak_TFCE_pFDR": None,
                    "Cluster_Peak_TFCE_pFWE": None,
                    "Connections": {}}
                line = line.strip("\n").split("TFCE = ")
                line_scores = line[1].split(" ")
                line_scores = [x for x in line_scores if len(x) != 0]

                # > 4 elements are expected (score, punc, pFDR, pFWE)
                if len(line_scores) != 4:
                    msg = "Unexpected nb of elements at line " + str(nb_line)
                    msg += " file " + conn_textfile
                    raise ValueError(msg)
                results[cluster_nb]["Cluster_Peak_TFCE_stat"] = line_scores[0]
                results[cluster_nb]["Cluster_Peak_TFCE_punc"] = line_scores[1]
                results[cluster_nb]["Cluster_Peak_TFCE_pFDR"] = line_scores[2]
                results[cluster_nb]["Cluster_Peak_TFCE_pFWE"] = line_scores[3]
            elif "Connection " in line:
                line = line.strip("\n")
                (conn, conn_name, pvals_info, df,
                    residuals) = parse_conn_roi_to_roi_output_conn_line(
                        line, statistic=statistic)
                if conn_name in results[cluster_nb]["Connections"].keys():
                    raise ValueError(
                        "Doublon connections : {0}".format(conn_name))
                results[cluster_nb]["Connections"][conn_name] = {}
                results[cluster_nb]["Connections"][conn_name]["ROI1"] = conn[0]
                results[cluster_nb]["Connections"][conn_name]["ROI2"] = conn[1]
                results[cluster_nb]["Connections"][conn_name][
                    "Conn_stat"] = pvals_info[0]
                results[cluster_nb]["Connections"][conn_name][
                    "Conn_punc"] = pvals_info[1]
                results[cluster_nb]["Connections"][conn_name][
                    "Conn_FDR"] = pvals_info[2]
            else:
                raise ValueError(
                    "Cannot recognize how to parse line " + str(nb_line) + ".")

    # Parametric Univariate Statistics (PUS) or All connections
    elif method == "4_PUS" or method == "0_ALL_CONNECTIONS":
        conn_roi_header = ["ROI1", "ROI2", "Unit"]
        header = []
        connections_header_cols = ["Conn_stat", "Conn_punc", "Conn_FDR"]
        unit = "No_unit"
        results = {}
        results[unit] = {}
        results[unit]["Connections"] = {}
        for line in lines[1:]:
            if len(line.strip("\n").strip(" ")) == 0:
                continue
            line = line.strip("\n")
            (conn, conn_name, pvals_info, df,
                residuals) = parse_conn_roi_to_roi_output_conn_line(
                    line, statistic=statistic)
            if conn_name in results[unit]["Connections"].keys():
                raise ValueError(
                    "Doublon connections : {0}".format(conn_name))
            results[unit]["Connections"][conn_name] = {}
            results[unit]["Connections"][conn_name]["ROI1"] = conn[0]
            results[unit]["Connections"][conn_name]["ROI2"] = conn[1]
            results[unit]["Connections"][conn_name][
                "Conn_stat"] = pvals_info[0]
            results[unit]["Connections"][conn_name][
                "Conn_punc"] = pvals_info[1]
            results[unit]["Connections"][conn_name][
                "Conn_FDR"] = pvals_info[2]

    # Other method
    else:
        raise ValueError(
            "Parsing not implemented for method {0}.".format(method))

    # Write output
    with open(outfile, "wt") as open_file:
        open_file.write(",".join(
            conn_roi_header + header + connections_header_cols))
        open_file.write("\n")
        for unit, unit_data in results.items():
            for conn, conn_data in unit_data["Connections"].items():
                line = [conn_data["ROI1"], conn_data["ROI2"], str(unit)]
                for col in header:
                    line.append(unit_data[col])
                for col in connections_header_cols:
                    line.append(conn_data[col].strip("\n"))
                open_file.write(",".join(line))
                open_file.write("\n")
    return outfile


def conn_seed_level_fdr_correction(conn_parsed_textfile, alpha=0.05,
                                   expected_nb_conns=None):
    """
    Imitate conn seed level FDR correction.

    Parameters
    ----------
    conn_parsed_textfile: str
        path to Conn ROI-to-ROI parsed textfile containing all connections
        with p-unc values (parsed using parse_conn_roi_to_roi_output_textfile).
        Input file header must contain ROI1 and ROI2 columns for connection
        roi names and Conn_punc col for the connection uncorrected pvalue.
    outfile: str
        path to output file with FDR correction seed levels
    outfile_threshold_alpha: str
        path to output file with FDR correction seed levels connections. Only
        connections under threshold are kept.

    Returns
    -------
    outfile: str
        path to output file with FDR correction seed levels
    outfile_threshold_alpha: str
        path to output file with FDR correction seed levels connections. Only
        connections under threshold are kept.
    """

    # Load data
    conn_data = pd.read_csv(conn_parsed_textfile)

    # Get all roi names
    all_rois = np.unique(list(conn_data["ROI1"]) + list(conn_data["ROI2"]))

    # For each roi get all connections pvalues and do fdr correction
    results = {}
    for roi in all_rois:
        nb_conns = 0
        if roi in results.keys():
            raise ValueError("ROI doublon : {0}".format(roi))
        results[roi] = {}
        results[roi]["pvalues"] = {}
        results[roi]["pvalues_fdr_corrected"] = {}

        target_rois_list = []
        target_rois_list_pvalues = []

        # > Get source roi to target connections
        roi_subdata = conn_data[conn_data["ROI1"] == roi]
        if roi_subdata.shape[0] != 0:

            for idx, row in roi_subdata.iterrows():
                if row["ROI2"] in results[roi]["pvalues"].keys():
                    raise ValueError(
                        "Double connection : {0} <-> {1}".format(
                            roi, row["ROI2"]))
                results[roi]["pvalues"][row["ROI2"]] = row["Conn_punc"]
                target_rois_list.append(row["ROI2"])
                target_rois_list_pvalues.append(row["Conn_punc"])
                nb_conns += 1

        roi_subdata = conn_data[conn_data["ROI2"] == roi]
        if roi_subdata.shape[0] != 0:
            for idx, row in roi_subdata.iterrows():
                if row["ROI1"] in results[roi]["pvalues"].keys():
                    raise ValueError(
                        "Double connection : {0} <-> {1}".format(
                            roi, row["ROI1"]))
                results[roi]["pvalues"][row["ROI1"]] = row["Conn_punc"]
                target_rois_list.append(row["ROI1"])
                target_rois_list_pvalues.append(row["Conn_punc"])
                nb_conns += 1

        # > Check if number of connections from source to target is expected
        if expected_nb_conns is not None:
            if nb_conns != expected_nb_conns:
                msg = "Not expected number of targets for roi {0} :".format(
                    roi)
                msg += "{1} instead of {2}".format(
                    str(nb_conns), str(expected_nb_conns))
                raise ValueError(msg)

        # > Apply fdr and save corrected pvalues
        rejected, pvalues_corrected = fdrcorrection(
            target_rois_list_pvalues, alpha=alpha)
        for idx_tg, target in enumerate(target_rois_list):
            if target in results[roi]["pvalues_fdr_corrected"].keys():
                raise ValueError(
                    "Multiple target {0} from source {1} for corr pval".format(
                        target, roi))
            results[roi]["pvalues_fdr_corrected"][target] = pvalues_corrected[
                idx_tg]

    # Save data
    # > All connections from seed to target and fdr corrected
    with open(outfile, "wt") as open_file:
        open_file.write(
            "Seed,Target,Conn_punc,Conn_pFDR_seed_level_corrected\n")
        for roi, roi_data in results.items():
            for target, pval_unc in results[roi]["pvalues"].items():
                line = [roi, target, str(pval_unc),
                        str(results[roi]["pvalues_fdr_corrected"][target])]
                open_file.write(",".join(line))
                open_file.write("\n")

    # > All connections from seed to target that survive threshold
    with open(outfile_threshold_alpha, "wt") as open_file:
        open_file.write(
            "Seed,Target,Conn_punc,Conn_pFDR_seed_level_corrected\n")
        for roi, roi_data in results.items():
            for target, pval_unc in results[roi]["pvalues"].items():
                corrected_pval = results[roi]["pvalues_fdr_corrected"][target]
                if corrected_pval < alpha:
                    line = [roi, target, str(pval_unc), str(corrected_pval)]
                    open_file.write(",".join(line))
                    open_file.write("\n")

    return outfile, outfile_threshold_alpha
