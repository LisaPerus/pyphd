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
import nibabel
import numpy as np
from pyphd.constants import CONN_INPUTS
from pyphd.fsl.utils import text2vest
import pandas as pd
from scipy.io import loadmat
from scipy.ndimage import label

try:
    from statsmodels.stats.multitest import fdrcorrection
    def conn_seed_level_fdr_correction(conn_parsed_textfile, outfile,
                                       outfile_threshold_alpha, alpha=0.05,
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
        alpha: float
            alpha level for fdr correction.
        expected_nb_conns: int
            expected number of targets for each seed.

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
                "ROI1,ROI2,Conn_punc,Conn_pFDR_seed(ROI1)_level_corrected\n")
            for roi, roi_data in results.items():
                for target, pval_unc in results[roi]["pvalues"].items():
                    line = [roi, target, str(pval_unc),
                            str(results[roi]["pvalues_fdr_corrected"][target])]
                    open_file.write(",".join(line))
                    open_file.write("\n")

        # > All connections from seed to target that survive threshold
        with open(outfile_threshold_alpha, "wt") as open_file:
            open_file.write(
                "ROI1,ROI2,Conn_punc,Conn_pFDR_seed(ROI1)_level_corrected\n")
            for roi, roi_data in results.items():
                for target, pval_unc in results[roi]["pvalues"].items():
                    corrected_pval = results[roi]["pvalues_fdr_corrected"][target]
                    if corrected_pval < alpha:
                        line = [roi, target, str(pval_unc), str(corrected_pval)]
                        open_file.write(",".join(line))
                        open_file.write("\n")

        return outfile, outfile_threshold_alpha
except ImportError:
    print("Could not import statsmodels")


def extract_connectivities(group_name, tp=None, center_name=None,
                           covariates=None, network=None,
                           conn_file_pattern=None, datafile=None, outdir=None,
                           conn_datapath=None, tp_name=None,
                           rename_conns=None):
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
    rename_conns: dict
        Dict to rename connection names.

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

    # Rename connections if needed
    if rename_conns is not None:
        delete_conns = []
        new_conn_results = {}
        for conn, conn_data in conn_results.items():
            delete_conns.append(conn)
            conn = conn.split(":")
            new_conn_name_first = rename_conns[conn[0]]
            new_conn_name_second = rename_conns[conn[1]]
            new_conn_results[
                ":".join([
                    new_conn_name_first, new_conn_name_second])] = conn_data
        delta_cols = set(conn_results.keys()) - set(delete_conns)
        for col in delta_cols:
            new_conn_results[col] = conn_results[col]
        conn_results = new_conn_results

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
    deleted_subjects: pd.DataFrame
        dataframe with lines of data deleted when extracting groups
    """

    # Read conn file
    conn_data = pd.read_csv(conn_file, header=0, dtype=str)
    outdata = conn_data.copy()

    # Select subjects in subgroups
    # Note 30/09/2020 : if there are nan values in subgroups it may delete
    # some lines
    copy_outdata_before_deletion = outdata.copy()
    for col, col_data in groups_info.items():
        outdata = outdata[outdata.iloc[:, col].isin(list(col_data.keys()))]

    # > Get eventual deleted subjects
    spe_index = list(set(list(copy_outdata_before_deletion.index)) - set(list(
        outdata.index)))
    if len(spe_index) == 0:
        deleted_subjects = {}
        for col in outdata.columns:
            deleted_subjects[col] = []
        deleted_subjects = pd.DataFrame(deleted_subjects)
    else:
        deleted_subjects = copy_outdata_before_deletion.iloc[spe_index, :]

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
    return (outfile, deleted_subjects)


def extract_clinical_data_for_mapt_metrics(
        datafile, mri_metrics_csv_file, outdir, subjects_na_var_delete=None,
        columns_na_var_delete=None, tp=None, center_name=None, covariates=None,
        cov_add_file_beginning=None, replace_cov_names={},
        delete_sid_col=False):
    """Extract and append clinical data information to MAPT MRI metrics csv
       file.

    Parameters
    ----------
    datafile : str
        path to csv datafile listing subjects clinical characteristics. Must
        have a 'Sid' column with subjects ids.
    mri_metrics_csv_file: str
        path to MRI metrics csv file. Must have a column with subject ids named
        'Sid'.
    outdir : str
        path to output directory.
    subjects_na_var_delete: list of str
        List of NA values. If these values are found for a subject in one
        column, this subject is deleted in mri_metrics_csv_file.
    columns_na_var_delete: list of str
        List of NA values in columns. If a column has only NA values it gets
        deleted.
    tp: str
        Timepoint
    center_name : str
        Center name. If set, keep only subjects from this center.
    covariates : list of str
        list of covariates. Names must be the same as those in datafile
        columns. The covariates will be added at the end of the file columns
        in the same order as they are in the list. There can be an exception
        for one covariate if the argument cov_add_file_beginning is set.
    cov_add_file_beginning: str
        name of one covariate to add to the beginning of the file.
    replace_cov_names: dict
        replace name of covariates added to file.
    delete_sid_col: bool
        If set to True delete sid column.

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # Read data
    with open(mri_metrics_csv_file, "rt") as open_file:
        mri_metrics_lines = open_file.readlines()
    with open(datafile, "rt") as open_file:
        clinical_lines = open_file.readlines()

    # Extract subjects by center if needed
    new_lines = []
    if center_name is not None:
        sid_col = None
        for idx, elt in enumerate(mri_metrics_lines[0].split(",")):
            if elt == "Sid":
                sid_col = idx
        new_lines = [mri_metrics_lines[0]]
        for line in mri_metrics_lines[1:]:
            sid = line.split(",")[sid_col].replace("sub-", "")
            sid_center = sid[0:4]
            if sid_center == center_name:
                new_lines.append(line)
    else:
        for line in mri_metrics_lines:
            new_lines.append(line)

    # Delete subjects if they have specific na values in one column
    if subjects_na_var_delete is not None:
        lines_no_na = []
        lines_no_na.append(new_lines[0])
        for line in new_lines[1:]:
            keep_subject = True
            for na_val in subjects_na_var_delete:
                line_split = line.split(",")
                if na_val in line_split:
                    keep_subject = False
            if keep_subject:
                lines_no_na.append(line)
    else:
        for line in new_lines:
            lines_no_na.append(line)

    # Delete columns that have only NA values
    lines_no_na_columns = []
    if columns_na_var_delete is not None:

        columns = lines_no_na[0].split(",")
        col_data = {}

        # > Get all elements for each column
        for line in lines_no_na[1:]:
            line = line.split(",")
            for idx, col in enumerate(columns):
                if col not in col_data.keys():
                    col_data[col] = []
                col_data[col].append(line[idx])

        # > If col has only NA values, store col idx to delete
        delete_cols_idx = []
        for idx, col in enumerate(columns):
            all_values = set(col_data[col])
            for val in columns_na_var_delete:
                if val in all_values:
                    all_values.remove(val)
            if len(all_values) == 0:
                delete_cols_idx.append(idx)

        # > Delete cols
        lines_no_na_columns = []
        for line in lines_no_na:
            line = line.split(",")
            new_line = []
            for idx, elt in enumerate(line):
                if idx not in delete_cols_idx:
                    new_line.append(elt)
            lines_no_na_columns.append(",".join(new_line))
    else:
        for line in lines_no_na:
            lines_no_na_columns.append(line)

    # Add data from clinical data
    if covariates is not None:

        # > Get subject col id in clinical data
        clinical_header = clinical_lines[0].split(",")
        clinical_sid_col = None
        for idx, col in enumerate(clinical_header):
            if col == "Sid":
                clinical_sid_col = idx

        # > Get for each covariate subjects values
        # >> Get covariates columns id in clinical data
        cov_clinical_ids = {}
        for cov in covariates:
            cov_id = None
            for idx, col in enumerate(clinical_header):
                if cov == col.strip("\n"):
                    cov_id = idx
            cov_clinical_ids[cov] = cov_id

        # >> Store for each cov subject val
        cov_subjects_clinical_data = {}
        for line in clinical_lines[1:]:
            line = line.split(",")
            if len(line) != len(clinical_header):
                raise ValueError("Bad parsing for clinical datafile.")
            sid = line[clinical_sid_col]
            for cov in covariates:
                if cov not in cov_subjects_clinical_data.keys():
                    cov_subjects_clinical_data[cov] = {}
                cov_subjects_clinical_data[cov][sid] = line[
                    cov_clinical_ids[cov]]

        # > Add cov values to new_lines
        sid_col_line_no_na_data = None
        for idx, elt in enumerate(
                lines_no_na_columns[0].strip("\n").split(",")):
            if elt == "Sid":
                sid_col_line_no_na_data = idx
        for cov in covariates:
            cov_data = cov_subjects_clinical_data[cov]
            tmp_line = lines_no_na_columns[0].strip("\n").split(",")
            if cov in replace_cov_names.keys():
                cov_name = replace_cov_names[cov]
            else:
                cov_name = cov
            if cov == cov_add_file_beginning:
                tmp_line.insert(0, cov_name)
            else:
                tmp_line.append(cov_name)
            lines_no_na_columns[0] = ",".join(tmp_line)
            for idx in range(1, len(lines_no_na_columns)):
                tmp_line = lines_no_na_columns[idx].strip("\n").split(",")
                sid = tmp_line[sid_col_line_no_na_data]
                if cov == cov_add_file_beginning:
                    tmp_line.insert(0, cov_data[sid])
                else:
                    tmp_line.append(cov_data[sid])
                lines_no_na_columns[idx] = ",".join(tmp_line)

            # >> Increment subject column id if covariate was added to the
            # beginning of the file.
            if cov == cov_add_file_beginning:
                sid_col_line_no_na_data += 1

    # Delete Sid column
    lines_final = []
    if delete_sid_col:
        sid_col_idx = None
        for idx, col in enumerate(lines_no_na_columns[0].split(",")):
            if col == "Sid":
                sid_col_idx = idx
        for line in lines_no_na_columns:
            tmp_line = line.split(",")
            new_line = []
            for idx, elt in enumerate(tmp_line):
                if idx != sid_col_idx:
                    new_line.append(elt)
            lines_final.append(",".join(new_line))
    else:
        lines_final = lines_no_na_columns

    # Write outfile
    outfile = os.path.join(outdir, os.path.basename(mri_metrics_csv_file))
    ext_file = os.path.splitext(outfile)[1]
    if center_name is not None:
        outfile = outfile.replace(
            ext_file, "_" + center_name + "_only" + ext_file)
    if subjects_na_var_delete is not None:
        outfile = outfile.replace(
            ext_file, "_nasubsdel" + ext_file)
    if columns_na_var_delete is not None:
        outfile = outfile.replace(
            ext_file, "_nacolsdel" + ext_file)
    if covariates is not None:
        outfile = outfile.replace(ext_file, "_".join(covariates) + ext_file)

    with open(outfile, "wt") as open_file:
        for line in lines_final:
            open_file.write(line)
            open_file.write("\n")

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
    if statistic == "T":
        residuals = "NA"
    else:
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


def split_roi_network_by_features(network_imfile, outfile_basename):
    """
    Split network roi images by rois.

    This function is useful if all rois in the network image are set to 1 and
    there is no way to differentiate them by intensity.

    This function is based on scipy.ndimage.label function.

    Parameters
    ----------
    network_imfile: str
        path to network image with multiple rois all set to the same intensity
        value.
    outfile_basename: str
        path to output file basename without extension.


    Returns
    -------
    outfiles: list of str
       List of all roi images
    """

    # Load network image
    network_im = nibabel.load(network_imfile)

    # Get labels images
    labelled_data, nb_components = label(network_im.get_data())

    # Check that there is a good superposition of labelled data with network
    # image if all values are set to 1
    test_labelled_data = labelled_data.copy()
    test_labelled_data[labelled_data > 0] = 1
    comparison = test_labelled_data != network_im.get_data()
    if comparison.all():
        raise ValueError(
            "Labelled data does not superpose with input network data.")

    # Create a new image for each component
    all_labels = np.unique(labelled_data)

    # > Get rid of background
    all_labels = [x for x in list(all_labels) if x != 0]

    outfiles = []
    for label_val in all_labels:
        label_data = labelled_data.copy()
        label_data[label_data != label_val] = 0
        label_data[label_data == label_val] = 1

        # > Save image
        label_im = nibabel.Nifti1Image(
            label_data, header=network_im.header, affine=network_im.affine)
        label_imfile = outfile_basename + "_" + str(int(label_val)) + ".nii.gz"
        nibabel.save(label_im, label_imfile)
        outfiles.append(label_imfile)

    return outfiles


def mapt_conn_input_filter_add_info(
    conn_input, output_file, filter_subjects, filter_colname,
        keep_cols_conn_input, add_file=None, add_file_delete_cols=["Sid"]):
    """
    Filter MAPT csv conn input datafile and add additional info.
    Function useful to select subgroup of subjects analysed in MAPT rsfmri
    study and do annex analyses (for example on cortical thickness).

    Parameters
    ----------
    conn_input: str
        Path to conn input csv
    output_file: str
        Path to output file.
    filter_subjects: dict
        Dictionnary containing subgroups of subjects with their characteristics
        in each column of Conn input.
        Example : if I want to keep subjects CDR0 who underwent omega-3 diet
        or did not underwent the diet.
        filter_colname = "CDR0_only_O3_no_O3"
        filter_subjects = {
            "CDR0_only_O3" : {"scoreCDR1" : [0],
                              "gpeMapt4c" : ["1:omega3+IM", "2:omega3"]},
            "CDR0_only_no_O3" : {"scoreCDR1 = [0],
                                 "gpeMapt4c" : ["3:IM", "4:ctrl"]}
        }

        This filtering will be saved in a new column of output_file.
        Ouput file will have a new column CDR0_only_O3_no_O3 such as

        Subject  CDR0_only_O3_no_O3  scoreCDR1   gpeMapt4c
        Sub001   CDR0_only_O3        0           2:omega3
        Sub002   CDR0_only_O3        0           1:omega3+IM
        ...
        Sub00n   CDR0_only_no_O3     0           4:ctrl

    filter_colname: str
        Name of new column of output_file representing characteristics of
        filtered subjects.
    keep_cols_conn_input: list of str
        List of columns to keep in conn_input.
    add_file: str
        Path to additional file from which to add information to output.
        This file first column must be named Sid and must contain subject
        matching with the Sid column from conn_input.
    add_file_delete_cols: list of str
        List of column from additional not be added in output_file.
    """

    # Load data
    data = pd.read_csv(conn_input, index_col=0)

    # Add center (works only for MAPT rsfmri conn inputs)
    center_col = []
    for sub in data.index:
        center_sub = sub.replace("sub-", "")[0:4]
        center_col.append(center_sub)
    data["Center"] = center_col

    # Filter MAPT Conn input data
    # > Keep subjects that fill all conditions
    subjects_kept = []
    subjects_kept_subgroups = []
    for idx_row, row in data.iterrows():
        subject_kept_for_one_subgroup = False
        for subgroup, subgroup_elts in filter_subjects.items():
            keep_subject = True
            for col, col_elts_to_keep in subgroup_elts.items():
                if row[col] not in col_elts_to_keep:
                        keep_subject = False
            if keep_subject and not subject_kept_for_one_subgroup:
                subject_kept_for_one_subgroup = True
                subjects_kept.append(row.name)
                subjects_kept_subgroups.append(subgroup)
            elif keep_subject and subject_kept_for_one_subgroup:
                raise ValueError(
                    "Subject {0} in multiple subgroups".format(row.name))

    # > Filter
    sub_data = data.loc[subjects_kept, :]
    sub_data[filter_colname] = subjects_kept_subgroups

    # Keep only some columns
    sub_data = sub_data[keep_cols_conn_input]

    # Add additional data
    if add_file is not None:
        add_data = pd.read_csv(add_file, index_col=0)
        add_data_info = {}
        for col in add_data.columns:
                add_data_info[col] = []
        for sub in sub_data.index:
            for col in add_data.columns:
                add_data_info[col].append(add_data.loc[sub, col])
        for col, col_data in add_data_info.items():
            if col not in add_file_delete_cols:
                sub_data[col] = col_data

    # Save data
    sub_data.to_csv(output_file, index=True)

    return output_file
