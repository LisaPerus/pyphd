# Useful functions for DTI data analyses
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

# System imports
import os

# Third party imports
import pandas as pd


def extract_dti_group_analysis_data(
    dti_datafile, clinical_datafile, dti_data_sid_col, clinical_data_sid_col,
        group_name, outdir, covariates=None, additional_outfile_name=None):
    """Outputs file with merged info about dti and clinical data.

    Parameters
    ----------
    dti_datafile: str
        Path to file containing path to dti scalar file and subjects ids.
        Comma-separated.
    clinical_datafile: str
        Path to file containing for subjects clinical data.
        Comma-separated.
    dti_data_sid_col: str
        Name of subjects IDs column in dti_datafile.
    clinical_data_sid_col: str
        Column with names of subjects IDs in clinical_datafile.
    group_name: str
        Add information from column group_name in clinical_datafile
        for each subject.
    outdir: str
        Path to output directory.
    covariates: list of str
        list of covariates to add to outfile.
    additional_outfile_name: str
        add additional name to outfile.

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # Read data
    dti_data = pd.read_csv(dti_datafile)
    clinical_data = pd.read_csv(clinical_datafile, sep=",")
    clinical_data.index = clinical_data[clinical_data_sid_col]

    # Get group values from clinical_data
    subjects = dti_data[dti_data_sid_col]
    group_values = []
    for sub in subjects:
        group_values.append(clinical_data.loc[sub, group_name])

    # Get covariates if necessary
    covariates_data = {}
    if covariates is not None:
        for cov in covariates:
            covariates_data[cov] = []
            for sub in subjects:
                covariates_data[cov].append(clinical_data.loc[sub, cov])

    # Merge data
    output_data = dti_data.copy()
    output_data[group_name] = group_values
    if len(covariates_data) != 0:
        for cov in covariates:
            output_data[cov] = covariates_data[cov]

    # Save data
    outfile = os.path.join(outdir, "dti_data.txt")
    if additional_outfile_name is not None:
        outfile = outfile.replace(
            ".txt", "_" + additional_outfile_name + ".txt")
    outfile = outfile.replace(".txt", "_{0}.txt".format(group_name))
    if covariates is not None:
        outfile = outfile.replace(
            ".txt", "_{0}.txt".format("_".join(covariates)))
    output_data.to_csv(outfile, index=False)

    return outfile
