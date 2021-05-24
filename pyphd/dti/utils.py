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
        dti_datafile, cols_to_keep, outdir, additional_outfile_name=None):
    """Extract from file with dti input and covariates covariate of interest
       and path to dti file.

    Parameters
    ----------
    dti_datafile: str
        Path to file containing path to dti scalar file and covariates.
        Comma-separated. Need to be read with pandas because some covariate
        have comma in one col.
        Each scalar col is <scalar>_<timepoint>.
    cols_to_keep: list of str
        List of columns to keep in dti_datafile.
    outdir: str
        path to output directory
    additional_outfile_name: str
        add additional name to outfile.

    Returns
    -------
    outfile: str
        Path to file subjects with all their connection
    """

    # Read data
    dti_data = pd.read_csv(dti_datafile, header=0, dtype=str)

    # Get col to keep
    dti_data = dti_data[cols_to_keep]

    # Save data
    outfile = os.path.join(outdir, "dti_data.txt")
    if additional_outfile_name is not None:
        outfile = outfile.replace(
            ".txt", "_" + additional_outfile_name + ".txt")
    outfile = outfile.replace(
        ".txt", "_{0}.txt".format("_".join(cols_to_keep)))
    dti_data.to_csv(outfile, index=False)

    return outfile
