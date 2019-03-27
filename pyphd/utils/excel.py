# Utility functions for data in xsl format.
# Copyright (C) 2018  Lisa Perus
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
import re

# Third-party imports
import xlrd


def compute_excel_index(col_name):
    """Returns column number when given an excel file column name.

    E.g : Column A returns 0.

    Parameters
    ----------
    col_name: str
        Column name.

    Returns
    -------
    index: int
     Column index. Begins to 0 for first column.
    excel_columns: list
     List of excel column names.
    """

    # Generate excel columns
    alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
                "Y", "Z"]
    excel_columns = alphabet.copy()
    for idx1, col1 in enumerate(alphabet):
        for idx2, col2 in enumerate(alphabet):
            val = "{0}{1}".format(col1, col2)
            excel_columns.append(val)

    # Select column
    for idx, col in enumerate(excel_columns):
        if col == col_name:
            return idx, excel_columns
    return None, excel_columns


def excel_to_csv(xls_file, outdir, skip_row=0, header=None):
    """Convert .xls file to one or multiple csv files (one file for one sheet).


    Parameters
    ----------
    xls_file: str
        Path to xls file.
    outdir: str
        Path to output directory.
    skip_row: int
        Rows with index inferior to skip_row are skipped.
    header: array of int
        Select a line or concatenate multiple lines as the header.

    Returns
    -------
    csv_files: array
        Array of pathes to csv files.
    """
    csv_files = []
    wb = xlrd.open_workbook(xls_file)
    delete_characters = [","]
    for sheet in wb.sheets():
        print("Converting '{0}' sheet...".format(sheet.name))

        # Transform xls lines to string
        lines_csv = []
        first_col = sheet.col(0)[skip_row:]
        for idx in range(len(first_col)):
            line = sheet.row(skip_row + idx)
            line_csv = []
            for elt in line:
                elt_value = elt.value

                # Delete inconvenient characters from elt_value
                elt_value = str(elt_value)
                for char in delete_characters:
                    elt_value = elt_value.replace(char, "")
                line_csv.append(elt_value)
            line_csv = ",".join(line_csv)
            lines_csv.append(line_csv)

        # Select a header
        header_csv = []
        header_csv_merge = []
        if header is not None and len(first_col) >= len(header):

            for idx in header:
                line = sheet.row(idx)
                line_header = []
                for elt in line:
                    elt_value = elt.value

                    # Delete inconvenient characters from elt_value
                    elt_value = str(elt_value)
                    for char in delete_characters:
                        elt_value = elt_value.replace(char, "")
                        line_header.append(str(elt.value))

                header_csv.append(line_header)

            # > Check all lines in header have the same length
            header_length = len(header_csv[0])
            for line in header_csv:
                if len(line) != header_length:
                    raise ValueError(
                        "Different length for header for sheet : {0}".format(
                            sheet.name))

            # > Merge different line to create one header
            if len(header_csv) > 1:
                for idx, elt in enumerate(header_csv[0]):
                    merge = header_csv[0][idx] + " " + " ".join(
                        [x[idx] for x in header_csv[1:]])
                    header_csv_merge.append(merge)
            else:
                header_csv_merge = header_csv[0]
        basename = os.path.basename(xls_file)
        basename = re.sub("\..*$", "", basename)
        basename = basename + "_" + sheet.name
        csv_file = os.path.join(outdir, basename + ".csv")
        with open(csv_file, "wt") as open_file:
            if header is not None:
                open_file.write(",".join(header_csv_merge))
                open_file.write("\n")
            for line in lines_csv:
                open_file.write(line)
                open_file.write("\n")
        csv_files.append(csv_file)
    return csv_files
