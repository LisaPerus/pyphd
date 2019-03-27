# Utility functions getting info from dicom header.
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
import re
import subprocess
from shutil import which

# Third-party imports
import numpy as np
import pydicom


def get_vox_mm_pos_in_slice(dcm_file, vox_pos=[0, 0]):
    """Returns from a dicom file a voxel coordinates in the frame's image plane
    in units of mm.

    Formula is derived from this :
    https://dicom.innolitics.com/ciods/cr-image/general-series/00102210

    Parameters
    ----------
    dcm_file: str
        Path to dicom image.
    vox_pos: int array
        Voxel coordinates in voxel units.

    Returns
    -------
    pos_mm: float array
        Voxel coordinates in mm units.
    """

    dcm = pydicom.read_file(dcm_file)
    image_position_patient = dcm["0x0020", "0x0032"].value  # S_xyz
    image_orientation_patient = dcm["0x0020", "0x0037"].value
    X_xyz = image_orientation_patient[:3]
    Y_xyz = image_orientation_patient[3:]
    pixel_spacing = dcm["0x0028", "0x0030"].value
    pixel_spacing_i = pixel_spacing[0]
    pixel_spacing_j = pixel_spacing[1]

    # Compute vox to mm matrix
    matrix = np.zeros((4, 4))

    matrix[0, 0] = X_xyz[0] * pixel_spacing_i
    matrix[1, 0] = X_xyz[1] * pixel_spacing_i
    matrix[2, 0] = X_xyz[2] * pixel_spacing_i

    matrix[0, 1] = Y_xyz[0] * pixel_spacing_j
    matrix[1, 1] = Y_xyz[1] * pixel_spacing_j
    matrix[2, 1] = Y_xyz[2] * pixel_spacing_j

    matrix[0, 3] = image_position_patient[0]
    matrix[1, 3] = image_position_patient[1]
    matrix[2, 3] = image_position_patient[2]
    matrix[3, 3] = 1

    pos_mm = np.dot(
        matrix, np.array([vox_pos[0], vox_pos[1], 0, 1]))
    return pos_mm


def get_csa_header(dcm_file, scanner):
    """Extract CSA header from dicom file (SIEMENS only) using gdcmdump tool.

    Parameters
    ----------
    dcm_file: str
        Path to dicom image.

    Returns
    -------
    csa_header: str
        CSA header
    """

    # Test scanner
    if "SIEMENS" not in scanner.upper():
        raise ValueError(
            "CSA header extraction only available for SIEMENS dicom.")

    # Test if gdcmdump is available
    cmd_exists = which('gdcmdump')
    if cmd_exists is None:
        print("Command gdcmdump unavailable, cannot read CSA header.")
        return None

    cmd = ['gdcmdump', '--csa', dcm_file]
    csa_header = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
    csa_header = csa_header.decode("utf-8")
    return csa_header


def get_slice_timing_from_csa_header(dcm_file, scanner):
    """Extract slices timings from CSA header of a dicom SIEMENS-only file.

    Parameters
    ----------
    dcm_file: str
        Path to dicom image.

    Returns
    -------
    slice_timings : array of float
        slices timings
    slice_order : array of int
        order of slice acquisitions
    """
    if "SIEMENS" not in scanner.upper():
        raise ValueError(
            "Slice timings extraction from CSA header only available for "
            "SIEMENS dicom.")

    csa_header = get_csa_header(dcm_file, scanner)
    csa_header = csa_header.split("\n")

    # Get line with attribute MosaicRefAcqTimes
    line_acq_times = None
    for idx, line in enumerate(csa_header):
        if "MosaicRefAcqTimes" in line:
            line_acq_times = line

    # Get data
    slice_timings = re.findall(r"Data .*$", line_acq_times)
    if len(slice_timings) == 0:
        print("No slice timings could be extracted from csa header. Exiting.")
        return None, None
    elif len(slice_timings) > 1:
        print("Multiple patterns for slice timings in csa header. Choosing the"
              " first one.")
        slice_timings = slice_timings[0]
    else:
        slice_timings = slice_timings[0]
    slice_timings = slice_timings.replace("Data", "")
    slice_timings = slice_timings.split("\\")
    slice_timings = [x.replace(" ", "") for x in slice_timings]
    slice_timings = [x.replace("\"", "") for x in slice_timings]
    slice_timings = [x.replace("\'", "") for x in slice_timings]
    slice_timings = [float(x) for x in slice_timings]

    # Compute slice order
    slice_order = []
    sorted_slice_timings = sorted(slice_timings)
    for st in sorted_slice_timings:
        for idx, st2 in enumerate(slice_timings):
            if st == st2:
                slice_order.append(idx + 1)
    return slice_timings, slice_order
