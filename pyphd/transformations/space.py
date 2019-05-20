# Utility functions for manipulating image space data.
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

# Third-party imports
import numpy as np
import pydicom
import nibabel
from scipy.ndimage.measurements import center_of_mass


def vox_coords_change_orientation(
        coords,
        im,
        init_orientation,
        orientation):
    """Transform voxel 3D coordinates from one orientation to another.

    E.g : Change voxel coordinates in LAS+ space to LPI+ space

    Parameters
    ----------
    coords: array
        Array of voxel coordinates.
    im: str
        path to image
    init_orientation: str
        Initial orientation. Can be LAS, RAS, LPI, RPI...
    orientation: str
        Orientation to change coordinates to. Can be LAS, RAS, LPI, RPI...

    Returns
    -------
    reoriented_coords: array
        Coordinates in new orientation.
    """

    # Change to upper case (in case)
    init_orientation = init_orientation.upper()
    orientation = orientation.upper()

    # Load image
    im = nibabel.load(im)
    x_max = im.shape[0]
    y_max = im.shape[1]
    z_max = im.shape[2]

    if ((init_orientation == "LAS" and orientation == "LPI") or
            (init_orientation == "LPI" and orientation == "LAS")):
        reoriented_coords = [coords[0], y_max - coords[1], z_max - coords[2]]
    elif ((init_orientation == "RAS" and orientation == "LPI") or
          (init_orientation == "LPI" and orientation == "RAS")):
        reoriented_coords = [x_max - coords[0], y_max - coords[1],
                             z_max - coords[2]]
    elif ((init_orientation == "PSR" and orientation == "LPI") or
          (init_orientation == "LPI" and orientation == "PSR")):
        reoriented_coords = [x_max - coords[0], y_max - coords[1],
                             z_max - coords[2]]
    else:
        raise NotImplementedError(
            "Case {0} -> {1} to be implemented yet...".format(
                init_orientation, orientation))
    return reoriented_coords


def get_im_orientation(im_file):
    """Get image orientation (LAS, RAS, etc.)

    Parameters
    ----------
    im_file: str
        Path to image file.

    Returns
    -------
    orientation: str
        Orientation.
    """
    im = nibabel.load(im_file)
    orientation = nibabel.aff2axcodes(im.affine)
    orientation = orientation[0] + orientation[1] + orientation[2]
    return orientation


def get_im_affine_info(im_file):
    """Get image information on affine matrix, sform, qform...

    Parameters
    ----------
    im_file: str
        Path to image file.

    Returns
    -------
    affine: numpy.ndarray
        Image affine matrix.
    sform_matrix: numpy.ndarray
        Image sform matrix
    qform_matrix: numpy.ndarray
        Image qform matrix
    sform_code: int
        The sform code values that specify which RAS+ space the sform affine
        refers to, with these interpretations:
        0 	unknown 	sform not defined
        1 	scanner 	RAS+ in scanner coordinates
        2 	aligned 	RAS+ aligned to some other scan
        3 	talairach 	RAS+ in Talairach atlas space
        4 	mni 	RAS+ in MNI atlas space
    """
    im = nibabel.load(im_file)
    aff = im.affine
    sform = im.header.get_sform()
    sform_code = im.header['sform_code']
    qform = im.header.get_qform()
    return aff, sform, qform, sform_code


def get_mm_to_vox_pos(coords, affine):
    """Get coordinates in voxels from coordinates in mm.

    Parameters
    ----------
    coords: list
        Point coordinates in mm.
    affine:
        image affine.

    Returns
    -------
    vox_coords: list
        Point coordinates in voxel coordinates.
        /!\ Does not return coordinates rounded /!\
    """

    # Inverse matrix
    inv_affine = np.linalg.inv(affine)

    # Multiply by mm coordinates
    coords = np.array(coords + [1])
    vox_coords = np.dot(inv_affine, coords).tolist()
    vox_coords = vox_coords[:3]

    return vox_coords


def get_vox_to_mm_pos(coords, affine):
    """Get coordinates in voxels from coordinates in mm.

    Parameters
    ----------
    coords: list
        Point coordinates in mm.
    affine:
        image affine.

    Returns
    -------
    mm_coords: list
        Point coordinates in mm coordinates.
        /!\ Does not return coordinates rounded /!\
    """
    coords = np.array(coords + [1])
    mm_coords = np.dot(affine, coords).tolist()
    mm_coords = mm_coords[:3]

    return mm_coords


def compute_im_centroid(im_file, vol=None):
    """Compute image center of mass in voxel and mm coordinates.

    Parameters
    ----------
    im_file: str
        Path to image.
    vol: int
        If 4D image, specify which volume should be used to look at the center
        of mass.

    Returns
    -------
    vox_coords: list
        Center of mass coordinates in voxel coordinates.
    mm_coords: list
        Center of mass coordinates in mm coordinates.
    """

    # Load data
    im = nibabel.load(im_file)
    im_data = im.get_data()
    if vol is not None:
        if len(im_data.shape) != 4:
            raise ValueError(
                "Please provide a 4D image if vol argument is set.")
        im_data = im_data[:, :, :, vol]

    # Compute center of mass in voxel coordinates
    vox_coords = center_of_mass(im_data)
    vox_coords = list(vox_coords)

    # Compute center of mass in mm coordinates
    aff = im.affine
    mm_coords = get_vox_to_mm_pos(vox_coords, aff)

    return vox_coords, mm_coords
