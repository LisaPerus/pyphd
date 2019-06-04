# Utils functions to compute statistics on nifti image.
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

# Third party imports
import numpy as np
import nibabel


def mean_brain_val(im_file, mask_file=None):
    """ Get mean value over a nifti file.

    Parameters
    ----------
    im_file: str
        Path to nifti file.
    mask_file: str
        Path to mask image over which to compute the mean.

    Returns
    -------
    mean_im_val: float
        mean image value.
    """

    # Load data
    im = nibabel.load(im_file)

    if mask_file is None:
        mean_im_val = np.mean(im.get_data())
    else:
        print("Averaging on mask : {0} ...".format(mask_file))
        im_masked_values = []
        mask_im = nibabel.load(mask_file)

        # Assume voxels > 0 are to be kept in the mask
        vox_pos_mask = np.where(mask_im.get_data() > 0)

        for idx in range(len(vox_pos_mask[0])):
            val = im.get_data()[
                vox_pos_mask[0][idx], vox_pos_mask[1][idx],
                vox_pos_mask[2][idx]]
            im_masked_values.append(val)
        im_masked_values = np.array(im_masked_values)
        mean_im_val = np.mean(im_masked_values)

    return mean_im_val


def substract_mean_brain_val(im_file, outdir, mask_file=None):
    """ Substract to the image its mean value.

    Parameters
    ----------
    im_file: str
        Path to nifti file.
    outdir: str
        Path to output directory
    mask_file: str
        Path to mask image over which to compute and substract the mean.

    Returns
    -------
    im_minus_mean_file: str
        Path to image with substracted mean value.
    """
    # Load data
    im = nibabel.load(im_file)
    new_im_data = im.get_data().copy()

    # Compute mean val
    mean_im_val = mean_brain_val(im_file, mask_file)

    if mask_file is None:
        new_im_data -= mean_im_val
    else:
        mask_im = nibabel.load(mask_file)
        new_im_data[mask_im.get_data() > 0] -= mean_im_val

    # Save data
    new_im = nibabel.Nifti1Image(new_im_data, im.affine)
    im_minus_mean_file = os.path.join(outdir, "m" + os.path.basename(im_file))
    nibabel.save(new_im, im_minus_mean_file)

    return im_minus_mean_file


def divide_brain_by_mean_val(im_file, outdir, mask_file=None):
    """ Divide the image by its mean value.

    Parameters
    ----------
    im_file: str
        Path to nifti file.
    outdir: str
        Path to output directory
    mask_file: str
        Path to mask image over which to compute and substract the mean.

    Returns
    -------
    im_minus_mean_file: str
        Path to image with substracted mean value.
    """
    # Load data
    im = nibabel.load(im_file)
    new_im_data = im.get_data().copy()

    # Compute mean val
    mean_im_val = mean_brain_val(im_file, mask_file)

    if mask_file is None:
        new_im_data /= mean_im_val
    else:
        mask_im = nibabel.load(mask_file)
        new_im_data[mask_im.get_data() > 0] /= mean_im_val

    # Save data
    new_im = nibabel.Nifti1Image(new_im_data, im.affine)
    im_minus_mean_file = os.path.join(outdir, "d" + os.path.basename(im_file))
    nibabel.save(new_im, im_minus_mean_file)

    return im_minus_mean_file
