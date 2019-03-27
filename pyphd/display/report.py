# Utility functions for generating reports.
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
import json
from collections import OrderedDict

# Third party imports
import nibabel
from nilearn import plotting


def get_image_snap(
        im_file,
        output_basename,
        modality="anat",
        title=None,
        coords=None,
        display_mode="ortho"):
    """
    Generate image snap using nilearn functions.

    Parameters
    ----------
    im_file: str
        Path to image file.
    output_basename: array
        output(s) basename.
    modality: str
        type of image (anat or epi)
    title: str
        title on image
    coords: int array
        display coordinates in mm.
    display_mode : str
        nilearn display_mode.
        Can be: {'ortho', 'x', 'y', 'z', 'yx', 'xz', 'yz'}

    Returns
    -------
    outputs: array
        array of path to outputs png
    """

    # Load data
    im = nibabel.load(im_file)
    outputs = []

    if modality != "anat":
        raise NotImplementedError(
            "Display for images that are not anatomical images has yet to be "
            "added.")
    out_file = output_basename + ".png"
    display = plotting.plot_anat(
        anat_img=im_file,
        cut_coords=coords,
        display_mode=display_mode,
        title=title)
    display.savefig(out_file)
    outputs.append(out_file)

    return outputs


def generate_pdf_struct_file(
        images,
        out_file,
        nb_im_per_pages=3,
        style="OneCol",
        type_pdf="triplanar"):
    """
    Generate pdf struct useful for pyconnectome generate_pdf function.

    Parameters
    ----------
    images: array
        array of images
    out_file : str
        path to json outfile
    nb_im_per_pages: int
        number of images per page.
    style: str
        style of each pdf page
    type_pdf: str
        type of pdf
    """
    data = OrderedDict()
    data["cover"] = {"type": "cover"}
    cpt = 0
    cpt_page = 1
    page_dict = OrderedDict()
    while 1:
        images_page = []
        for i in range(nb_im_per_pages):
            if cpt == (len(images) - 1):
                with open(out_file, "wt") as open_file:
                    json.dump(page_dict, open_file, sort_keys=True,
                              check_circular=True, indent=4)
                    return
            images_page.append([images[cpt]])
            cpt += 1
        page_dict["page{0}".format(cpt_page)] = {
            "type": type_pdf,
            "style": style,
            "images": images_page,
            "texts": [],
            "topmargin": 0.1,
            "linecount": 120
        }
        cpt_page += 1
