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
import os
import json
import subprocess
from collections import OrderedDict

# Third party imports
import numpy as np
import nibabel
from nilearn import plotting


def get_image_snap(
        im_file,
        output_basename,
        modality="anat",
        mask=None,
        vmin=None,
        vmax=None,
        title=None,
        coords=None,
        cmap=None,
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
        type of image (anat or func)
    vmin: float
        min intensity value.
    vmax: float
        max intensity value.
    title: str
        title on image
    coords: int array
        display coordinates in mm.
    cmap: str
        matplotlib color map.
    display_mode : str
        nilearn display_mode.
        Can be: {'ortho', 'x', 'y', 'z', 'yx', 'xz', 'yz'}

    Returns
    -------
    outfile: str
        path to output png
    """

    # Check errors
    MODALITIES = ["anat", "func"]
    if modality not in MODALITIES:
        raise NotImplementedError(
            "Display for images that are not {0} images has yet to be "
            "added.".format(" or ".join(MODALITIES)))

    # Load data
    im = nibabel.load(im_file)
    if modality == "anat":
        max_vox_value = np.max(im.get_data())
    elif modality == "func":
        max_vox_value = np.max(im.get_data()[:, :, :, 0])
    else:
        raise NotImplementedError(
            "Display for images that are not {0} images has yet to be "
            "added.".format(" or ".join(MODALITIES)))

    # Set if needed vmin and vmax value
    # if vmin is not None:
    #     vmin = vmin * max_vox_value
    # if vmax is not None:
    #     vmax = vmax * max_vox_value

    # Create snap
    out_file = output_basename + ".png"
    if modality == "anat":

        # > Mask if needed
        if mask is not None:
            mask_im = nibabel.load(mask)
            im = nibabel.Nifti1Image(im.get_data()[mask > 0], im.affine)

        display = plotting.plot_anat(
            anat_img=im,
            cut_coords=coords,
            display_mode=display_mode,
            title=title,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap)

    elif modality == "func":

        # Display first volume
        im = nibabel.Nifti1Image(im.get_data()[:, :, :, 0], im.affine)

        # > Mask if needed
        if mask is not None:
            mask_im = nibabel.load(mask)
            im = nibabel.Nifti1Image(im.get_data()[mask > 0], im.affine)

        # Even if it is an epi sequence we still use plot_anat function
        # as plot_epi function does not serve our purpose here
        display = plotting.plot_anat(
            anat_img=im,
            cut_coords=coords,
            display_mode=display_mode,
            title=title,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap)
    else:
        raise NotImplementedError(
            "No case implemented for {0} modality.".format(modality))
    display.savefig(out_file)

    return out_file


def generate_pdf_struct_file(
        images,
        out_file,
        texts=None,
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
    texts: list
        texts for each page (multiplied by the number of image per page)
        e.g : for 3 images/page [page1, page1, page1, page2, page2, page2]
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
        page_text = []
        for i in range(nb_im_per_pages):
            if cpt == (len(images) - 1):
                with open(out_file, "wt") as open_file:
                    json.dump(page_dict, open_file, sort_keys=True,
                              check_circular=True, indent=4)
                    return
            images_page.append([images[cpt]])
            if texts is not None:
                page_text = [texts[cpt]]
            cpt += 1
        page_dict["page{0}".format(cpt_page)] = {
            "type": type_pdf,
            "style": style,
            "images": images_page,
            "texts": page_text,
            "topmargin": 0.1,
            "linecount": 120
        }
        cpt_page += 1


def fsleyes_render(
    im_file,
    out_file,
    overlay_files=[],
    overlay_types=[],
    overlay_mask_colors=[],
        mm_pos=None):
    """
    Create png snap of image using fsleyes render.

    Parameters
    ----------
    im_file: str
        path to input image file.
    out_file : str
        path to output file, without png extension.
    overlay_files: list of str
        list of overlay images
    overlay_types: list of str
        type of overlays. fsleyes render overlay types : label, linevector,
        mask, mesh, mip, rgbvector, sh, tensor, volum.
    overlay_mask_colors: list of list of int
        list of RGB color for each overlay. For mask only.
    mm_pos: list of int
        MNI pos to center snap on.

    Returns
    -------
    output_png: str
        path to output png

    /!\ Warnings : fsleyes render cannot be wrapped with pyconnectome
        FSLWrapper, FSLenv must be set before /!\
    """
    cmd = ["fsleyes", "render"]

    # Center on mni coords if needed
    if mm_pos is not None:
        cmd += ["-wl"] + [str(x) for x in mm_pos]

    cmd += ["-of", out_file]
    cmd += [im_file]

    # Add overlays
    if len(overlay_files) != 0:
        if len(overlay_files) != len(overlay_types):
            raise ValueError(
                "An overlay type must be specified for each overlay file.")

        for idx, overlay_file in enumerate(overlay_files):
            cmd += [overlay_file, "-ot", overlay_types[idx]]

            if overlay_types[idx] == "mask":
                if len(overlay_mask_colors) != 0:
                    cmd += ["--maskColour"]
                    cmd += [str(x) for x in overlay_mask_colors[idx]]

    # Run cmd
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command {0} failed : {1}".format(" ".join(
                cmd), stderr))

    # Check that file
    output_png = out_file + ".png"
    if not os.path.isfile(output_png):
        raise ValueError("Could not find {0}.".format(output_png))
    return output_png
