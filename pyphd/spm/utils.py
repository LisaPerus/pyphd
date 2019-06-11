# Wrapper for spm functions not available through nipype
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
import math
import shutil
import subprocess


def spm_standalone_reorient(
    ims, sid, origin_coords, outdir, spm_sh, spm_mcr, delete_mfile=False,
        delete_mat_file=False):
    """ Wraps SPM commands to change image origin.
    ----------------------------------------------
    Works only with SPM standalone.

    Parameters
    ----------
    ims: array of str
        array of pathes to images whose origin has to be changed.
    sid: str
        subject ID.
    origin_coords: float array
        array with new origin coordinates.
    spm_sh: str
        path to SPM sh file.
    spm_mcr: str
        path to SPM MCR directory.
    delete_mfile: bool, default False
        Delete matlab script generated to reorient image.
    delete_mat_file:: bool, default False
        Delete .mat transformation file.

    TODO:
    Use nipype spm interface with SPMCommand
    e.g:
    from nipype.interfaces import spm
    matlab_cmd = '/path/to/run_spm8.sh /path/to/Compiler_Runtime/v713/ script'
    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
    """

    for im in ims:

        # Check image type
        if not im.endswith(".nii"):
            raise ValueError(
                "spm_standalone_reorient : Please provided  a .nii file")

        # Write script
        script = "path = '{0}'\n".format(im)
        script += "matlabbatch{1}.spm.util.reorient.transform.transprm = "
        script += "[{0} {1} {2} 0 0 0 1 1 1 0 0 0];\n".format(
            -origin_coords[0], -origin_coords[1], -origin_coords[2])
        script += "matlabbatch{1}.spm.util.reorient.srcfiles = "
        script += "cellstr(path);\n"
        script += "spm_jobman('run',matlabbatch);"
        script_file = os.path.join(outdir, "{0}_reorient_acpc.m".format(sid))
        with open(script_file, "wt") as open_file:
            open_file.write(script)

        # Create cmd
        cmd = [spm_sh, spm_mcr, "script", script_file]

        # Launch script
        proc = subprocess.Popen(cmd,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode == 1:
            raise ValueError(
                "Setting origin via SPM '{0}' failed : {1}".format(" ".join(
                    cmd), stderr))

        # Delete m file
        if delete_mfile:
            shutil.remove(script_file)

        # Delete mat file
        mat_file = im.replace(".nii", ".mat")
        if delete_mat_file and os.path.isfile(mat_file):
            print("Deleting file {0}...".format(mat_file))
            os.remove(mat_file)


def get_slice_order(nb_slices, order, scanner):
    """ Returns slice order for SPM slice timing correction
    -------------------------------------------------------

    Scanner modes and subsequent slice orders were found on:
    https://en.wikibooks.org/wiki/SPM/Slice_Timing

    For SIEMENS and interleaved slice order, if the number of slices is even
    slice order will be 2 4 6 1 3 5, and 1 3 5 2 4 6 if it is uneven : see
    https://practicalfmri.blogspot.com/2012/07/siemens-slice-ordering.html

    Parameters
    ----------
    nb_slices: int
        number of slices in one fmri volume.
    order: str
        slice order. Can be :
            - For SIEMENS : Ascending sequential, Ascending sequential
                            reversed, Descending sequential,
                            Descending sequential reversed,
                            Ascending interleaved,
                            Ascending interleaved reversed (Descending)
            - For PHILIPS : Default Single package, Default Two packages,
                            Default Multi-packages (>2), Ascending Single
                            package, Ascending Multi-packages, Decending
                            Single package, Decending Multi-packages,
                            Central Single package?,
                            Reverse Central Single package?,
                            Interleaved Single package?
    scanner: str
        scanner.

    Returns:
    --------
    slices_order: array of int
        slices order array.
    """
    slices_order = []

    if scanner.upper() == "SIEMENS":
        if order == "Ascending sequential":
            slices_order = [x for x in range(1, nb_slices + 1)]
        elif order == "Descending sequential":
            slices_order = [x for x in range(nb_slices, 0, -1)]

        # Even first
        elif order == "Ascending interleaved" and nb_slices % 2 == 0:
            slices_order.extend(
                [x for x in range(2, nb_slices + 1, 2)])
            slices_order.extend(
                [x for x in range(1, nb_slices + 1, 2)])

        # Odd-first
        elif order == "Ascending interleaved" and nb_slices % 2 != 0:
            slices_order.extend(
                [x for x in range(1, nb_slices + 1, 2)])
            slices_order.extend(
                [x for x in range(2, nb_slices + 1, 2)])
        elif order in [
            "Ascending sequential reversed", "Descending sequential reversed",
                "Ascending interleaved reversed"]:
            # TODO: There is a difference between SIEMENS Magnetom machines and
            # other machines for reversed sequences (see SPM wiki).
            # To be implemented later
            raise NotImplementedError(
                "Case for SIEMENS scanner and reversed sequences has not been "
                "implemented yet.")
        else:
            raise NotImplementedError(
                "Unknown sequence for SIEMENS '{0}'".format(order))
    elif scanner.upper() == "PHILIPS":
        if order == "Interleaved Single package?":
            # NOTE : For Interleaved Single Package? sequence, the slice order
            # begins to one and is spaced with an increment that is the root
            # square of the number of slices
            increment = math.sqrt(nb_slices)
            increment = math.ceil(increment)
            for idx in range(1, nb_slices + 1):
                if idx not in slices_order:
                    slices_order.extend(
                        [x for x in range(idx, nb_slices + 1, increment)])
        elif order in (
            "Default Single package", "Default Two packages",
            "Default Multi-packages (>2)", "Ascending Single package",
            "Ascending Multi-packages", " Decending Single package",
            "Decending Multi-packages", "Central Single package?",
                "Central Single package?", "Reverse Central Single package?"):
            raise NotImplementedError(
                "Lazy developper has not implemented yet sequence '{0}' for "
                "scanner PHILIPS".format(order))
    else:
        raise NotImplementedError(
            "No case implemented for scanner '{0}', only for SIEMENS and "
            "PHILIPS.".format(scanner))
    return slices_order


def st_get_ref_slice(ref_slice, slice_order):
    """ Returns slice index for slice timing correction reference slice.
    --------------------------------------------------------------------

    Returns 'First', 'Middle', or 'Last' slice index based on slice order.

    Parameters
    ----------
    ref_slice: str
        Indicates if user wants to use first, midlle (temporal) or last slice
        for slice timing correction.
    slice_order: int array
        slice order.

    Returns:
    --------
    ref_slice_idx: int
        Reference slice index.
    """

    if ref_slice == "First":
        ref_slice_idx = slice_order[0]
    elif ref_slice == "Last":
        ref_slice_idx = slice_order[-1]
    elif ref_slice == "Middle":
        if len(slice_order) % 2 == 0:
            ref_slice_idx = slice_order[int(len(slice_order) / 2) - 1]
        else:
            ref_slice_idx = slice_order[math.floor(len(slice_order) / 2)]
    else:
        raise ValueError("Ref slice must be either First, Middle or Last."
                         " Not : {0}.".format(ref_slice))
    return ref_slice_idx
