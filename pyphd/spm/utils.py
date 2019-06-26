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

# Third-party imports
import nibabel


def spm_standalone_reorient(
    im, sid, origin_coords, outdir, spm_sh, spm_mcr, delete_mfile=False,
        delete_mat_file=False):
    """ Wraps SPM commands to change image origin.
    ----------------------------------------------

    /!\ Directly change input image origin! /!\

    Works only with SPM standalone.

    Parameters
    ----------
    im: str
       path to image whose origin has to be changed.
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

    Returns:
    --------
    im: str
        path to image whose origin has been changed.

    TODO:
    Use nipype spm interface with SPMCommand
    e.g:
    from nipype.interfaces import spm
    matlab_cmd = '/path/to/run_spm8.sh /path/to/Compiler_Runtime/v713/ script'
    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
    """

    # System imports (needed for nipype)
    import os
    import subprocess

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

    return im


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


def get_stc_parameters(im_file, slice_order, st_ref_slice, scanner):
    """ Returns all needed parameters for slice timing correction.
    --------------------------------------------------------------

    Returns number of slices, ta, tr, reference slice index and slice order.

    Parameters
    ----------
    im_file: str
        path to functional image.
    st_ref_slice: str
        Indicates if user wants to use first, midlle (temporal) or last slice
        for slice timing correction.
    slice_order: str
        slice order.
    scanner: str
        scanner.

    Returns:
    --------
    nslices: int
        Number of slices.
    ta: float
        Acquisition time.
    tr: float
        Repetition time.
    ref_slice_idx: int
        Reference slice index.
    slice_order_list:
        list of slice indexes by order of acquisition.
    """

    im = nibabel.load(im_file)
    nslices = im.shape[2]
    tr = im.header.get_zooms()[3]
    if tr == 0 or tr == 1:
        raise ValueError(
            "Suspicious TR value in NIFTI header : {0}".format(tr))
    ta = tr - (tr / nslices)
    slice_order_list = get_slice_order(
        nslices, slice_order, scanner)

    # > Compute slice of reference index
    ref_slice_idx = st_get_ref_slice(
        ref_slice=st_ref_slice,
        slice_order=slice_order_list)

    return nslices, ta, tr, ref_slice_idx, slice_order_list


def dartel_normalize_to_mni(template, flowfield, infile, fwhm,
                            voxel_size, outdir, spm_sh, spm_mcr, bb=None,
                            modulate=False):
    """ Wraps SPM Dartel function Normalize to MNI
    ----------------------------------------------

    TO READ : As for now there is a little issue with nipype wrapping of
    Normalize to MNI, hence why this wrapper was rewritten.
    TODO: Check if problem has been solved in nipype.

    Returns normalized file and normalization matrix.

    Parameters
    ----------
    template: str
        path to template generated by DARTEL.
    flowfield: str
        path to subject flowfield file (u_rc1*.nii)
    infile: str
        path to file to register.
    fwhm: array of int
        fwhm for smoothing.
        /!\ a fwhm of 0 0 0 (no smoothing) can create aliasing /!\
    voxel_size: array of int
        voxel size.
    outdir: str
        path to output directory
    spm_sh: str
        path to SPM sh file.
    spm_mcr: str
        path to SPM MCR dir.
    bb: 2x3 float array
        bounding box array. Can be set to None.
    modulate: bool
        Modulate out images - no modulation preserves concentrations. By
        default False for fMRI data.

    Returns:
    --------
    registered_input_file: str
        path to registered input file.
    template_mat_file: str
        path to template registration mat file.
    m_file: str
        path to m script used for mni normalisation.
    """

    # Setup script
    SCRIPT = """
        fprintf(1,'Executing %s at %s:\n',mfilename(),datestr(now));
        ver,
        try,

        jobs{1}.spm.tools.dartel.mni_norm.fwhm(1) = <fwhm_1>;
        jobs{1}.spm.tools.dartel.mni_norm.fwhm(2) = <fwhm_2>;
        jobs{1}.spm.tools.dartel.mni_norm.fwhm(3) = <fwhm_3>;
        jobs{1}.spm.tools.dartel.mni_norm.data.subjs.images = {...
        {...
        '<infile>';...
        };
        };
        jobs{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = {...
        '<flowfield>';...
        };


        jobs{1}.spm.tools.dartel.mni_norm.template = {...
        '<template>';...
        };
        jobs{1}.spm.tools.dartel.mni_norm.preserve = <preserve>;
        jobs{1}.spm.tools.dartel.mni_norm.vox(1) = <voxel_1>;
        jobs{1}.spm.tools.dartel.mni_norm.vox(2) = <voxel_2>;
        jobs{1}.spm.tools.dartel.mni_norm.vox(3) = <voxel_3>;

        jobs{1}.spm.tools.dartel.mni_norm.bb = <bounding_box>;

        spm_jobman('run', jobs);


        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;

        ,catch ME,
        fprintf(2,'MATLAB code threw an exception:\n');
        fprintf(2,'%s\n',ME.message);
        if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
        end;
        """

    SCRIPT = SCRIPT.replace("<fwhm_1>", str(fwhm[0]))
    SCRIPT = SCRIPT.replace("<fwhm_2>", str(fwhm[1]))
    SCRIPT = SCRIPT.replace("<fwhm_3>", str(fwhm[2]))
    SCRIPT = SCRIPT.replace("<infile>", infile)
    SCRIPT = SCRIPT.replace("<flowfield>", flowfield)
    SCRIPT = SCRIPT.replace("<template>", template)
    if modulate:
        SCRIPT = SCRIPT.replace("<preserve>", "1")
    else:
        SCRIPT = SCRIPT.replace("<preserve>", "0")
    SCRIPT = SCRIPT.replace("<voxel_1>", str(voxel_size[0]))
    SCRIPT = SCRIPT.replace("<voxel_2>", str(voxel_size[1]))
    SCRIPT = SCRIPT.replace("<voxel_3>", str(voxel_size[2]))
    if bb is not None:
        bb_str = "[{0} {1} {2}; {3} {4} {5}]".format(
            bb[0][0], bb[0][1], bb[0][2], bb[1][0], bb[1][1], bb[1][2])
    else:
        bb_str = "[NaN NaN NaN; NaN NaN NaN]"
    SCRIPT = SCRIPT.replace("<bounding_box>", bb_str)

    # Write script
    m_file = os.path.join(outdir, "pyphd_dartel_normalize_to_mni.m")
    with open(m_file, "wt") as open_file:
        open_file.write(SCRIPT)

    # Execute m-file
    cmd = [spm_sh, spm_mcr, "script", m_file]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Dartel normalization to MNI '{0}' failed : {1}".format(" ".join(
                cmd), stderr))

    # Find outputs
    # > Find template .mat file
    template_mat_file = template.replace(".nii", ".mat")

    # > Find registered input file
    if fwhm == [0, 0, 0]:
        registered_input_file = os.path.join(
            os.path.dirname(infile), "w" + os.path.basename(infile))
    else:
        registered_input_file = os.path.join(
            os.path.dirname(infile), "sw" + os.path.basename(infile))
    for output in [template_mat_file, registered_input_file]:
        if not os.path.isfile(output):
            raise ValueError("Could not find output : {0}".format(output))

    return registered_input_file, template_mat_file, m_file
