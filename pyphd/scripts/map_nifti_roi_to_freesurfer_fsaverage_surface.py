# System import
import os
import argparse
import textwrap
import glob
import re
import json
import shutil
import subprocess
import time
import progressbar
from datetime import datetime
from argparse import RawTextHelpFormatter
from datetime import datetime
from pprint import pprint

# Third party imports
from pyfreesurfer.wrapper import FSWrapper

# Script documentation
DOC = """
Map ROI in nifti format to freesurfer fsaverage surface
-------------------------------------------------------

python3.5 $SCRIPT_DIR/GIT_REPOS/pyphd/pyphd/scripts/map_nifti_roi_to_freesurfer_fsaverage_surface.py \
    -a ~/PHD/DATA/Functional_ROIs/Cambridge_atlas/template_cambridge_basc_multiscale_nii_asym/template_cambridge_basc_multiscale_asym_scale036_regions/resampled_to_fsl_6.0.0_MNI152_T1_2mm_brain/MNI152_T1_2mm_brain.nii \
    -r ~/PHD/DATA/Functional_ROIs/Cambridge_atlas/template_cambridge_basc_multiscale_nii_asym/template_cambridge_basc_multiscale_asym_scale036_regions/resampled_to_fsl_6.0.0_MNI152_T1_2mm_brain/resliced_to_fsl6.0.0_MNI152_T1_2mm_braintemplate_cambridge_basc_multiscale_asym_scale036_cluster*.nii \
    -s $NFS_CATI/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/fsaverage/surf \
    -f $HOME/setup_freesurfer_i2bm.sh \
    -o $NFS_CATI/MAPT/MAPT_T1MRI/database_freesurfer/freesurfer_v6_longitudinal/03/roi_to_fsaverage_surface \
    -V 2


Notes:

1) If a ROI span over two hemisphere (for example Cambridge R36 region 15)
do I get the same result if I split the ROI into left and right and then
transform it into a freesurfer surface or if I directly transform the whole
ROI to a surface (on both hemisphere)

Answer : yes! You get exactly the same result.

2) If the ROI is not one component but is formed by multiple components
 (it is more of a network) do I get the same results if I project each
 component or the whole ROI

Answer : yes! To get the same results between the two different inputs, do not
forget to weight each component thickness by its size, in case you split the
roi into multiple components.
"""


def is_file(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError("File does not exist: %s" % filepath)
    return filepath


def is_directory(dirarg):
    """ Type for argparse - checks that directory exists.
    """
    if not os.path.isdir(dirarg):
        raise argparse.ArgumentError(
            "The directory '{0}' does not exist!".format(dirarg))
    return dirarg


# Parse input arguments
def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python map_nifti_roi_to_freesurfer_fsaverage_surface.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--anat-im", type=is_file, required=True,
        help="Anatomical image in the same space as the ROI image.")
    required.add_argument(
        "-r", "--roi-im", type=is_file, required=True, nargs="+",
        help="ROI image in the same space as the anatomical image.")
    required.add_argument(
        "-s", "--freesurfer-fsaverage-surf-dir", type=is_directory,
        required=True, help="Path to freesurfer fsaverage surf dir.")
    required.add_argument(
        "-f", "--fs-sh", type=is_file, help="Path to Freesurfer sh file.")
    required.add_argument(
        "-o", "--outdir",
        type=is_directory, required=True, metavar="<path>",
        help="Path to the output directory.")

    # Optional argument
    parser.add_argument(
        "-V", "--verbose",
        type=int, choices=[0, 1, 2], default=2,
        help="Increase the verbosity level: 0 silent, [1, 2] verbose.")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")
    return kwargs, verbose


"""
Parse the command line.
"""
inputs, verbose = get_cmd_line_args()
runtime = {
    "timestamp": datetime.now().isoformat(),
    "tool": "map_nifti_roi_to_freesurfer_fsaverage_surface.py"
}
outputs = {}
if verbose > 0:
    pprint("[info] Map nifti roi to freesurfer fsaverage surface...")
    pprint("[info] Runtime:")
    pprint(runtime)
    pprint("[info] Inputs:")
    pprint(inputs)

# Chdir to fsaverage surf dir
# cd $SUBJECTS_DIR/fsaverage/surf
os.chdir(inputs["freesurfer_fsaverage_surf_dir"])

# Register anatomical image to fsaverage
reg_file = os.path.join(
    inputs["outdir"], "{0}_to_fsaverage.dat".format(os.path.basename(
        inputs["anat_im"])))
cmd = ["fslregister", "--s", "fsaverage", "--mov", inputs["anat_im"], "--reg",
       reg_file]
fscmd = FSWrapper(cmd, inputs["fs_sh"])
fscmd()
if not os.path.isfile(reg_file):
    raise ValueError("Cannot find T1 to fsaverage reg file {0}...".format(
        reg_file))
outputs["Transformation {0} to fsaverage".format(
        os.path.basename(inputs["anat_im"]))] = reg_file

# Map the ROI-mask(s) to the fsaverage surface, to create an fsaverage-ROI
# surface overlay
for roi_im in inputs["roi_im"]:
    _, roi_im_ext = os.path.splitext(roi_im)
    roi_name = os.path.basename(roi_im).replace(roi_im_ext, "")
    for side in ["lh", "rh"]:
        outfile = os.path.join(
            inputs["outdir"], "{0}.fsaverage.{1}.mgh".format(side, roi_name))
        cmd = ["mri_vol2surf",
               "--mov", roi_im,
               "--reg", reg_file,
               "--projdist-max", "0", "1", "0.1",
               "--interp", "nearest",
               "--hemi", side,
               "--out", outfile]
        fscmd = FSWrapper(cmd, inputs["fs_sh"])
        fscmd()
        if not os.path.isfile(outfile):
            raise ValueError(
                "Cannot find roi fsaverage surface file {0}...".format(
                    outfile))
        outputs["{0} fsaverage surface file".format(roi_name)] = outfile


"""
Write outputs
"""
logdir = os.path.join(inputs["outdir"], "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
for name, final_struct in [("inputs", inputs), ("outputs", outputs),
                           ("runtime", runtime)]:
    log_file = os.path.join(
        logdir,
        "map_nifti_roi_to_freesurfer_fsaverage_surface_{0}.json".format(name))
    with open(log_file, "wt") as open_file:
        json.dump(final_struct, open_file, sort_keys=True, check_circular=True,
                  indent=4)

if verbose > 0:
    pprint("[Outputs]:")
    pprint(outputs)
