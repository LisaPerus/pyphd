# Constants file.
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

# Constant for fmri data preprocessing script adapted from Clara Manesco
# pipeline.
RSFMRI_PREPROC_CLARA_MANESCO = {
    "NewSegment":
        {
            "_comment": ["Save only native, not dartel map.",
                         "channel_info contains : bias regularization, FWHM,",
                         "which maps to save (Field, Corrected)"],
            "tissues": [[["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 1], 1, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 2], 1, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 3], 2, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 4], 3, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 5], 4, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 6], 2, [True, False], [False, False]]
                        ],
            "write_deformation_fields" : [False, False],
            "channel_info" : [0.001, 60.0, [False, False]],
            "affine_regularization" : "mni",
            "warping_regularization" : [0, 0.001, 0.5, 0.05, 0.2],
            "sampling_distance" : 3.0
        },
    "NewSegmentDartel":
        {
            "tissues": [[["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 1], 1, [True, True], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 2], 1, [True, True], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 3], 2, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 4], 3, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 5], 4, [True, False], [False, False]],
                         [["/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii", 6], 2, [True, False], [False, False]]
                        ],
            "write_deformation_fields" : [True, True],
            "channel_info" : [0.001, 60.0, [False, False]],
            "affine_regularization" : "mni",
            "warping_regularization" : [0, 0.001, 0.5, 0.05, 0.2],
            "sampling_distance" : 3.0
        },
    "Normalize12":
        {
            "affine_regularization_type": "mni",
            "bias_regularization": 0.0001,
            "tpm": "/home/lp254663/opt/spm12/spm12/spm12_mcr/spm12/tpm/TPM.nii",
            "warping_regularization" : [0, 0.001, 0.5, 0.05, 0.2],
            "bias_fwhm": 60,
            "smoothness": 0,
            "sampling_distance": 3.0,
            "write_bounding_box": [[-78, -112, -70],
                                   [78, 76, 85]],
            "write_voxel_sizes": [2, 2, 2],
            "write_interp": 4
        },
    "Smooth":
        {
            "fwhm": [6.0, 6.0, 6.0],
            "data_type": 0,
            "implicit_masking": False,
            "out_prefix": "s"
        },
    "SliceTiming":
        {
            "out_prefix" : "a"
        },
    "Realign":
        {
            "quality" : 0.9,
            "separation" : 4.0,
            "fwhm" : 5.0,
            "register_to_mean" : True,
            "interp" : 2,
            "wrap" : [0, 0, 0],
            "weight_img" : None,
            "write_which" : [2, 1],
            "write_interp" : 4,
            "write_wrap" : [0, 0, 0],
            "write_mask" :  True,
            "out_prefix" : "r",
            "jobtype" : "estwrite"
        },
    "Coregister":
        {
            "cost_function" : "nmi",
            "separation" : [4.0, 2.0],
            "tolerance" : [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01,
                           0.01, 0.001, 0.001, 0.001],
            "fwhm" : [7.0, 7.0],
            "write_interp" : 4,
            "write_mask" :  False,
            "write_wrap" : [0, 0, 0]
        },
    "DartelNormalize2MNI":
    {
        "bounding_box" : [[-78, -112, -70],
                          [78, 76, 85]],
        "fwhm" : {"func" : [6.0, 6.0, 6.0],
                  "anat" : [1.0, 1.0, 1.0]},
        "modulate" : False,
        "voxel_size" : (2.0, 2.0, 2.0)
    }
}

SCIMAGOJR_RANKING_DATA = {
    "data" : os.path.join(os.path.dirname(os.path.realpath(__file__)), "ressources", "scimagojr_2018.csv"),
    "date" : "2019-12-11",
    "link" : "https://www.scimagojr.com/journalrank.php?out=xls",
    "version" : "2018"
}
