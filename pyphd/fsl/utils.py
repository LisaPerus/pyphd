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

# System import
import os
import glob
import subprocess
import pandas as pd
import numpy as np

# Pyconnectome imports
from pyconnectome.wrapper import FSLWrapper


def palm(indata, design_file, contrast_file, f_contrast, output_basename,
         nb_permutations, twotail=False, saveparametric=False,
         alternate_palm_bin=None, singularity_cmd=None):
    """ Wraps FSL PALM command.
    ---------------------------

    Parameters
    ----------
    indata: str
        Path to csv or image used as input
    design_file: str
        Path to .mat design file.
    contrast_file: str
        Path to contrast .con file.
    f_contrast: str
        F contrast fts file.
    output_basename: str
        Path to output basename.
    nb_permutations: int
        Number of permutation (default : 1000).
    twotail: bool
        Run two-tailed tests for all the t-contrasts instead of
        one-tailed.
    saveparametric: bool
        Option to save parametric pval.
    alternate_palm_bin: str
        Path to alternate palm bin file : useful for cluster.
    singularity_cmd : list of str
        singularity command to use with palm if palm is contained in a
        singularity image.

    Returns
    -------
    stat_val: array of str
        path to csv files listing test values.
    pval_unc: array of str
        path to csv files listing uncorrected p-values for each contrast.
    pval_fwe: array of str
        path to csv files listing FWE-corrected p-values for each contrast.
    """
    if alternate_palm_bin is None:
        cmd = ["palm"]
    else:
        cmd = [alternate_palm_bin]
    if singularity_cmd is not None:
        cmd = singularity_cmd + cmd
    cmd += ["-i", indata, "-d", design_file, "-t", contrast_file, "-o",
            output_basename, "-n", str(nb_permutations)]
    if f_contrast is not None:
        cmd += ["-f", f_contrast]
    if twotail:
        cmd.append("-twotail")
    if saveparametric:
        cmd.append("-saveparametric")
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
                "Command '{0}' failed : {1} + {2}".format(
                    " ".join(cmd), stderr, stdout))
    stat_val = glob.glob(os.path.join(output_basename + "*dat_tstat_c*.csv"))
    pval_unc = glob.glob(os.path.join(output_basename + "*uncp*.csv"))
    p_fwe = glob.glob(os.path.join(output_basename + "*fwep*.csv"))
    return stat_val, pval_unc, p_fwe


def text2vest(indata, outdata, fsl_sh):
    """ Wraps FSL Text2Vest command.
    ---------------------------

    Parameters
    ----------
    indata: str
        Path to text file input data. Can be a text file for design matrix
        or contrast.
    outdata: str
        Path to output data (.mat or .con file).
    fsl_sh: str
        Path to fsl init sh file.
    """
    cmd = ["Text2Vest", indata, outdata]
    fslprocess = FSLWrapper(cmd, shfile=fsl_sh)
    fslprocess()


def create_group_design_matrix(
        datafile, gpe_col, outdir, covariates=None, merge_cols=None,
        demean=True):
    """
    Create a group design matrix file from a csv file containing a column with
    group of subjects and other information.

    E.g:
    input_data:
    Alzheimer  age sex
    AD           70  0
    Healthy      80  1
    AD           75  1

    output:
    Alzheimer_1_for_AD   Alzheimer_1_for_Healthy age sex
    1                    0                       70   0
    0                    1                       80   1
    1                    0                       75   1

    This function is similar to one part of the script
    pyphd/pyphd/scripts/pyphd_prepare_permutation_connectivities.py

    Parameters
    ----------
    datafile: str
        path to datafile with group column and additional columns with
        covariates.
    gpe_col: str
        name of column in datafile with group information.
    outdir: str
        path to output directory
    covariates: list of str
        list of columns with covariates
    merge_cols: list of str
        contains names of columns that need to be merged. If set, all columns
        will be merged and the merged column will become the new gpe_col.
        WARNING :
    demean: bool
        if true, demean covariates.
    nb_subjects[name]: dict
        nb of subjects/variables per gpe_col level

    Returns
    -------
    design_file: str
        path to design file.
    design_file_with_header: str
        path to design file with additional header.
    nb_gpe_cols:
        number of columns concerning group information.
    """

    # Read data
    data = pd.read_csv(datafile)

    # Merge columns if necessary
    if merge_cols is not None:
        new_col = "_".join(merge_cols)
        new_col_value = []
        for idx, row in data.iterrows():
            sub_val = row[merge_cols]
            sub_val = "_".join(sub_val)
            new_col_value.append(sub_val)
        data[new_col] = new_col_value
        gpe_col = new_col
    group_values = set(list(data[gpe_col]))

    # If group values come from a merged column
    if merge_cols is not None:
        group_values = []
        if len(merge_cols) != 2:
            raise NotImplementedError(
                "Did not implement merge for more than 2 cols")
        first_col_vals = set(data[merge_cols[0]])
        second_col_vals = set(data[merge_cols[1]])
        for col in first_col_vals:
            for scol in second_col_vals:
                val = col + "_" + scol
                group_values.append(val)
        for col in merge_cols:
            del data[col]
    else:
        group_values = sorted(list(group_values))

    nb_gpe_cols = len(group_values)

    # Create columns for groups
    groups_dict = {}
    nb_subjects = {}
    outdata = pd.DataFrame()
    for val in group_values:
        name = gpe_col + "_1_for_{0}".format(str(val))
        nb_subjects[name] = 0
        col = []
        for elt in list(data[gpe_col]):
            if elt == val:
                nb_subjects[name] += 1
                col.append("1")
            else:
                col.append("0")
        outdata[name] = col

    # Add covariates
    if covariates is not None:
        for cov in covariates:
            col_values = list(data[cov])
            if demean:
                mean_val_col = np.mean([float(x) for x in col_values])
                col_values = [float(x) - mean_val_col for x in col_values]
            outdata[cov] = col_values

    # Write output
    design_file_with_header = os.path.join(outdir, "design.csv")
    design_file = os.path.join(outdir, "design.txt")
    outdata.to_csv(design_file_with_header, index=False)
    outdata.to_csv(design_file, index=False, header=False, sep=" ")

    return design_file, design_file_with_header, nb_gpe_cols, nb_subjects


def randomise(
    input_file, output_file, design_file, contrast_file, fsl_sh, mask=None,
        f_contrast=None, tfce=True, tfce_2d_opt=False, nb_perms=None,
        demean=False, fstat_only=False, raw_statistic_im=False):
    """ Wraps FSL randomise command.
    ---------------------------

    For the moment implements only TFCE correction.

    Parameters
    ----------
    input_file: str
        path to 4D input image.
    output_file: str
        path to output file
    design_file: str
        path to design file.
    contrast_file: str
        path to contrast file
    fsl_sh: str
        path to fsl init sh file.
    mask: str
        path to mask file
    f_contrast:
        path to f contrast file if f contrast needed
    tfce: bool
        apply tfce correction
    tfce_2d_opt: bool
        apply tfce correction with 2D optimisation (for TBSS data)
        see command doc.
    nb_perms: int
        nb of permutations. If none, default of 5000 permutations is set.
    demean: bool
        demean data temporally before model fitting
        ( demean model as well if required )
    fstat_only: bool
        calculate f-statistics only
    raw_statistic_im: bool
        outputs raw statistic image. For tfce will output
        for example <output>_randomise_tfce_tstat1.nii.gz. For tfce, this
        raw statistic image can be generated using fslmaths with:
        fslmaths tstat1.nii.gz -tfce 2 1 26 tfce_tstat1.nii.gz

    Returns
    -------
    stat_files: list of str
        list of stat files
    tfce_files: list of str
        list of tfce fwe corrected pvalues files
    cmd: str
        randomise command that was run.
    """
    cmd = ["randomise", "-i", input_file, "-o", output_file,
           "-d", design_file, "-t", contrast_file]
    if demean:
        cmd += ["-D"]
    if mask is not None:
        cmd += ["-m", mask]
    if f_contrast is not None:
        cmd += ["-f", f_contrast]
    if tfce and not tfce_2d_opt:
        cmd += ["-T"]
    if tfce_2d_opt:
        cmd += ["--T2"]
    if nb_perms is not None:
        cmd += ["-n", str(nb_perms)]
    if fstat_only:
        cmd += ["--fonly"]
    if raw_statistic_im:
        cmd += ["-R"]
    print("Running randomise...")
    fslprocess = FSLWrapper(cmd, shfile=fsl_sh)
    fslprocess()
    cmd = " ".join(cmd)

    # Get stat files
    stat_files = glob.glob(output_file + "_tstat*.nii.gz")

    # Get tfce files
    tfce_files = glob.glob(output_file + "_tfce_corrp_tstat*.nii.gz")

    return stat_files, tfce_files, cmd


def find_contrast(nb_gpe_cols, nb_covariates, model, one_sided=True,
                  multiple_levels_nb=None):
    """ Find pre-made contrast files for specific group comparison models.
    Parameters
    ----------
    nb_gpe_cols: int
        Number of groups for group comparison.
    nb_covariates: int
        Number of covariates to adjust for.
    model: str
        Name of the statistic models used. Can be glm, ttest, anova, ancova.
    one_sided: bool
        If True outputs both opposite contrasts for t contrast.
    gpe_col: str
        name of column in datafile with group information.
    outdir: str
        path to output directory
    covariates: list of str
        list of columns with covariates
    demean: bool
        if true, demean covariates.
    multiple_levels_nb: list of int
        if model is complex (e.g : interaction between variable with multiple
        level) each variable number of level is specified.

    Returns
    -------
    contras_file: str
        path to contrast text file.
    fts_file: str
        path to f contrast text file.
    """
    contrast_file = None
    fts_file = None
    if model == "ttest":
        if one_sided:
            contrast_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                "ressources", "contrasts", "ttest_one_sided.txt")
        else:
            contrast_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                "ressources", "contrasts", "ttest_two_sided.txt")
    elif model == "glm":
        if nb_covariates == 4:
            if one_sided:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "glm_2gpes_4covs_one_sided.txt")
            else:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "glm_2gpes_4covs_two_sided.txt")
        elif nb_covariates == 5:
            if one_sided:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts",
                    "glm_2gpes_5covs_one_sided.txt")
            else:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "glm_2gpes_5covs_two_sided.txt")
        else:
            raise NotImplementedError(
                "Does not provide pre-made contrast file for glm without 4 "
                "or 5covs")
    elif model == "anova":
        if nb_gpe_cols == 4:
            contrast_file = os.path.join(
                os.path.dirname(os.path.dirname(
                    os.path.realpath(__file__))),
                "ressources", "contrasts", "anova_4gpes.txt")
            fts_file = os.path.join(
                os.path.dirname(os.path.dirname(
                    os.path.realpath(__file__))),
                "ressources", "contrasts", "anova_4gpes_fcontrast.txt")
        else:
            raise NotImplementedError(
                "Does not provide pre-made contrast file for anova with more "
                "or less than 4 groups")
    elif model == "ancova":
        if nb_gpe_cols == 4:
            fts_file = os.path.join(
                os.path.dirname(os.path.dirname(
                    os.path.realpath(__file__))),
                "ressources", "contrasts", "anova_4gpes_fcontrast.txt")
            if nb_covariates == 4:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "ancova_4gpes_4covs.txt")
            elif nb_covariates == 5:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "ancova_4gpes_5covs.txt")
            else:
                raise NotImplementedError(
                    "Does not provide pre-made contrast file for ancova with "
                    "4 or 5 covariates")
        else:
            raise NotImplementedError(
                "Does not provide pre-made contrast file for ancova with more "
                "or less than 4 groups")
    elif model == "two_way_anova":
        if multiple_levels_nb[0] == 2 and multiple_levels_nb[1] == 4:
            if nb_covariates == 4:
                contrast_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts", "two_way_anova_2x4_4covs.txt")
                fts_file = os.path.join(
                    os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))),
                    "ressources", "contrasts",
                    "two_way_anova_2x4_4covs_fcontrast.txt")
            else:
                raise ValueError(
                    "Does not provide pre-made contrast file for 2x4 anova "
                    "with more or less than 4 covariates")
        else:
            raise ValueError(
                "Does not provide pre-made contrast file for two_way anova "
                "other than 2x4 design")
    else:
        raise NotImplementedError(
            "Do not provide contrast file other than for glm, ttest, anova, "
            "ancova and two_way_anova.")
    return contrast_file, fts_file
