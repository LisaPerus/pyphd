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

# Pyconnectome imports
from pyconnectome.wrapper import FSLWrapper


def palm(indata, design_file, contrast_file, f_contrast, output_basename,
         nb_permutations, twotail=False, alternate_palm_bin=None):
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
    alternate_palm_bin: str
        Path to alternate palm bin file : useful for cluster.

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
    cmd += ["-i", indata, "-d", design_file, "-t", contrast_file, "-o",
           output_basename, "-n", str(nb_permutations)]
    if f_contrast is not None:
        cmd += ["-f", f_contrast]
    if twotail:
        cmd.append("-twotail")
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
