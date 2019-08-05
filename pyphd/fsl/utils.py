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


def palm(indata, design_file, contrast_file, output_basename, nb_permutations):
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
    output_basename: str
        Path to output basename.
    nb_permutations: int
        Number of permutation (default : 1000).

    Returns
    -------
    stat_val: array of str
        path to csv files listing test values.
    pval_unc: array of str
        path to csv files listing uncorrected p-values for each contrast.
    pval_fwe: array of str
        path to csv files listing FWE-corrected p-values for each contrast.
    """
    cmd = ["palm", "-i", indata, "-d", design_file, "-t", contrast_file, "-o",
           output_basename, "-n", str(nb_permutations)]
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode == 1:
        raise ValueError(
            "Command '{0}' failed : {1} + {2}".format(" ".join(cmd), stderr,
            stdout))

    stat_val = glob.glob(os.path.join(output_basename + "*dat_tstat_c*.csv"))
    pval_unc = glob.glob(os.path.join(output_basename + "*uncp*.csv"))
    p_fwe = glob.glob(os.path.join(output_basename + "*fwep*.csv"))
    return stat_val, pval_unc, p_fwe
