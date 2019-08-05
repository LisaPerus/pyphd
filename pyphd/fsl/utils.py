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

# Pyconnectome imports
from pyconnectome.wrapper import FSLWrapper


def palm(indata, design_file, contrast_file, output_basename, fsl_sh,
         nb_permutations):
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
    fsl_sh: str
        Path to FSL sh init file.
    nb_permutations: int
        Number of permutation (default : 1000).

    Returns:
    --------


    """
    cmd = ["palm", "-i", indata, "-d", design_file, "-t", contrast_file, "-o",
           output_basename, "-n", nb_permutations]
    fslprocess = FSLWrapper(cmd, shfile=fsl_sh)
    fslprocess()
