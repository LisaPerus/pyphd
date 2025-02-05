# Utility functions for generating figures and reports.
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
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
import matplotlib.patches as mpatches
import pyphd.display.venn as venn

# Third-party imports
from nilearn import plotting
import numpy as np
from pyfreesurfer.wrapper import FSWrapper


def venn_diagram(
        groups, outdir, groups_names=None, png_basename="venn_diagram",
        title=None):
    """ Create a Venn diagram for two or three groups

    Parameters:
    -----------
    groups: array of array
        array with all the lists to use for the venn diagram
    groups_names: tuple of str
        names of the groups
    outdir: str
        path to output file
    png_basename: str
        basename for png outfile
    title: str
        figure title

    Returns:
    --------
    out_png: str
        path to output png
    """
    for group in groups:
        if len(group) == 0:
            raise ValueError("Empty group! Venn diagram cannot be drawn")

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # If no name has been defined for the groups, set a random name
    if groups_names is None:
        groups_names = []
        for idx in range(len(groups)):
            groups_names.append("Set {0}".format(idx + 1))
        groups_names = tuple(groups_names)
    groups_sets = [set(x) for x in groups]

    # Plot and save img
    if len(groups_sets) > 6:
        raise ValueError("Cannot plot Venn diagram for more than 6 groups.")
    elif len(groups_sets) == 2:
        venn2(groups_sets, groups_names)
    elif len(groups_sets) == 3:
        venn3(groups_sets, groups_names)
    elif len(groups_sets) == 4:
        labels = venn.get_labels(groups)
        fig, ax = venn.venn4(labels, names=groups_names)
    elif len(groups_sets) == 5:
        labels = venn.get_labels(groups)
        fig, ax = venn.venn5(labels, names=groups_names)
    else:
        labels = venn.get_labels(groups)
        fig, ax = venn.venn6(labels, names=groups_names)

    # > Add title
    if title is not None:
        plt.title(title)
    out_png = os.path.join(outdir, png_basename + ".png")
    plt.savefig(out_png)
    plt.close()

    return out_png


def histogram(data, outdir, png_basename):
    """ Create an histogram using pyplot.hist.

    Parameters:
    -----------
    data: array
        array of elements
    outdir: str
        path to output file
    png_basename: str
        basename for png outfile

    Returns:
    --------
    out_png: str
        path to output png
    """
    plt.hist(data)
    out_png = os.path.join(outdir, png_basename + ".png")
    plt.savefig(out_png)
    plt.close()

    return out_png


def histogram_bars(data, outdir, png_basename, xtitle=None, ytitle=None,
                   color="b"):
    """ Create an histogram using pyplot.bar.

    Parameters:
    -----------
    data: array
        array of elements
    outdir: str
        path to output file
    xtitle: str
        xtitle
    ytitle: str
        ytitle
    png_basename: str
        basename for png outfile

    Returns:
    --------
    out_png: str
        path to output png
    """

    # Count unique elements
    unique_elts = np.unique(np.array(data))
    elt_counts = []
    for elt in unique_elts:
        elt_counts.append(data.count(elt))

    # Plot elements repartition
    plt.bar(unique_elts.tolist(), height=elt_counts, color=color)
    if xtitle is not None:
        plt.xlabel(xtitle)
    if ytitle is not None:
        plt.ylabel(ytitle)
    out_png = os.path.join(outdir, png_basename + ".png")
    plt.savefig(out_png)
    plt.close()

    return out_png
