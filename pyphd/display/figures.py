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
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

# Third-party imports
from nilearn import plotting
import numpy as np


def venn_diagram(
        groups, outdir, groups_names=None, png_basename="venn_diagram"):
    """ Create a Venn diagram for two groups

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

    Returns:
    --------
    out_png: str
        path to output png
    """

    # Checks
    if len(groups) > 5:
        raise ValueError(
            "Venn diagram with more than 5 groups is visually horrible."
            "Please find another representation for your data.")

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
    venn2(groups_sets, groups_names)
    out_png = os.path.join(outdir, png_basename + ".png")
    plt.savefig(out_png)

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


def plot_connectome(adjacency_matrix, coords, outname, patches=None):
    """
    Plot connectome using nilearn functions

    Parameters:
    -----------
    adjacency_matrix: numpy ndarray
        square matrix of connections values
    coords: list of list of int
        List of coordinates for each element of the matrix.
        Must be in the same order as elements in the square matrix.
    outname: str
        output name.
    patches: list of list of str
        list of list for patches in the form : [["label", "color"]]
        Must be in the same order as elements in the square matrix.
    """
    coords = np.array(coords)

    # Plot
    if patches is None:
        display = plotting.plot_connectome(adjacency_matrix, coords)
    else:
        handle_patches = []
        node_colors = []
        for node_info in patches:
            label, color = node_info
            node_patch = mpatches.Patch(color=color, label=label)
            handle_patches.append(node_patch)
            node_colors.append(color)
        display = plotting.plot_connectome(
            adjacency_matrix, coords, node_color=node_colors)
        plt.legend(
            handles=handle_patches, loc="upper left", bbox_to_anchor=(-0.1, 1))

    # Save plot
    display.savefig(outname)
    display.close()
