# Utility functions about system files.
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


def get_file_size(path_file):
    """Returns the size of a file in bytes.

    Parameters
    ----------
    path_file: str
        Path to file.

    Returns
    -------
    size: float
        File size in bytes.
    """
    size = os.stat(path_file).st_size
    return size


def check_dir_copy(indir, copy_dir):
    """Check if one directory is a copy of another, if the
       two directories have the same contents.

    Parameters
    ----------
    indir: str
        Path to input directory content that was copied.
    copy_dir: str
        Path to copy directory.

    Returns
    -------
    exact_copy: bool
        True if copy_dir content is the same as indir content.
    diff: list of str
        content in indir that is not in copy_dir
    """

    diff = []

    # List of files and directories in indir
    to_check_files = []
    to_check_dirs = []

    # Walk through the directories
    print("Cd to {0} ...".format(indir))
    os.chdir(indir)
    for root, dirnames, filenames in os.walk("."):
        to_check_files += [os.path.join(root, x) for x in filenames]
        to_check_dirs += [os.path.join(root, x) for x in dirnames]

    # Go to copy dir and check that file and dirs have been copied correctly
    print("Cd to {0} ...".format(copy_dir))
    os.chdir(copy_dir)
    exact_copy = True
    for directory in to_check_dirs:
        if not os.path.isdir(directory):
            exact_copy = False
            diff.append(directory)
    for fid in to_check_files:
        if not os.path.isfile(fid):
            exact_copy = False
            diff.append(fid)
    return exact_copy, diff
