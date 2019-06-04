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
