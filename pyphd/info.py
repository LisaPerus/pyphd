#! /usr/bin/env python
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

# Module current version
version_major = 1
version_minor = 0
version_micro = 0

# Expected by setup.py: string of form "X.Y.Z"
__version__ = "{0}.{1}.{2}".format(version_major, version_minor, version_micro)

# Expected by setup.py: the status of the project
CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Environment :: X11 Applications :: Qt",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering",
               "Topic :: Utilities"]

# Project descriptions
description = """
Toolbox gathering scripts useful for my PHD work.
"""
SUMMARY = """
.. container:: summary-carousel

    pyphd is a Python module that gathers scripts and functions used during my
    PhD.

    * statistics tools
"""
long_description = """
============
pyphd
============

pyphd is a Python module that gathers scripts and functions used during my PhD.
"""

# Main setup parameters
NAME = "pyphd"
ORGANISATION = "University of Montpellier"
MAINTAINER = "Lisa Perus"
MAINTAINER_EMAIL = "lisa.perus@etu.umontpellier.fr"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = ""
DOWNLOAD_URL = ""
LICENSE = "GPLv3"
CLASSIFIERS = CLASSIFIERS
AUTHOR = "Lisa Perus"
AUTHOR_EMAIL = "lisa.perus@etu.umontpellier.fr"
PLATFORMS = "OS Independent"
ISRELEASE = True
VERSION = __version__
PROVIDES = ["pyphd"]
REQUIRES = [
    "numpy>=1.15.4",
    "nibabel>=2.3.1",
]
EXTRA_REQUIRES = {
    "obsolete": {},
    "standalone": {}
}
SCRIPTS = []
