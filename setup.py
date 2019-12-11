#!/usr/bin/env python
# -*- coding: utf-8 -*-

# System imports
import os
from setuptools import setup, find_packages

# Get module info
release_info = {}
infopath = os.path.join(os.path.dirname(__file__), "pyphd", "info.py")
with open(infopath) as open_file:
    exec(open_file.read(), release_info)
setup(
    name=release_info["NAME"],
    version=release_info["VERSION"],
    packages=find_packages(),
    author=release_info["AUTHOR"],
    author_email=release_info["AUTHOR_EMAIL"],
    description=release_info["DESCRIPTION"],
    long_description=release_info["LONG_DESCRIPTION"],
    url=release_info["URL"],
    classifiers=release_info["CLASSIFIERS"],
    license=release_info["LICENSE"],
    include_package_data=True
)


