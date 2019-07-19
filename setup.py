#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name="dzutils",
    version="0.1",
    description="Dmitri Zorine's utility scripts from the Baker Lab",
    author="Dmitri Zorine",
    author_email="dzorine@gmail.com",
    packages=[
        "dzutils",
        "dzutils/pyrosetta_utils",
        "dzutils/pyrosetta_utils/phos_binding",
        "dzutils/pyrosetta_utils/geometry",
        "dzutils/pyrosetta_utils/test_tools",
        "dzutils/pyrosetta_utils/secstruct",
    ],
)
