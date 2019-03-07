#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name="dzutils",
    version="0.1",
    description="Dmitri Zorine's utility scripts from the Baker Lab",
    author="Dmitri Zorine",
    author_email="dzorine@gmail.com",
    packages=["pyrosetta_utils"],
    py_modules=["stringcheese", "func_utils", "pdb_file_utils"],
)
