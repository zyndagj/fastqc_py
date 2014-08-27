#!/usr/bin/python

"""
Setup script for fastqc-py
"""

from distutils.core import setup

setup(name = "fastqc_py",
	version = "0.1",
	author = "Greg Zynda",
	author_email="gzynda@tacc.utexas.edu",
	license="GNU",
	descritpion="Quality control for fastq files",
	requires=[
		'numpy',
		'matplotlib'],
	packages = ["fastqc_py"])
