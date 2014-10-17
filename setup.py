#!/usr/bin/python

"""
Setup script for fastqc-py
"""

from distutils.core import setup, Extension

mod1 = Extension('cIO', sources=["fastqc_py/cIO.c"])

setup(name = "fastqc_py",
	version = "0.2",
	author = "Greg Zynda",
	author_email="gzynda@tacc.utexas.edu",
	license="GNU",
	descritpion="Quality control for fastq files",
	requires=[
		'numpy',
		'matplotlib'],
	packages = ["fastqc_py"],
	ext_modules=[mod1],
	package_data={'fastqc_py':['adapter_sequences.fa']})
