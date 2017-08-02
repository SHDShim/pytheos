#!/usr/bin/env python

from setuptools import setup
import pypandoc

"""
The lines below are for having README.md as the standard readme and
generate rst version automatically for the pypi.org (pip)
"""
try:
    long_description = pypandoc.convert('README.md', 'rst')
    long_description = long_description.replace("\r", "")  # Do not forget this line
except OSError:
    print("Pandoc not found. Long_description conversion failure.")
    import io
    # pandoc is not installed, fallback to using raw contents
    with io.open('README.md', encoding="utf-8") as f:
        long_description = f.read()

setup(
    setup_requires=['pbr>=1.9', 'setuptools>=17.1'],
    pbr=True,
    long_description=long_description,
)
