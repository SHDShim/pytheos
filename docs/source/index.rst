.. pytheos documentation master file, created by
   sphinx-quickstart on Sat Jun  3 22:25:29 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PYTHEOS documentation
======================

.. image:: https://zenodo.org/badge/93273486.svg
   :target: https://zenodo.org/badge/latestdoi/93273486

PYTHEOS
=======

.. image:: https://zenodo.org/badge/93273486.svg
   :target: https://zenodo.org/badge/latestdoi/93273486

Overview
--------

Pytheos provides a tool set for a wide range of tasks in high pressure
science:

- calculate pressure from a number of built-in pressure scales

- convert pressure scales

- propagate uncertainties properly using the uncertainties package

- fit pressure-volume and pressure-volume-temperature data sets using the scipy and lmfit packages

- fit with a wide range of different equations and their combinations

Install
-------

Pytheos is a pure python package.

It can be installed from source or::

  pip install pytheos


Anaconda users may try::

  conda install -c shdshim pytheos


Contact
-------

Please contact Dan Shim (shdshim@gmail.com) for bug reports, comments, and
suggestions.  I am happy to include new pressure scales or other pressure
scales in pytheos as well.

Examples and Tutorials
----------------------

The pytheos package includes examples in Jupyter Notebook (under the examples
folder), which demonstrate a range of operations, calculations, and fittings
you can do with pytheos. Pytheos is designed to support data
analysis using Jupyter Notebook as well as python scripts.

Documentation is available at: https://shdshim.github.io/pytheos-docs/.

How to cite
-----------

S.-H. Shim (2017) Pytheos - a python tool set for equations of state.
Zenodo. http://doi.org/10.5281/zenodo.802392

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
