PYTHEOS
=======

|DOI|

Overview
--------

``Pytheos`` provides a tool set for a wide range of tasks in high
pressure science:

-  calculate pressure from a number of built-in pressure scales

-  convert pressure scales

-  propagate uncertainties properly using the ``uncertainties`` package

-  fit pressure-volume and pressure-volume-temperature data sets using
   the ``scipy`` and ``lmfit`` packages

-  fit with a wide range of different equations and their combinations

Install
-------

``Pytheos`` is a pure python package.

It can be installed from source or,

::

    $ pip install pytheos

For ``anaconda`` users,

::

    $ conda install -c shdshim pytheos

Instruction for mineral physicists
----------------------------------

I provide environment files for mineral physicists. The environment file
includes useful packages for mineral physicists. So those who are
interested only ``pytheos`` may not want to use the environment files.

If you decide to use the environment, you may install the environment
by:

::

    $ conda env create -f py35ds-osx.yml

For window users, use: ``py35ds-win64.yml``.

The command above will create a new environment: ``py35ds``.

If you need more detailed one, please find my instruction in `my gist
site <https://gist.github.com/SHDShim/4f5987e4e1693b10dfa025baa9ab6f9d>`__.

Contact
-------

Please contact Dan Shim (shdshim@gmail.com) for bug reports, comments,
and suggestions. I am happy to include new pressure scales or other
pressure scales in ``pytheos`` as well.

Examples and Tutorials
----------------------

The ``pytheos`` package includes examples in Jupyter Notebook (under the
examples folder), which demonstrate a range of operations, calculations,
and fittings you can do with ``pytheos``. ``Pytheos`` is designed to
support data analysis using Jupyter Notebook as well as python scripts.

Documentation is available at: https://shdshim.github.io/pytheos-docs/.

How to cite
-----------

S.-H. Shim (2017) Pytheos - a python tool set for equations of state.
Zenodo. http://doi.org/10.5281/zenodo.802392

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.802392.svg
   :target: https://doi.org/10.5281/zenodo.802392
