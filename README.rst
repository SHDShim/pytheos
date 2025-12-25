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

Installation
============

This section describes how to download **pytheos** from GitHub and install it
using **conda**. The procedure is suitable for macOS and Linux.

1. Clone the GitHub repository

Clone the repository and enter the source directory::

    git clone https://github.com/SHDShim/pytheos.git
    cd pytheos

To work from a specific branch or tagged release (recommended for reproducibility)::

    git checkout master
    # or, for a tagged release
    # git checkout v0.0.2

2. Create a dedicated conda environment

Create a clean conda environment to avoid dependency conflicts::

    conda create -n pytheos python=3.11 -y
    conda activate pytheos

Python version 3.9 or newer is recommended.

3. Install core dependencies

Install numerical and scientific dependencies from ``conda-forge``::

    conda install -c conda-forge numpy scipy matplotlib sympy pandas jupyter -y

If you plan to run the example notebooks, also install Jupyter support::

    conda install -c conda-forge ipykernel nbconvert -y
    python -m ipykernel install --user --name pytheos

4. Install pytheos

From the repository root (where ``setup.py`` or ``pyproject.toml`` is located),
install pytheos::

    pip install .

Verify the installation::

    python - <<EOF
    import pytheos
    print(pytheos.__version__)
    EOF

5. Test the installation

Run a simple test to confirm that core functionality works::

    python - <<EOF
    from pytheos import bm3_p, vinet_p
    print(bm3_p(10.0, 160.0, 4.0, 4.0))
    print(vinet_p(10.0, 160.0, 4.0))
    EOF

6. Updating pytheos

To update pytheos to the latest version from GitHub::

    cd pytheos
    git pull origin master
    pip install .

7. Clean removal

To completely remove pytheos and its conda environment::

    conda deactivate
    conda remove -n pytheos --all


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
