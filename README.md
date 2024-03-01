# SPARSKIT2_F95
A fully refactored version of SPARSKIT2, developed by Prof. Youssef Saad, from Fortran77 to Fortran95. The refactoring was done using an academic license from plusFORT.

![GitHub License](https://img.shields.io/github/license/jaadt7/SPARSKIT2_F95) ![GitHub Release](https://img.shields.io/github/v/release/jaadt7/SPARSKIT2_F95)

# Package Dependencies
This package requires only make and gfortran as a compiler.

# Using the Legacy version
If you still wish to use the fortran 77 version for whatever reason, that version was left intact in the legacy folder. Simply change directory and either do a make if you want 
a static library, or ./dotests to run the developed tests

# Using the Refactored Version
To use the fully refactored version, simply type in make in you CLI and the static library will be compiled. To run the tests, you can simply type ./run_tests.sh . You might need to chage
the executable permissions for the bash script after cloning the repository so it would run.

# Available routines
The README found in src/libary contains the list of routines and functionality coded in each source file

Authors
-------

- Jaad A. Tannous <jtannou@g.clemson.edu>
- Bradley S. Meyer <mbradle@g.clemson.edu>

Contribute
----------

- Issue Tracker: `<https://github.com/jaadt7/issues/>`_
- Source Code: `<https://github.com/jaadt7/SPARKSKIT2_F95/>`_

License
-------

The project is licensed under the GNU Lesser General Public License v2.1 (or later).

.. |license| image:: https://img.shields.io/github/license/jaadt7/SPARSKIT2_F95
    :alt: GitHub
