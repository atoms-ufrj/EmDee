EmDee: A Molecular Dynamics Laboratory
======================================

[![Build Status](https://travis-ci.org/atoms-ufrj/EmDee.svg?branch=master)](https://travis-ci.org/atoms-ufrj/EmDee)

[![codecov.io](http://codecov.io/github/atoms-ufrj/EmDee/coverage.svg?branch=master)](http://codecov.io/github/atoms-ufrj/EmDee?branch=master)

EmDee is a minimalist library that can perform, in a very efficient way, the most computationally
demanding tasks of classical molecular dynamics simulations:

1. Calculation of pairwise interaction forces, energies, and virials.
2. Building and manipulation of neighbor lists.
3. Scaling and/or translation of atom coordinates and conjugate momenta.

The purpose of EmDee is to serve as a platform to test simulation methods that involve the mechanics
of many atoms, such as different thermostats and barostats, hybrid Monte Carlo algorithms, etc.

Contents
--------

This repository contains:

* Source code (fortran) and build instructions (Makefile) for the EmDee library
* A test program in 3 versions:
  * A fortran source code + build instructions (Makefile)
  * A C source code + build instructions (Makefile)
  * A julia script
* An example input data file for the test program

Standard compilation and testing
--------------------------------

### Dependencies (considering a Ubuntu 16.04.1 LTS fresh install):

#### Basic for compiling the library

* gfortran

> **Tested with:**
> - GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609
>
> **Can usually be installed via apt:**
>
>      sudo apt-get install gfortran

#### Further dependencies for running the test

##### To run the julia script test:

* julia

> **Tested with:**
> - julia version 0.4.5
>
> **Can usually be installed via apt:**
>
>      sudo apt install julia

### Compiling the library:

* Clone the repository

        git clone https://github.com/craabreu/EmDee

  This will create a local copy of the repository on your device.

* Execute the Makefile in the root directory of the repository tree

        make

  This will build the shared library (lib/libemdee.so).

### Running the tests:

* Use building option **test**

        make test

  This will compile the c and fortran test programs and create a copy of the default julia
testscript inside the `test` directory.

* Run each test using the input file provided in the `example` directory.

        test/testfortran <N_threads> examples/data.inp
        test/testc <N_threads> examples/data.inp
        test/testjulia <N_threads> examples/data.inp

  In the commands above, N_threads (optional) is the number of parallel threads to be used for
the calculations.

## Installing the library in the system path:

* Use building option **install**

        make install

  This will copy the static library (lib/libemdee.so) to `/usr/local/lib` and the C header file
(include/emdee.h) and fortran module (include/emdee.f03) to `/usr/local/include`.
