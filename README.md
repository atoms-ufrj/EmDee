EmDee: A Molecular Dynamics Laboratory
======================================

EmDee is a minimalist library that can perform, in a very efficient way, the most computationally
demanding tasks of classical molecular dynamics simulations:

1. Calculation of pairwise interaction forces, energies, and virials.
2. Building and manipulation of neighbor lists.
3. Scaling and/or translation of atom coordinates and conjugate momenta.

The purpose of EmDee is to serve as a platform to test simulation methods that involve the mechanics
of many atoms, such as different thermostats and barostats, hybrid Monte Carlo algorithms, etc.

Contents
------------

This repository contains:

* Source code (fortran) and build instructions (Makefile) for the EmDee library
* A test program in 3 versions:
  * A fortran source code + build instructions (Makefile)
  * A c source code + build instructions (Makefile)
  * A julia script
* An example input data file for the test program

Standard installation and testing
------------

### Dependencies (considering a Ubuntu 16.04.1 LTS fresh install):

#### Basic for compiling the library

* gfortran

> **Tested with:**
>- GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609
>
> **Can usually be installed via apt:**
>
>      sudo apt-get install gfortran

#### Further dependencies for running the test

* to run the julia script test:

        sudo apt-get install julia

### Compiling the library

* Clone the repository

        git init
        git remote add upstream https://github.com/craabreu/EmDee
         git pull upstream master
  This will create a local copy of the repository on your device.

* Execute the Makefile in the root directory of the repository tree

        make
  This will build the shared and static libraries (libemdee.so and libemdee.a).

### Running the tests:

* Run the label **install** from the Makefile with administrator permissions:

        sudo make install
  This will copy the library, c header, and fortran 2003 header-like files) to your system folders (`/usr/local/lib/` and `/usr/local/include/`) and run `ldconfig`.

* Run the label **test**

        make test
  This will compile the c and fortran test programs and  create a copy of the default julia test script inside the `test` directory.

* Run each test using the input file provided in the `example` directory.

        ./test/testfortran examples/data.inp
        ./test/testc examples/data.inp
        julia ./test/testjulia examples/data.inp
