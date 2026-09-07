# SPECFEM2D

### This fork is extended from the original source (https://github.com/geodynamics/specfem2d.git) to couple specfem2D with another numerical framework to model dynamic earthquake rupture. An acceleration time series can be given explicitly as the source of a specfem2D simulation.

**See <http://htmlpreview.github.io/?https://github.com/kura-okubo/specfem2d/blob/master/EXAMPLES/Validation/note/Note_of_Specfem2DCoupling.html>**

Ported from v7.0.0 to upstream **v8.1.0**; see [PORTING_v7_to_v8.md](PORTING_v7_to_v8.md)
for what changed and how it was validated.

[![DOI](https://zenodo.org/badge/14293189.svg)](https://zenodo.org/badge/latestdoi/14293189)<br>

SPECFEM2D allows users to perform 2D and 2.5D (i.e., axisymmetric) simulations
of acoustic, elastic, viscoelastic, and poroelastic seismic wave propagation.
The package can also be used for full waveform imaging (FWI) or adjoint tomography.


SPECFEM2D was founded by Dimitri Komatitsch and Jeroen Tromp, and is now being developed by a large, collaborative, and inclusive community. A complete list of authors can be found at
https://specfem2d.readthedocs.io/en/latest/authors/

> 2026.09.07: Updated the version originally developed in 2019 based on v7.0.0 to incorporate updates from v8.1.0.


## Installation

Instructions on how to install and use `SPECFEM2D` are
available in the

- PDF manual located in directory: [doc/USER_MANUAL](doc/USER_MANUAL).

- HTML manual (latest version): [specfem2d.readthedocs.io](http://specfem2d.readthedocs.io/)

For a quick test, run the default example with these commands:
```
./configure FC=gfortran
make all
./bin/xmeshfem2D
./bin/xspecfem2D
```
and check the output files in `./OUTPUT_FILES/`

>__NOTE__: Do not modify the 'configure' script directly. Please modify the 
    'configure.ac' file instead, and generate a new 'configure' script with 
    the command: `autoreconf -i`


## Development

[![Actions Status](https://github.com/kura-okubo/specfem2d/actions/workflows/CI.yml/badge.svg)](https://github.com/kura-okubo/specfem2d/actions)
[![codecov](https://codecov.io/gh/SPECFEM/specfem2d/branch/devel/graph/badge.svg)](https://codecov.io/gh/SPECFEM/specfem2d)
[![Documentation Status](https://readthedocs.org/projects/specfem2d/badge/?version=latest)](https://specfem2d.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

* Actions tests: [github actions kura-okubo/specfem2d](https://github.com/kura-okubo/specfem2d/actions)


Development is hosted on GitHub in the
[SPECFEM/specfem2d repository](https://github.com/SPECFEM/specfem2d).


SPECFEM2D is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive and
participatory community so we are happy that you are interested in
participating! We have collected a set of guidelines and advice on how to get
involved in the community and keep them in the specfem3d github wiki:
[specfem3d wiki](https://github.com/SPECFEM/specfem3d/wiki)


## Computational Infrastructure for Geodynamics (CIG)

SPECFEM2D is part of the software that is hosted by the Computational Infrastructure for Geodynamics (CIG). It is available on the CIG website [here (SPECFEM2D)](https://geodynamics.org/resources/specfem2d).
