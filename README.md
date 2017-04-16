2D Lattice Boltzmann -Discrete Element Method
=============================================
> Cambridge Berkeley - Geomechanics

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/lbm-dem/master/license.md)
[![CircleCI](https://circleci.com/gh/cb-geo/2d-lbm-dem.svg?style=svg)](https://circleci.com/gh/cb-geo/2d-lbm-dem)
[![](https://img.shields.io/github/issues-raw/cb-geo/2d-lbm-dem.svg)](https://github.com/cb-geo/2d-lbm-dem/issues)

## Prerequisites
OpenACC v1.0 or higher


## Compile and Run

0. Run `mkdir build && cd build && cmake ..`

1. Run `make clean && make -jN` (where N is the number of cores)

3. Run lbmdem `./bin/lbmdem /<path-to-inputfile>/inputfile.dat`


