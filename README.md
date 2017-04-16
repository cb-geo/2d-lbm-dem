2D Lattice Boltzmann -Discrete Element Method
==============================================
> Cambridge-Berkeley Computational Geomechanics

## Author(s)
* Krishna Kumar, Department of Engineering, University of Cambridge
* Jean-Yves Delenne, INRA, France
* Kenichi Soga, University of California, Berkeley

## Prerequisites
OpenACC v1.0 or higher


## Compile and Run

0. Run `mkdir build && cd build && cmake ..`

1. Run `make clean && make -jN` (where N is the number of cores)

3. Run lbmdem `./bin/lbmdem /<path-to-inputfile>/inputfile.dat`


