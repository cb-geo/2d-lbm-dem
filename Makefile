## OpenACC
# CC=pgcc
# CFLAGS=-acc -openmp -Minfo -fast -Msafeptr -ta=nvidia,cc20 -tp=nehalem-64 -fpic 

## GCC
CC=gcc
CFLAGS= -g -Wall -O3 -fopenmp -lm

# Find where we are
CURDIR := $(shell pwd)

# Define destination directory
ROOTFS := ./bin

# Make sure destination directory exists before invoking any tags
$(shell [ -d "$(ROOTFS)" ] || mkdir -p $(ROOTFS))
BINDIR := $(shell pwd)

lbmdem: 
	$(CC) $(CFLAGS) src/main.c -o $(ROOTFS)/lbmdem

.PHONY: clean cleanall

clean:
	rm $ ./bin/lbmdem 

cleanall:
	rm -f $ bin/lbmdem
	rm -f $ bin/DEM*
	rm -f $ bin/stats.data
	rm -f $ bin/Fig*
	rm -f $ bin/mgp*	
	rm -f $ bin/densities*
	rm -f $ bin/LBM*
	rm -f $ bin/press*
