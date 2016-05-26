#**************************************************************************
#**************************************************************************

# File Locations ##########################################################
BIN_LOC    := ${shell pwd}/bin
INCL_LOC   := $(shell pwd)/include
LIB_LOC    := $(shell pwd)/lib

export BIN_LOC
export INCL_LOC
export LIB_LOC

#CXX = g++-4.8 -m32 # for use with KIM
#CXX = g++-4.8 -m64 # for use with FEAP
CXX = g++
OPTIMIZE   = -O
DEBUG      = -g
#PROFILE    = -pg
ANSI       = -pedantic
#SOLVE      = -DSOLVE_PLU
SOLVE      =
TERMINAL =
#TERMINAL = -DUNIX_TERMINAL
PERL = perl

export CXX
export OPTIMIZE
export DEBUG
export PROFILE
export ANSI
export SOLVE
export TERMINAL
export PERL

USE_NEO  = -DUSE_NEO
USE_ELA  = -DUSE_ELA
#USE_KIM  = -DUSE_KIM
#USE_FEAP = -DUSE_FEAP
#USE_QC   = -DUSE_QC

export USE_NEO
export USE_ELA
export USE_KIM
export USE_FEAP
export USE_QC

.PHONY: all install clean

all:
	$(MAKE) -e -C MyMath
	$(MAKE) -e -C LinearAlgebra
	$(MAKE) -e -C LatticeStatics

install:
	$(MAKE) -e -C MyMath install
	$(MAKE) -e -C LinearAlgebra install
	$(MAKE) -e -C LatticeStatics install

uninstall:
	$(MAKE) -e -C MyMath uninstall
	$(MAKE) -e -C LinearAlgebra uninstall
	$(MAKE) -e -C LatticeStatics uninstall

clean:
	$(MAKE) -e -C MyMath clean
	$(MAKE) -e -C LinearAlgebra clean
	$(MAKE) -e -C LatticeStatics clean
