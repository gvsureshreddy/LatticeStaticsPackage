#**************************************************************************
#**************************************************************************

# File Locations ##########################################################
BIN_LOC    := ${shell pwd}/bin
INCL_LOC   := $(shell pwd)/include
LIB_LOC    := $(shell pwd)/lib

export BIN_LOC
export INCL_LOC
export LIB_LOC

.PHONY: all install clean

all: 
	$(MAKE) -e -C MyMath
	$(MAKE) -e -C LinearAlgebra
	$(MAKE) -e -C LatticeStatics
	$(MAKE) -e -C LatticeStatics/Utility -f MakeUtilityProgs

install: 
	$(MAKE) -e -C MyMath install
	$(MAKE) -e -C LinearAlgebra install
	$(MAKE) -e -C LatticeStatics install
	$(MAKE) -e -C LatticeStatics/Utility -f MakeUtilityProgs install

clean: 
	$(MAKE) -e -C MyMath clean
	$(MAKE) -e -C LinearAlgebra clean
	$(MAKE) -e -C LatticeStatics clean
	$(MAKE) -e -C LatticeStatics/Utility -f MakeUtilityProgs clean