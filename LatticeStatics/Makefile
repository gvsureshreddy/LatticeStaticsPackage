#**************************************************************************
#**************************************************************************

VPATH = bin lib depend Lattices Modes SolutionMethods Potentials Utility

UNAME      = $(shell uname)

EXE        = LatticeStatics
SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
BUILD_DATE = \"`date +%c`\"

# Compiling options, etc. #################################################
#
gCC        = g++

ifeq ($(UNAME),Darwin)
   # gCC compiler
   CC         = $(gCC)
   OPTIMIZE   = -O
   ANSI       = -pedantic -Wall
endif

ifeq ($(UNAME),Linux)
   # PGI compiler
   #CC         = pgCC
   #OPTIMIZE   = -fast -Minline=levels:10 --no_exceptions
   #ANSI       = -A -Wall
   CC	= $(gCC)
   OPTIMIZE = -O
   ANSI	= -pedantic -Wall
endif

# general options
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
#BOUNDS     = -DCHECK_BOUNDS
TERMINAL   = -DUNIX_TERMINAL
DATE_DEF   = -DBUILD_DATE="$(BUILD_DATE)"
SOLVE      = -DSOLVE_SVD
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE) $(TERMINAL)


# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LIBS       = -lLinearAlgebra -lMyMath -lm
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(DATE_DEF)
SUBDIRS    = Lattices Modes SolutionMethods Potentials Utility

# Defines the source files needed #########################################
include $(addsuffix /Makefile, $(SUBDIRS))
SOURCES  = ./LatticeStatics.cpp\
	   $(addprefix Lattices/,$(LATTICES))\
	   $(addprefix Modes/,$(MODES))\
	   $(addprefix SolutionMethods/,$(SOLUTIONMETHODS))\
	   $(addprefix Potentials/,$(POTENTIALS))\
	   $(addprefix Utility/,$(UTILITY))


# Defines the object files to be created using the SYSTEM name ############
STRIPEDSOURCES = $(notdir $(SOURCES))
OBJECTS  = $(STRIPEDSOURCES:$(SRC_EXT)=$(OBJ_EXT))

# Defines the Program Target (linking command) ############################
all:     $(EXE) Utilities

Utilities:
	$(MAKE) -e -C ./Utility  -f MakeUtilityProgs


$(EXE):  $(OBJECTS)
	 $(CC) $(CCFLAGS) $(DATE_DEF) -o lib/builddate.o -c Utility/builddate.cpp
	 $(PURE) $(CC) $(LFLAGS) `perl -MExtUtils::Embed -e ldopts` lib/builddate.o \
		 -L$(LIB_LOC) -o ./bin/$(EXE) \
		 $(addprefix lib/,$(OBJECTS)) $(LIBS)
	 /bin/rm -f lib/builddate.o

# include the dependicies
include $(addprefix depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
PerlObject.o: Makefile
	$(CC) $(CCFLAGS) `perl -MExtUtils::Embed -e ccopts` \
        -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        -o ./lib/$(notdir $@) \
        -c Utility/PerlObject.cpp

$(filter-out PerlObject.o, $(OBJECTS)): Makefile
	$(CC) $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        -o ./lib/$(notdir $@) \
        -c $(filter %/$(*),$(basename $(SOURCES)))$(SRC_EXT)

KnownPotentials.o: Makefile
	echo $(patsubst %.h,%,$(PTNLS))

.PHONY: install uninstall
# Define install
install: $(BIN_LOC)/$(EXE) Utilities
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs install

$(BIN_LOC)/$(EXE):	$(EXE)
	/bin/cp bin/$(EXE) $(BIN_LOC)
	/bin/cp gpl_scripts/AddRank1Convex $(BIN_LOC)
	/bin/cp gpl_scripts/MultiLattice2gpl $(BIN_LOC)
	/bin/cp gpl_scripts/BlochWaveRecalc $(BIN_LOC)
	/bin/cp gpl_scripts/CoarsenOutput $(BIN_LOC)
	/bin/cp gpl_scripts/CoarsenGPL $(BIN_LOC)

# Define uninstall
uninstall:
	/bin/rm -f $(BIN_LOC)/$(EXE)
	/bin/rm -f $(BIN_LOC)/AddRank1Convex
	/bin/rm -f $(BIN_LOC)/MultiLattice2gpl
	/bin/rm -f $(BIN_LOC)/BlochWaveRecalc
	/bin/rm -f $(BIN_LOC)/CoarsenOutput
	/bin/rm -f $(BIN_LOC)/CoarsenGPL
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs uninstall


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix lib/,$(OBJECTS)) \
		  $(addsuffix /*~,$(SUBDIRS))\
		  $(addprefix lib/,$(OBJECTS:$(OBJ_EXT)=.ti))\
                  lib/$(EXE).* bin/$(EXE) core *~ \
		  depend/*.d
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs clean

# Define rule for dependicy files
depend/$(notdir %).d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCC) -MM -I. $(addprefix -I ,$(SUBDIRS)) \
	      -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"

