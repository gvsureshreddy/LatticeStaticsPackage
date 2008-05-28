# -*- mode: makefile -*-
#**************************************************************************
#**************************************************************************

VPATH = bin lib depend

UNAME      = $(shell uname)

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
   #ANSI       = -A -Minform=inform
   CC = $(gCC)
   OPTIMIZE = -O
   ANSI = -pedantic -Wall
endif

# general options
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
#BOUNDS     = -DCHECK_BOUNDS
TERMINAL   = -DUNIX_TERMINAL
DATE_DEF   = -DBUILD_DATE="$(BUILD_DATE)"
SOLVE      = -DSOLVE_SVD
CC_PERL     = `perl -MExtUtils::Embed -e ccopts`
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE) $(TERMINAL) $(CC_PERL)
LIBS_PERL   = `perl -MExtUtils::Embed -e ldopts`
LIBS       = -lLinearAlgebra -lMyMath -lm
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(LIBS_PERL)


# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
SUBDIRS = ../Utility ../Lattices ../Modes ../Potentials
#
# Defines the source files needed #########################################
include ../Utility/Makefile
include ../Lattices/Makefile
include ../Modes/Makefile
include ../Potentials/Makefile

SOURCES = ./energy.cpp
#

# Defines the object files to be created using the SYSTEM name ############
OBJECTS = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))
#

PROGRAMS = $(OBJECTS)

all: $(PROGRAMS)

# include the dependicies
include $(addprefix depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        -o ./lib/$(notdir $@) \
        -c $(filter %/$(*),$(basename $(SOURCES)))$(SRC_EXT)


# Defines the Program Target (linking command) ############################
#Rank1Convex:  $(RK1OBJECTS) $(BSCOBJECTS)
#	      $(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/builddate.o -c builddate.cpp
#	      $(PURE) $(CC) $(LFLAGS) \
#	      -L$(LIB_LOC) -o ./bin/Rank1Convex \
#	      $(addprefix lib/,$(RK1OBJECTS)) \
#	      $(addprefix lib/,$(BSCOBJECTS)) $(LIBS) ./lib/builddate.o
#	      /bin/rm -f ./lib/builddate.o

# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	-/bin/rm -f $(addprefix lib/,$(OBJECTS:$(OBJ_EXT)=.*)) \
                  core *~ \
                  depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCC) -MM -I. $(addprefix -I ,$(SUBDIRS)) \
	      -I $(INCL_LOC)/ $(CCFLAGS) $< \
	      | sed '\''s/\($(notdir $*)\)\.o[ :]*/\1.o $(notdir $*).d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"

