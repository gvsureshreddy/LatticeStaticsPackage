#**************************************************************************
#**************************************************************************

VPATH = bin lib depend Lattices Restrictions SolutionMethods Potentials Utility

UNAME      = $(shell uname)

EXE        = LatticeStatics
LIBRARY    = libBFBSymPac.a
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
#SOLVE      = -DSOLVE_SVD
CC_PERL     = `perl -MExtUtils::Embed -e ccopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE) $(TERMINAL) $(CC_PERL)
LIBS_PERL   = `perl -MExtUtils::Embed -e ldopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
FEAP_OK := $(shell if test -e $(FEAPHOME8_3)/Feap8_3.a; then echo YES; else echo NO; fi)
ifeq (YES,$(strip $(FEAP_OK)))
  LIBS_FEAP = $(FEAPHOME8_3)/Feap8_3.a -lgfortran -L/usr/X11R6/lib -lX11
else
  LIBS_FEAP = 
endif
LIBS       = -lLinearAlgebra -lMyMath -lm
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(DATE_DEF) $(LIBS_PERL) $(LIBS_FEAP)


# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
SUBDIRS    = Lattices Restrictions SolutionMethods Potentials Utility gpl_scripts

# Defines the source files needed #########################################
include $(addsuffix /Makefile, $(SUBDIRS))
ifeq (YES,$(strip $(FEAP_OK)))
  LATTICES := $(filter-out FEAP-stub.cpp,$(LATTICES))
endif
SOURCES  = ./LatticeStatics.cpp\
	   $(addprefix Lattices/,$(LATTICES))\
	   $(addprefix Restrictions/,$(RESTRICTIONS))\
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
	 $(PURE) $(CC) $(LFLAGS) lib/builddate.o \
		 -L$(LIB_LOC) -o ./bin/$(EXE) \
		 $(addprefix lib/,$(OBJECTS)) $(LIBS)
	 /bin/rm -f lib/builddate.o

# include the dependicies
include $(addprefix depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        -o ./lib/$(notdir $@) \
        -c $(filter %/$(*),$(basename $(SOURCES)))$(SRC_EXT)

KnownPotentials.o: Makefile
	echo $(patsubst %.h,%,$(PTNLS))

.PHONY: install uninstall
# Define install
install: $(BIN_LOC)/$(EXE) $(addprefix $(BIN_LOC)/,$(GPLSCRIPTS)) $(LIB_LOC)/$(LIBRARY) Utilities
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs install

$(LIB_LOC)/$(LIBRARY): $(filter-out LatticeStatics.o,$(OBJECTS)) Utility/BFBSymPacCInterface.cpp
	$(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/builddate.o -c Utility/builddate.cpp
	$(CC) $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/  -o ./lib/BFBSymPacCInterface.o -c Utility/BFBSymPacCInterface.cpp
	ar -r $(LIB_LOC)/$(LIBRARY) $(addprefix ./lib/,$(filter-out QC-stub.o LatticeStatics.o,$(OBJECTS))) ./lib/builddate.o ./lib/BFBSymPacCInterface.o
	/bin/rm -f ./lib/builddate.o
	/bin/rm -f ./lib/BFBSymPacCInterface.o
	ranlib $(LIB_LOC)/$(LIBRARY)

$(BIN_LOC)/$(EXE):	$(EXE)
	/bin/cp bin/$(EXE) $(BIN_LOC)

$(addprefix $(BIN_LOC)/,$(GPLSCRIPTS)):	$(BIN_LOC)/%: gpl_scripts/%
	/bin/cp gpl_scripts/$(@F) $(BIN_LOC)

# Define uninstall
uninstall:
	/bin/rm -f $(BIN_LOC)/$(EXE)
	/bin/rm -f $(LIB_LOC)/$(LIBRARY)
	/bin/rm -f $(addprefix $(BIN_LOC)/,$(GPLSCRIPTS))
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
	      -I $(INCL_LOC)/ $(CCFLAGS) $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"

