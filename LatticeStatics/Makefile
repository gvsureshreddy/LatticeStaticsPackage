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

# Lattice include options (comment out to exclude the Lattice type)
#USE_KIM  = -DUSE_KIM
#USE_FEAP = -DUSE_FEAP
#USE_QC   = -DUSE_QC

# various compile options
CC_PERL     = `$(PERL) -MExtUtils::Embed -e ccopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
KIM_CONFIG_HELPER = kim-api-build-config
CC_KIM     = $(shell $(KIM_CONFIG_HELPER) --includes)
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE) $(TERMINAL) $(CC_PERL) $(CC_KIM)
LIBS_PERL   = `$(PERL) -MExtUtils::Embed -e ldopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
FEAP_OK := $(shell if test -e $(FEAPHOME8_3)/Feap8_3.a; then echo YES; else echo NO; fi)
ifneq (,$(strip $(USE_FEAP)))
  ifeq (NO,$(strip $(FEAP_OK)))
    $(error You must set the environment variable FEAPHOME8_3 and compile the library)
  endif
  LIBS_FEAP = $(FEAPHOME8_3)/Feap8_3.a -lgfortran -L/usr/X11R6/lib -lX11
else
  LIBS_FEAP = 
endif
ifneq (,$(strip $(USE_KIM)))
  ifeq ($(shell $(KIM_CONFIG_HELPER) --version),)
    $(error $(KIM_CONFIG_HELPER) utility is not available.  Something is wrong with your KIM API package setup)
  endif
  LIBS_KIM = $(shell $(KIM_CONFIG_HELPER) --ldflags) $(shell $(KIM_CONFIG_HELPER) --ldlibs)
else
  LIBS_KIM =
endif
LIBS       = -lLinearAlgebra -lMyMath -lm $(LIBS_PERL) $(LIBS_FEAP) $(LIBS_KIM) -lgfortran
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(DATE_DEF) $(LIBS_PERL)


# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
SUBDIRS    = Lattices Restrictions SolutionMethods Potentials Utility gpl_scripts

# Defines the source files needed #########################################
include $(addsuffix /Makefile, $(SUBDIRS))
ifeq (,$(strip $(USE_FEAP)))
  LATTICES := $(filter-out FEAP.cpp,$(LATTICES))
endif
ifeq (,$(strip $(USE_KIM)))
  LATTICES := $(filter-out MultiLatticeKIM.cpp,$(LATTICES))
endif
ifeq (,$(strip $(USE_QC)))
  LATTICES := $(filter-out QC.cpp,$(LATTICES))
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
        $(USE_FEAP) $(USE_KIM) $(USE_QC) \
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
                  lib/$(EXE).* bin/$(EXE) lib/FEAP-stub.o core *~ \
		  depend/*.d
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs clean

# Define rule for dependicy files
depend/$(notdir %).d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCC) -MM -I. $(addprefix -I ,$(SUBDIRS)) \
	      -I $(INCL_LOC)/ $(CCFLAGS) $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"

