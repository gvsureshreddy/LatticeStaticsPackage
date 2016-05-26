#**************************************************************************
#**************************************************************************

VPATH = bin lib depend Lattices Restrictions SolutionMethods Potentials Utility

UNAME      = $(shell uname)

EXE        = LatticeStatics
LIBRARY    = libBFBSymPac.a
SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
BUILD_DATE = \"`date +%c;` from git commit `git rev-parse HEAD`\"

# Compiling options, etc. #################################################
#
gCXX        = g++

ifeq ($(UNAME),Darwin)
   # gCC compiler
   CXX         = $(gCXX)
   OPTIMIZE   = -O
   ANSI       = -pedantic -Wall
endif

ifeq ($(UNAME),Linux)
   # PGI compiler
   #CXX         = pgCC
   #OPTIMIZE   = -fast -Minline=levels:10 --no_exceptions
   #ANSI       = -A -Minform=inform
   CXX	= $(gCXX)
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
#USE_NEO  = -DUSE_NEO
#USE_ELA  = -DUSE_ELA
#USE_KIM  = -DUSE_KIM
#USE_FEAP = -DUSE_FEAP
#USE_QC   = -DUSE_QC

# various compile options
CXX_PERL     = `$(PERL) -MExtUtils::Embed -e ccopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
LIBS_PERL   = `$(PERL) -MExtUtils::Embed -e ldopts | sed -e 's/-arch[^-]*//g' | sed -e 's/-Wdeclaration-after-statement//g'`
ifneq (,$(strip $(USE_NEO)))
  LFLAGS_NEO = -Lneo-build
  LIBS_NEO = -lneo-hookean
else
  LFLAGS_NEO =
  LIBS_NEO =
endif
ifneq (,$(strip $(USE_ELA)))
  LFLAGS_ELA = -Lelastica-build
  LIBS_ELA = -lelastica-beam
else
  LFLAGS_ELA =
  LIBS_ELA =
endif
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
  KIM_CONFIG_HELPER = kim-api-build-config
  ifeq ($(shell $(KIM_CONFIG_HELPER) --version),)
    $(error $(KIM_CONFIG_HELPER) utility is not available.  Something is wrong with your KIM API package setup)
  endif
  CXX_KIM     = $(shell $(KIM_CONFIG_HELPER) --includes)
  LIBS_KIM = $(shell $(KIM_CONFIG_HELPER) --ldflags) $(shell $(KIM_CONFIG_HELPER) --ldlibs)
else
  CXX_KIM =
  LIBS_KIM =
endif
LIBS       = -lLinearAlgebra -lMyMath -lm $(LIBS_PERL) $(LIBS_FEAP) $(LIBS_KIM) $(LIBS_NEO) $(LIBS_ELA)
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(DATE_DEF) $(LIBS_PERL) $(LFLAGS_NEO) $(LFLAGS_ELA)
CXXFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE) $(TERMINAL) $(CXX_PERL) $(CXX_KIM)


# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
SUBDIRS    = Lattices Restrictions SolutionMethods Potentials Utility gpl_scripts

# Defines the source files needed #########################################
include $(addsuffix /Makefile, $(SUBDIRS))
ifeq (,$(strip $(USE_NEO)))
  LATTICES := $(filter-out NeoHookean2D.cpp,$(LATTICES))
endif
ifeq (,$(strip $(USE_ELA)))
  LATTICES := $(filter-out ElasticaBeam1D.cpp,$(LATTICES))
endif
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


$(EXE):  $(OBJECTS) $(if $(USE_NEO),neo-build,) $(if $(USE_ELA),elastica-build,)
	 $(CXX) $(CXXFLAGS) $(DATE_DEF) -o lib/builddate.o -c Utility/builddate.cpp
	 $(PURE) $(CXX) $(LFLAGS) lib/builddate.o \
		 -L$(LIB_LOC) -o ./bin/$(EXE) \
		 $(addprefix lib/,$(OBJECTS)) $(LIBS)
	 /bin/rm -f lib/builddate.o

# include the dependicies
include $(addprefix depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CXX) $(CXXFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        $(USE_FEAP) $(USE_KIM) $(USE_QC) $(USE_NEO) $(USE_ELA)\
        -o ./lib/$(notdir $@) \
        -c $(filter %/$(*),$(basename $(SOURCES)))$(SRC_EXT)

KnownPotentials.o: Makefile
	echo $(patsubst %.h,%,$(PTNLS))
	
.PHONY: elastica-build

elastica-build: elastica-beam/elastica-beam.cpp elastica-beam/CMakeLists.txt
	if test ! -d elastica-build; then mkdir elastica-build; fi
	cd elastica-build && cmake -DCREATE_LIB:BOOL=ON ../elastica-beam && make

.PHONY: neo-build
	
neo-build: neo-hookean/neo-hookean.cpp neo-hookean/CMakeLists.txt
	if test ! -d neo-build; then mkdir neo-build; fi
	cd neo-build && cmake -DCREATE_LIB:BOOL=ON ../neo-hookean && make

.PHONY: install uninstall
# Define install
install: $(BIN_LOC)/$(EXE) $(addprefix $(BIN_LOC)/,$(GPLSCRIPTS)) $(LIB_LOC)/$(LIBRARY) Utilities
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs install

$(LIB_LOC)/$(LIBRARY): $(filter-out LatticeStatics.o,$(OBJECTS)) Utility/BFBSymPacCInterface.cpp
	$(CXX) $(CXXFLAGS) $(DATE_DEF) -o ./lib/builddate.o -c Utility/builddate.cpp
	$(CXX) $(CXXFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/  -o ./lib/BFBSymPacCInterface.o -c Utility/BFBSymPacCInterface.cpp
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
	-/bin/rm -rf neo-build
	-/bin/rm -rf elastica-build
	$(MAKE) -e -C ./Utility -f MakeUtilityProgs clean

# Define rule for dependicy files
depend/$(notdir %).d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCXX) -MM -I. $(addprefix -I ,$(SUBDIRS)) \
	      -I $(INCL_LOC)/ $(CXXFLAGS) $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"
