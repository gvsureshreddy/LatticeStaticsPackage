#**************************************************************************
#**************************************************************************

VPATH = bin lib depend Lattices Modes SolutionMethods Potentials Utility

EXE        = LatticeStatics
SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d

# Compiling options, etc. #################################################
#
CC         = g++
ANSI       = -pedantic
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
OPTIMIZE   = -O
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LIBS       = -lLinearAlgebra -lMath -lm
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE)
SUBDIRS    = Lattices Modes SolutionMethods Potentials Utility

# Defines the source files needed #########################################
include $(addsuffix /Makefile, $(SUBDIRS))
SOURCES  = LatticeStatics.cpp\
	   $(addprefix Lattices/,$(LATTICES))\
	   $(addprefix Modes/,$(MODES))\
	   $(addprefix SolutionMethods/,$(SOLUTIONMETHODS))\
	   $(addprefix Potentials/,$(POTENTIALS))\
	   $(addprefix Utility/,$(UTILITY))


# Defines the object files to be created using the SYSTEM name ############
STRIPEDSOURCES = $(notdir $(SOURCES))
OBJECTS  = $(STRIPEDSOURCES:$(SRC_EXT)=$(OBJ_EXT))

# Defines the Program Target (linking command) ############################
$(EXE):  $(OBJECTS)
	 $(PURE) $(CC) $(LFLAGS) -L$(LIB_LOC) -o ./bin/$(EXE) \
		 $(addprefix lib/,$(OBJECTS)) $(LIBS) 

# include the dependicies
include $(addprefix depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) -I $(INCL_LOC)/ \
        -o ./lib/$(notdir $@) \
        -c $(filter %$(*),$(basename $(SOURCES)))$(SRC_EXT)


.PHONY: install uninstall
# Define install
install: $(EXE)
	/bin/cp bin/$(EXE) $(BIN_LOC)/
	/bin/cp Utility/UniDef_Symm_Solution $(BIN_LOC)/
	/bin/cp gpl_scripts/NiTi2gpl $(BIN_LOC)/

# Define uninstall
uninstall:
	/bin/rm -f $(BIN_LOC)/$(EXE)
	/bin/rm -f $(BIN_LOC)/UniDef_Symm_Solution
	/bin/rm -f $(BIN_LOC)/NiTi2gpl



# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix lib/,$(OBJECTS)) \
		  $(addsuffix /*~,$(SUBDIRS))\
                  bin/$(EXE) core *~ \
		  depend/*.d

# Define rule for dependicy files
depend/$(notdir %).d: %$(SRC_EXT)
	@$(SHELL) -ec '$(CC) -MM $(CCFLAGS) -I. $(addprefix -I ,$(SUBDIRS)) \
	      -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' \
	      > $@; [ -s $@ ] || rm -f $@'
	@echo "Building $@"
