#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend

SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d

# Compiling options, etc. #################################################
#
CC         = g++
ANSI       = -ansi -pedantic
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
OPTIMIZE   = -O
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE)

# Defines the source files needed #########################################
SOURCES  = Matrix.cpp Vector.cpp Vector3D.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS)
.PHONY: all

matrixtest : $(OBJECTS)
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o matrixtest ./matrixtest.cpp ./lib/Matrix.o ./lib/Vector.o

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: install uninstall
# Define install
install: $(EXE)
	ar -r $(LIB_LOC)/libLinearAlgebra.a ./lib/Vector.o ./lib/Vector3D.o ./lib/Matrix.o 
	/bin/cp Vector.h Vector3D.h Matrix.h $(INCL_LOC)

# Define uninstall
uninstall:
	ar -d $(LIB_LOC)/libLinearAlgebra.a Vector.o Vector3D.o Matrix.o
	/bin/rm -f $(INCL_LOC)/Vector.h $(INCL_LOC)/Vector3D.h $(INCL_LOC)/Matrix.h


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix ./lib/,$(OBJECTS)) core *~ \
		  matrixtest ./depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(CC) -MM $(CCFLAGS) -I. -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'
	@echo "Building $@"

