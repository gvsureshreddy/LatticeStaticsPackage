#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend ./bin

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
BOUNDS     = -DCHECK_BOUNDS
OPTIMIZE   = -O
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE)

# Defines the source files needed #########################################
SOURCES  = Matrix.cpp Vector.cpp Vector3D.cpp CMatrix.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS) EigenSys matrixtest
.PHONY: all

matrixtest : $(OBJECTS) matrixtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o bin/matrixtest ./matrixtest.cpp \
	$(addprefix ./lib/,$(OBJECTS))

EigenSys : $(OBJECTS) EigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o bin/EigenSys ./EigenSys.cpp ./lib/Matrix.o

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: install uninstall
# Define install
install: $(EXE)
	ar -r $(LIB_LOC)/libLinearAlgebra.a ./lib/Vector.o ./lib/Vector3D.o ./lib/Matrix.o ./lib/CMatrix.o; ranlib $(LIB_LOC)/libLinearAlgebra.a 
	/bin/cp Vector.h Vector3D.h Matrix.h CMatrix.h $(INCL_LOC)
	/bin/cp bin/EigenSys $(BIN_LOC)

# Define uninstall
uninstall:
	ar -d $(LIB_LOC)/libLinearAlgebra.a Vector.o Vector3D.o Matrix.o CMatrix.o; ranlib $(LIB_LOC)/libLinearAlgebra.a
	/bin/rm -f $(INCL_LOC)/Vector.h $(INCL_LOC)/Vector3D.h $(INCL_LOC)/Matrix.h $(INCL_LOC)/CMatrix.h
	/bin/rm -f $(BIN_LOC)/EigenSys


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix ./lib/,$(OBJECTS)) core *~ \
		  bin/matrixtest bin/EigenSys ./depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(CC) -MM $(CCFLAGS) -I. -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'
	@echo "Building $@"

