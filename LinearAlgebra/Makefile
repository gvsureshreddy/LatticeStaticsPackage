#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend ./bin

SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
HDR_EXT    = .h

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
SOURCES  = Matrix.cpp Vector.cpp Vector3D.cpp CMatrix.cpp CVector.cpp

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
	ar -r $(LIB_LOC)/libLinearAlgebra.a $(addprefix ./lib/,$(OBJECTS)); ranlib $(LIB_LOC)/libLinearAlgebra.a 
	/bin/cp $(SOURCES:$(SRC_EXT)=$(HDR_EXT)) $(INCL_LOC)
	/bin/cp bin/EigenSys $(BIN_LOC)

# Define uninstall
uninstall:
	ar -d $(LIB_LOC)/libLinearAlgebra.a $(OBJECTS); ranlib $(LIB_LOC)/libLinearAlgebra.a
	/bin/rm -f $(addprefix $(INCL_LOC)/, $(SOURCES:$(SRC_EXT)=$(HDR_EXT)))
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

