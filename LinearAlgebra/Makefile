#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend ./bin

SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
HDR_EXT    = .h
BUILD_DATE = \"`date +%c`\"

# Compiling options, etc. #################################################
#
CC         = g++
ANSI       = -ansi -pedantic
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
BOUNDS     = -DCHECK_BOUNDS
DATE_DEF   = -DBUILD_DATE="$(BUILD_DATE)"
OPTIMIZE   = -O
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE) $(DATE_DEF)
#
# the below lines are for OS X -- a bug in the libraries...
#
UNAME      = $(shell uname)
#
ifeq ($(UNAME),Darwin)
   LFLAGS     := $(LFLAGS) `ar -x /usr/lib/libstdc++.a lib/dcomio.o lib/dcomplex.o` lib/dcomio.o lib/dcomplex.o
endif
#
#
# Defines the source files needed #########################################
SOURCES  = Matrix.cpp Vector.cpp Vector3D.cpp CMatrix.cpp CVector.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS) EigenSys matrixtest
.PHONY: all

matrixtest : $(OBJECTS) matrixtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -o bin/matrixtest ./matrixtest.cpp \
	$(addprefix ./lib/,$(OBJECTS))

EigenSys : $(OBJECTS) EigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -o bin/EigenSys ./EigenSys.cpp ./lib/Matrix.o

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: install uninstall
# Define install
install: all 
	$(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/LinearAlgebraBuildDate.o -c LinearAlgebraBuildDate.cpp
	ar -r $(LIB_LOC)/libLinearAlgebra.a $(addprefix ./lib/,$(OBJECTS)) ./lib/LinearAlgebraBuildDate.o
	/bin/rm -f ./lib/LinearAlgebraBuildDate.o
	ranlib $(LIB_LOC)/libLinearAlgebra.a 
	/bin/cp $(SOURCES:$(SRC_EXT)=$(HDR_EXT)) $(INCL_LOC)
	/bin/cp bin/EigenSys $(BIN_LOC)

# Define uninstall
uninstall:
	/bin/rm -f $(LIB_LOC)/libLinearAlgebra.a
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

