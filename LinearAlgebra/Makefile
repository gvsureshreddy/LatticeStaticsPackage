#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend ./bin

UNAME      = $(shell uname)

SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
HDR_EXT    = .h
BUILD_DATE = \"`date +%c`\"

# Compiling options, etc. #################################################
#
gCC        = g++

ifeq ($(UNAME),Darwin)
   # gCC compiler
   CC         = $(gCC)
   OPTIMIZE   = -O
   ANSI       = -pedantic
endif

ifeq ($(UNAME),Linux)
   # PGI compiler
   CC         = pgCC
   OPTIMIZE   = -fast -Minline=levels:10 --no_exceptions
   #ANSI       = -A
endif

# general options 
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
BOUNDS     = -DCHECK_BOUNDS
DATE_DEF   = -DBUILD_DATE="$(BUILD_DATE)"
SOLVE      = -DSOLVE_SVD
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS) $(SOLVE)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LIBS       = -lMyMath -lm
LFLAGS     = $(DATE_DEF)

# Defines the source files needed #########################################
SOURCES  = Matrix.cpp Vector.cpp Vector3D.cpp CMatrix.cpp CVector.cpp SparseMatrix.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

PROGRAMS = EigenSys HermiteEigenSys PrincipalSys

all : $(OBJECTS) $(PROGRAMS) matrixtest
.PHONY: all

matrixtest : $(OBJECTS) matrixtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./matrixtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -L$(LIB_LOC) -o ./bin/matrixtest \
	./lib/matrixtest.o $(addprefix ./lib/,$(OBJECTS)) $(LIBS)

EigenSys : $(OBJECTS) EigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./EigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -o ./bin/EigenSys \
	./lib/EigenSys.o ./lib/Matrix.o 

HermiteEigenSys : $(OBJECTS) HermiteEigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./HermiteEigenSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -L$(LIB_LOC) -o ./bin/HermiteEigenSys \
	./lib/HermiteEigenSys.o ./lib/CMatrix.o ./lib/Matrix.o $(LIBS)

PrincipalSys : $(OBJECTS) PrincipalSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./PrincipalSys.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ $(LFLAGS) -o ./bin/PrincipalSys \
	./lib/PrincipalSys.o ./lib/Matrix.o ./lib/Vector.o

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: uninstall
# Define install
install: all $(LIB_LOC)/libLinearAlgebra.a $(addprefix $(BIN_LOC)/,$(PROGRAMS))

$(LIB_LOC)/libLinearAlgebra.a: $(OBJECTS)
	$(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/LinearAlgebraBuildDate.o -c LinearAlgebraBuildDate.cpp
	ar -r $(LIB_LOC)/libLinearAlgebra.a $(addprefix ./lib/,$(OBJECTS)) ./lib/LinearAlgebraBuildDate.o
	/bin/rm -f ./lib/LinearAlgebraBuildDate.o
	ranlib $(LIB_LOC)/libLinearAlgebra.a 
	/bin/cp $(SOURCES:$(SRC_EXT)=$(HDR_EXT)) $(INCL_LOC)

$(addprefix $(BIN_LOC)/,$(PROGRAMS)): $(BIN_LOC)/%: %
	/bin/cp ./bin/$(notdir $@) $(BIN_LOC)

# Define uninstall
uninstall:
	/bin/rm -f $(LIB_LOC)/libLinearAlgebra.a
	/bin/rm -f $(addprefix $(INCL_LOC)/, $(SOURCES:$(SRC_EXT)=$(HDR_EXT)))
	/bin/rm -f $(BIN_LOC)/EigenSys
	/bin/rm -f $(BIN_LOC)/HermiteEigenSys
	/bin/rm -f $(BIN_LOC)/PrincipalSys


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix ./lib/,$(OBJECTS)) \
	$(addprefix ./lib/,$(OBJECTS:$(OBJ_EXT)=.ti)) core *~ \
	lib/matrixtest.* lib/EigenSys.* lib/HermiteEigenSys.* lib/PrincipalSys.* \
	bin/matrixtest  bin/EigenSys bin/HermiteEigenSys bin/PrincipalSys ./depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCC) -MM -I. -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'
	@echo "Building $@"

