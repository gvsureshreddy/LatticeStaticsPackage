#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend

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
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DEBUG) $(PROFILE) $(OPTIMIZE)

# Defines the source files needed #########################################
SOURCES  = Polynomial.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS)
.PHONY: all

Mathtest : $(OBJECTS) Mathtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./bin/Mathtest ./Mathtest.cpp \
		./lib/Polynomial.o 

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: install uninstall
# Define install
install: $(EXE)
	ar -r $(LIB_LOC)/libMath.a $(addprefix ./lib/,$(OBJECTS)); ranlib $(LIB_LOC)/libMath.a
	/bin/cp Math.h $(INCL_LOC)

# Define uninstall
uninstall:
	ar -d $(LIB_LOC)/libMath.a Math.o; ranlib $(LIB_LOC)/libMath.a
	/bin/rm -f $(INCL_LOC)/Math.h


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix ./lib/,$(OBJECTS)) core *~ \
		  Mathtest ./depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(CC) -MM $(CCFLAGS) -I. -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'
	@echo "Building $@"

