#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend

UNAME      = $(shell uname)

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
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DATE_DEF)

# Defines the source files needed #########################################
SOURCES  = Polynomial.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS)
.PHONY: all

Mathtest : $(OBJECTS) Mathtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./Mathtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./bin/Mathtest ./lib/Mathtest.o \
		./lib/Polynomial.o 

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: install uninstall
# Define install
install: all
	$(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/MyMathBuildDate.o -c MyMathBuildDate.cpp
	ar -r $(LIB_LOC)/libMyMath.a $(addprefix ./lib/,$(OBJECTS)) ./lib/MyMathBuildDate.o
	/bin/rm -f ./lib/MyMathBuildDate.o
	ranlib $(LIB_LOC)/libMyMath.a
	/bin/cp MyMath.h $(INCL_LOC)/

# Define uninstall
uninstall:
	/bin/rm -f $(LIB_LOC)/libMyMath.a
	/bin/rm -f $(INCL_LOC)/MyMath.h


# Removes the current systems set of object files and executable file
.PHONY: clean

clean:
	 -/bin/rm -f $(addprefix ./lib/,$(OBJECTS)) core *~ \
		  $(addprefix ./lib/,$(OBJECTS:$(OBJ_EXT)=.ti)) \
		  ./lib/Mathtest.* ./bin/Mathtest ./depend/*.d

# Define rule for dependicy files
depend/%.d: %$(SRC_EXT)
	@$(SHELL) -ec '$(gCC) -MM -I. -I $(INCL_LOC)/ $< \
	      | sed '\''s/\($*\)\.o[ :]*/\1.o $*.d : /g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'
	@echo "Building $@"

