#**************************************************************************
#**************************************************************************

VPATH = ./lib ./depend

UNAME      = $(shell uname)

SRC_EXT    = .cpp
OBJ_EXT    = .o
DEP_EXT    = .d
BUILD_DATE = \"`date +%c;` from git commit `git rev-parse HEAD`\"

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
   #ANSI       = -A
   CC         = $(gCC)
   OPTIMIZE   = -O
   ANSI       = -pedantic -Wall
endif

# general options
#PURE       = purify
#PROFILE    = -pg
#DEBUG      = -g
#BOUNDS     = -DCHECK_BOUNDS
DATE_DEF   = -DBUILD_DATE="$(BUILD_DATE)"
CCFLAGS    = $(ANSI) $(DEBUG) $(PROFILE) $(OPTIMIZE) $(BOUNDS)

# File Locations ##########################################################
BIN_LOC    = ${HOME}/bin
INCL_LOC   = ${HOME}/include
LIB_LOC    = ${HOME}/lib
LFLAGS     = $(DATE_DEF)

# Defines the source files needed #########################################
SOURCES  = Polynomial.cpp MyComplexDouble.cpp

# Defines the object files to be created using the SYSTEM name ############
OBJECTS  = $(SOURCES:$(SRC_EXT)=$(OBJ_EXT))

all : $(OBJECTS)
.PHONY: all

Mathtest : $(OBJECTS) Mathtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@.o -c ./Mathtest.cpp
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./bin/Mathtest ./lib/Mathtest.o \
		./lib/Polynomial.o ./lib/MyComplexDouble.o

# include the dependicies
include $(addprefix ./depend/,$(OBJECTS:$(OBJ_EXT)=$(DEP_EXT)))

# Define object file dependicies and make rule
$(OBJECTS): Makefile
	$(CC) $(CCFLAGS) -I. -I $(INCL_LOC)/ -o ./lib/$@ -c $*$(SRC_EXT)

.PHONY: uninstall
# Define install
install: $(LIB_LOC)/libMyMath.a $(INCL_LOC)/MyMath.h $(INCL_LOC)/MyComplexDouble.h

$(INCL_LOC)/MyMath.h: MyMath.h
	/bin/cp MyMath.h $(INCL_LOC)/
$(INCL_LOC)/MyComplexDouble.h: MyComplexDouble.h
	/bin/cp MyComplexDouble.h $(INCL_LOC)/
$(LIB_LOC)/libMyMath.a: $(OBJECTS)
	$(CC) $(CCFLAGS) $(DATE_DEF) -o ./lib/MyMathBuildDate.o -c MyMathBuildDate.cpp
	ar -r $(LIB_LOC)/libMyMath.a $(addprefix ./lib/,$(OBJECTS)) ./lib/MyMathBuildDate.o
	/bin/rm -f ./lib/MyMathBuildDate.o
	ranlib $(LIB_LOC)/libMyMath.a

# Define uninstall
uninstall:
	/bin/rm -f $(LIB_LOC)/libMyMath.a
	/bin/rm -f $(INCL_LOC)/MyMath.h
	/bin/rm -f $(INCL_LOC)/MyComplexDouble.h


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
