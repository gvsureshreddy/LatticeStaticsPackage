#ifndef __KnownModes
#define __KnownModes

#include "UniDefTemp2DExpand.h"
#include "ML_Expand.h"
#include "ML_NaCl.h"
#include "ML_Rhombo.h"

#define KNOWNMODES 4

LatticeMode *InitializeMode(Lattice *Lat,char *datafile,const char *prefix);

#endif
