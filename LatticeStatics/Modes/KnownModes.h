#ifndef __KnownModes
#define __KnownModes

#include "UniDefTemp2DExpand.h"
#include "ML_Expand.h"
#include "ML_NaCl.h"
#include "ML_Rhombo.h"
#include "ML_B19.h"

#define KNOWNMODES 5

LatticeMode *InitializeMode(Lattice *Lat,char *datafile,const char *prefix);

#endif
