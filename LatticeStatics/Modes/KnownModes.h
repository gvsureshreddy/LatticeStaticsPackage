#ifndef __KnownModes
#define __KnownModes

#include "ML_Expand.h"
#include "ML_NaCl.h"
#include "ML_Rhombo.h"
#include "ML_B19.h"

#define KNOWNMODES 4

LatticeMode *InitializeMode(Lattice *Lat,char *datafile,const char *prefix);

#endif
