#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "RadiiMorse.h"
#include "TempMorse.h"

#define NOPOTENTIALS 2

PairPotentials *InitializePairPotential(char *datafile,const char *prefix,int i,int j);

#endif
