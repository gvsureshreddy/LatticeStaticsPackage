#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "LJ.h"
#include "RadiiMorse.h"
#include "TempMorse.h"

#define NOPOTENTIALS 3

PairPotentials *InitializePairPotential(char *datafile,const char *prefix,int i,int j);
void UpdatePairPotential(char *datafile,const char *prefix,int i,int j,
			 PairPotentials *Potential);

#endif
