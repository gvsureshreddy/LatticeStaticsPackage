#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "LJ.h"
#include "RadiiMorse.h"
#include "TempMorse.h"
#include "Dobson.h"

#define NOPOTENTIALS 4

PairPotentials *InitializePairPotential(char *datafile,const char *prefix,int i,int j);
void UpdatePairPotential(char *datafile,const char *prefix,int i,int j,
			 PairPotentials *Potential);

#endif
