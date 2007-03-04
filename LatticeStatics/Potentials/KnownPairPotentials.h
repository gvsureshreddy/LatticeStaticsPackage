#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "LJ.h"
#include "LJCutoff.h"
#include "RadiiMorse.h"
#include "RadiiMorseCutoff.h"
#include "TempMorse.h"
#include "Dobson.h"

#define NOPOTENTIALS 6

PairPotentials *InitializePairPotential(char *datafile,const char *prefix,int i,int j);
void UpdatePairPotential(char *datafile,const char *prefix,int i,int j,
			 PairPotentials *Potential);

#endif
