#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "PerlInput.h"
#include "LJ.h"
#include "LJCutoff.h"
#include "RadiiMorse.h"
#include "RadiiMorse2.h"
#include "RadiiMorseCutoff.h"
#include "RadiiMorseCutoff2.h"
#include "TempMorse.h"
#include "Dobson.h"

#define LINELENGTH 600

PairPotentials *InitializePairPotential(PerlInput::HashStruct ParentHash,PerlInput &Input,
                                        int i,int j);
PairPotentials *InitializePairPotential(char *HashName,PerlInput &Input,int i,int j);
void UpdatePairPotential(PerlInput::HashStruct ParentHash,PerlInput &Input,
                         int i,int j,PairPotentials *Potential);
void UpdatePairPotential(char *HashName,PerlInput &Input,int i,int j,PairPotentials *Potential);

#endif
