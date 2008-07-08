#ifndef __KnownPairPotentials
#define __KnownPairPotentials

#include "PerlInput.h"
#include "LJ.h"
#include "LJCutoff.h"
#include "RadiiMorse.h"
#include "RadiiMorseCutoff.h"
#include "RadiiMorseCutoff2.h"
#include "TempMorse.h"
#include "Dobson.h"

#define LINELENGTH 600

PairPotentials* InitializePairPotential(PerlInput::HashStruct const& ParentHash,
                                        PerlInput const& Input,int const& i,int const& j);
PairPotentials* InitializePairPotential(char const* const HashName,PerlInput const& Input,
                                        int const& i,int const& j);
void UpdatePairPotential(PerlInput::HashStruct const& ParentHash,PerlInput const& Input,
                         int const& i,int const& j,PairPotentials* const Potential);
void UpdatePairPotential(char const* const HashName,PerlInput const& Input,int const& i,
                         int const& j,PairPotentials* const Potential);

#endif
