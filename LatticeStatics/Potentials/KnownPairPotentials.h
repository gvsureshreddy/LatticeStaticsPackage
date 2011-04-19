#ifndef RSE__KnownPairPotentials
#define RSE__KnownPairPotentials

#include "PerlInput.h"
#include "RandomAlloy.h"
#include "LJ.h"
#include "LJConstCutoff.h"
#include "LJLinearCutoff.h"
#include "LJQuadraticCutoff.h"
#include "LJSplineCutoff.h"
#include "LJDobson.h"
#include "RadiiMorse.h"
#include "RadiiMorseCutoff.h"
#include "RadiiMorseCutoff2.h"
#include "RadiiMorseOriginal.h"
#include "GVMorse.h"
#include "TempMorse.h"
#include "Dobson.h"

PairPotentials* InitializePairPotential(PerlInput::HashStruct const& ParentHash,
                                        PerlInput const& Input, int const& i, int const& j);
PairPotentials* InitializePairPotential(char const* const HashName, PerlInput const& Input,
                                        int const& i, int const& j);
void UpdatePairPotential(PerlInput::HashStruct const& ParentHash, PerlInput const& Input,
                         int const& i, int const& j, PairPotentials* const Potential);
void UpdatePairPotential(char const* const HashName, PerlInput const& Input, int const& i,
                         int const& j, PairPotentials* const Potential);

#endif

