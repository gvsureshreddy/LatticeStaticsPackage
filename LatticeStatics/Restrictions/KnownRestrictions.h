#ifndef RSE__KnownRestrictions
#define RSE__KnownRestrictions

#include "PerlInput.h"
#include "RestrictToTranslatedSubSpace.h"
#include "RestrictToTranslatedSubSpaceOld.h"
#include "NoRestriction.h"

Restriction* const InitializeRestriction(Lattice* const Lat,PerlInput const& Input);

#endif
