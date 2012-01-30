#ifndef RSE__KnownSolutionMethod
#define RSE__KnownSolutionMethod

#include "PerlInput.h"
#include "RefineEqbmSolution.h"
#include "ArcLengthSolution.h"
#include "ScanningSolution.h"
#include "NewtonPCSolution.h"
#include "NEBSolution.h"
#include "ODSolution.h"

SolutionMethod* InitializeSolution(Restriction* const Restrict, PerlInput const& Input,
                                   Lattice* const Lat, ostream& out, int const& Width,
                                   int const& Echo);

#endif

