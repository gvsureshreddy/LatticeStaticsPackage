#ifndef RSE__KnownSolutionMethod
#define RSE__KnownSolutionMethod

#include "PerlInput.h"
#include "ArcLengthSolution.h"
#include "ScanningSolution.h"
#include "NewtonPCSolution.h"

SolutionMethod* InitializeSolution(LatticeMode* const Mode,PerlInput const& Input,
                                   Lattice* const Lat,fstream& out,int const& Width,
                                   int const& Echo);

#endif
