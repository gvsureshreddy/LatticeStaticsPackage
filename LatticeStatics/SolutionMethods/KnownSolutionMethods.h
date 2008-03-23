#ifndef __KnownSolutionMethod
#define __KnownSolutionMethod

#include "PerlInput.h"
#include "ArcLengthSolution.h"
#include "ScanningSolution.h"
#include "NewtonPCSolution.h"
#include "NewtonUpdatePCSolution.h"
#include "NewtonQRUpdatePCSolution.h"

SolutionMethod *InitializeSolution(LatticeMode *Mode,PerlInput &Input,Lattice *Lat,
                                   fstream &out,int Width,int Echo);

#endif
