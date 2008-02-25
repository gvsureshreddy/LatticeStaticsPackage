#ifndef __KnownSolutionMethod
#define __KnownSolutionMethod

#include "ArcLengthSolution.h"
#include "ScanningSolution.h"
#include "NewtonPCSolution.h"
#include "NewtonUpdatePCSolution.h"
#include "NewtonQRUpdatePCSolution.h"

#define KNOWNSOLUTIONMETHODS 4


SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
                                   char *startfile,Lattice *Lat,fstream &out,int Width,
                                   int Echo);

#endif
