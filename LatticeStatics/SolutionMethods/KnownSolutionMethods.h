#ifndef __KnownSolutionMethod
#define __KnownSolutionMethod

#include "ArcLengthSolution.h"
#include "ScanningSolution.h"
#include "NewtonPCSolution.h"

#define KNOWNSOLUTIONMETHODS 3


SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width,
				   int Echo);

#endif
