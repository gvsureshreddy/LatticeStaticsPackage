#ifndef __KnownSolutionMethod
#define __KnownSolutionMethod

#include "ArcLengthSolution.h"
#include "ScanningSolution.h"

#define KNOWNSOLUTIONMETHODS 2


SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width,
				   int Echo);

#endif
