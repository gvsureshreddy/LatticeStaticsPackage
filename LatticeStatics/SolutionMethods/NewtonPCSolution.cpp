#include "NewtonPCSolution.h"


#include "UtilityFunctions.h"

using namespace std;

#define CLOSEDDEFAULT 30

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,
				   const char *prefix,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{   
}

int NewtonPCSolution::AllSolutionsFound()
{
   return (CurrentSolution_ >= NumSolutions_);
}

double NewtonPCSolution::FindNextSolution(int &good)
{
}

int NewtonPCSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
				   int Width,fstream &out)
{
}
