#include "NewtonPCSolution.h"


#include "UtilityFunctions.h"

using namespace std;

#define CLOSEDDEFAULT 30

NetwonPCSolution::NetwonPCSolution(LatticeMode *Mode,char *datafile,
				   const char *prefix,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{   
}

int NetwonPCSolution::AllSolutionsFound()
{
   return (CurrentSolution_ >= NumSolutions_);
}

double NetwonPCSolution::FindNextSolution(int &good)
{
}

int NetwonPCSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
				   int Width,fstream &out)
{
}
