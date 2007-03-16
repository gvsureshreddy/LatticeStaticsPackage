#ifndef __NewtonPCSolution
#define __NewtonPCSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class NewtonPCSolution : public SolutionMethod
{
private:
   LatticeMode *Mode_;
   
   int CurrentSolution_;
   int Echo_;
   int NumSolutions_;

public:
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,int Echo=1);
   ~NewtonPCSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out);
   
};

#endif
