#ifndef __NetwonPCSolution
#define __NetwonPCSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class NetwonPCSolution : public SolutionMethod
{
private:
   int CurrentSolution_;
   int Echo_;
   int NumSolutions_;

public:
   NetwonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,int Echo=1);
   ~NetwonPCSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out);
   
};

#endif
