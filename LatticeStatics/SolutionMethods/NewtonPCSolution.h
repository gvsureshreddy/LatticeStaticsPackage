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
   double h;
   
public:
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,const Vector &one,
		    int Echo=1);
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,char *startfile,
		    fstream &out,int Echo);
   ~NewtonPCSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual Vector Force(Vector &D);
   virtual double Norm(Vector &D);
   virtual Vector tang(Vector &D);
   virtual Matrix MPI(Vector &D);
   virtual int BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out);
   
};

#endif
