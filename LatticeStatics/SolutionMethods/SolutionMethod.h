#ifndef __SolutionMethod
#define __SolutionMethod

#include <fstream>

using namespace std;

class Lattice;

class SolutionMethod
{
public:
   virtual ~SolutionMethod() {}
   
   virtual int AllSolutionsFound() = 0;
   virtual double FindNextSolution(int &good) = 0;
   virtual int FindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				 char *datafile,const char *prefix,int Width,fstream &out) = 0;
};

#endif
