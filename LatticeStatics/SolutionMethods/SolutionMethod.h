#ifndef __SolutionMethod
#define __SolutionMethod

#include <fstream.h>

class Lattice;

class SolutionMethod
{
public:
   virtual ~SolutionMethod() {}

   virtual int AllSolutionsFound() = 0;
   virtual double FindNextSolution(int &good) = 0;
   virtual int BisectAlert(Lattice *Lat,int Width,fstream &out) = 0;

};

#endif
