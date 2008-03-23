#ifndef __SolutionMethod
#define __SolutionMethod

#include <fstream>
#include "PerlInput.h"

using namespace std;

class Lattice;

class SolutionMethod
{
public:
   virtual ~SolutionMethod() {}
   
   virtual int AllSolutionsFound() = 0;
   virtual int FindNextSolution() = 0;
   virtual int FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out) = 0;
};

#endif
