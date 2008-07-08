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
   
   virtual int AllSolutionsFound() const = 0;
   virtual int FindNextSolution() = 0;
   virtual int FindCriticalPoint(Lattice* const Lat,PerlInput const& Input,int const& Width,
                                 fstream& out) = 0;
};

#endif
