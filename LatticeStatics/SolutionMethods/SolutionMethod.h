#ifndef RSE__SolutionMethod
#define RSE__SolutionMethod

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
   virtual void FindCriticalPoint(Lattice* const Lat,int& TotalNumCPCrossings,
                                  PerlInput const& Input,int const& Width,fstream& out) = 0;
};

#endif
