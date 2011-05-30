#ifndef RSE__SolutionMethod
#define RSE__SolutionMethod

#include "PerlInput.h"

using namespace std;

class Lattice;

class SolutionMethod
{
public:
   virtual ~SolutionMethod()
   {
   }

   virtual int AllSolutionsFound() const = 0;
   virtual int FindNextSolution() = 0;
   virtual void FindCriticalPoint(Lattice* const Lat, int* const TotalNumCPCrossings,
                                  PerlInput const& Input, int const& Width, ostream& out) = 0;
   virtual char const* const Type() const = 0;
};

#endif
