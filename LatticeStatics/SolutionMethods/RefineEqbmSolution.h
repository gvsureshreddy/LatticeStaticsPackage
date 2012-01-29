#ifndef RSE__RefineEqbmSolution
#define RSE__RefineEqbmSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "Restriction.h"

using namespace std;

class Lattice;

class RefineEqbmSolution : public SolutionMethod
{
public:
   enum ConvergeType {Both, Force, Displacement};
private:
   Restriction* Restrict_;

   int Echo_;
   int SolutionFound_;
   int NumSolutions_;
   Vector* Guesses_;

   double Converge_;            // Convergence criteria
   ConvergeType ConvergeType_;  // Quantities to check for convergence

public:
   RefineEqbmSolution(Restriction* const Restrict, Vector const& one,
                      double const& Converge, ConvergeType CnvrgTyp, int const& Echo = 1);
   RefineEqbmSolution(Restriction* const Restrict, PerlInput const& Input, Vector const& one,
                      int const& Echo = 1);
   RefineEqbmSolution(Restriction* const Restrict, PerlInput const& Input, int const& Echo);
   ~RefineEqbmSolution()
   {
      delete[] Guesses_;
   }

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const
   {
      return SolutionFound_ >= NumSolutions_;
   }

   virtual int FindNextSolution(PerlInput const& Input, int const& Width, ostream& out);
   virtual void FindCriticalPoint(Lattice* const Lat, int* const TotalNumCPCrossings,
                                  PerlInput const& Input, int const& Width, ostream& out)
   {
   }

   virtual char const* const Type() const
   {
      return "RefineEqbmSolution";
   }
};

#endif
