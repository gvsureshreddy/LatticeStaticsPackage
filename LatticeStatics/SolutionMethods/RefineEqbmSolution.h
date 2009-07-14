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
   enum ConvergeType {Both,Force,Displacement};
private:
   Restriction *Restrict_;
   
   int Echo_;
   int SolutionFound_;

   double Converge_;            //Convergence criteria
   ConvergeType ConvergeType_;  //Quantities to check for convergence
   
public:
   RefineEqbmSolution(Restriction* const Restrict,Vector const& one,
                    double const& Converge,ConvergeType CnvrgTyp,int const& Echo=1);
   RefineEqbmSolution(Restriction* const Restrict,PerlInput const& Input,Vector const& one,
                    int const& Echo=1);
   RefineEqbmSolution(Restriction* const Restrict,PerlInput const& Input,int const& Echo);
   ~RefineEqbmSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const {return SolutionFound_;}
   virtual int FindNextSolution();
   virtual void FindCriticalPoint(Lattice* const Lat,int& TotalNumCPCrossings,
                                  PerlInput const& Input,int const& Width,ostream& out) {}
   virtual char const* const Type() const {return "RefineEqbmSolution";}
};

#endif
