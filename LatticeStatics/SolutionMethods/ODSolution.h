#ifndef RSE__ODSolution
#define RSE__ODSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "Restriction.h"

using namespace std;

class Lattice;

class ODSolution : public SolutionMethod
{
public:
   enum ConvergeType {Both, Force, Displacement};
private:
   Restriction* Restrict_;

   int Echo_;
   int SolutionFound_;
   int DimDOFS_;
   int Width_;
   double dt_;
   double t_;
   double Epsilon_;
   double ForceNorm_;
   double Potential_;
   double Potential0_;
   double Potential1_;
   Vector Tangent_;
   Vector X0_;
   Vector X_;
   Vector K1_;
   Vector K2_;
   Vector K3_;
   Vector K4_;
   Vector Delta_;
   Vector CurrentSolution_;

   double Converge_;            // Convergence criteria
   ConvergeType ConvergeType_;  // Quantities to check for convergence

public:
   ODSolution(Restriction* const Restrict, PerlInput const& Input, int const& Echo);
   ~ODSolution()
   {
   }

   // Functions required by SolutionMethod
   virtual Vector const& Solution();
   virtual Vector const& X0();
   virtual Vector const& X1();
   virtual int FindNextSolution(PerlInput const& Input, int const& Width, ostream& out);
   virtual int AllSolutionsFound() const;
   virtual double Time() const;
   virtual double dt() const;
   virtual double Energy() const;
   virtual double Energy0() const;
   virtual double Energy1() const;
   virtual void FindCriticalPoint(Lattice* const Lat, int* const TotalNumCPCrossings,
                                  PerlInput const& Input, int const& Width, ostream& out)
   {
   }
   virtual char const* const Type() const
   {
      return "ODSolution";
   }
};

#endif
