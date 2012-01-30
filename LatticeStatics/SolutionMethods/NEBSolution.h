#ifndef RSE__NEBSolution
#define RSE__NEBSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "Restriction.h"

using namespace std;

class Lattice;

class NEBSolution : public SolutionMethod
{
public:
   enum ConvergeType {Both, Force, Displacement};
private:
   Restriction* Restrict_;

   int Echo_;
   int SolutionFound_;
   int NumSolutions_;
   int NumReplicas_;
   int NumInterStates_;
   int DimDOFS_;
   int DimReplicas_;
   int Width_;
   int VSeed_;
   double VScaling_;
   double SpringK_;
   double MaxElement_;
   double MinElement_;
   double Deltat_;
   double tFinal_;
   double QDMass_;
   double InitV_;
   double TotalReplicaForce_;
   Vector EnergyBarrier_;
   Vector InitState_;
   Vector FinalState_;
   Vector DOF_;
   Vector InitDOF_;
   Vector FinalDOF_;
   Vector FiTot_;
   Vector iTangent_;
   Vector TempState_;
   Matrix QuenchedState_;
   Matrix InterStates_;
   Matrix StateMatrix_;
   Matrix Replicas_;
   Matrix rReplica1_;
   Matrix rReplica2_;
   Matrix vReplica1_;
   Matrix vReplica2_;

   double Converge_;            // Convergence criteria
   ConvergeType ConvergeType_;  // Quantities to check for convergence

public:
   NEBSolution(Restriction* const Restrict, PerlInput const& Input, int const& Echo);
   ~NEBSolution()
   {
   }

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const
   {
      return SolutionFound_ >= NumSolutions_;
   }
   virtual int FindNextSolution(PerlInput const& Input, int const& Width, ostream& out);
   virtual Matrix const& StateMatrix()
   {
      return StateMatrix_;
   }
   virtual Vector const& EnergyBarrier()
   {
      return EnergyBarrier_;
   }
   virtual Vector const& RefineState(Vector const& CurrentState);
   virtual Vector const& NEBForce(Vector const& iMinusReplica, Vector const& iPlusReplica, Vector const& iReplica);
   virtual Vector const& NEBTangent(Vector const& iMinusReplica, Vector const& iPlusReplica, Vector const& iReplica);
   virtual Matrix const& NEBQuenchedDynamics(Vector const& InitDOF, Vector const& FinalDOF, Matrix const& Replicas_, double const& Deltat, double const& tFinal, double const& M);
   virtual Matrix const& NEBQDStep(Vector const& InitDOF, Vector const& FinalDOF, double const& Deltat, double const& M);
   virtual void FindCriticalPoint(Lattice* const Lat, int* const TotalNumCPCrossings,
                                  PerlInput const& Input, int const& Width, ostream& out)
   {
   }
   virtual char const* const Type() const
   {
      return "NEBSolution";
   }
};

#endif
