#ifndef RSE__Lattice
#define RSE__Lattice

#include <Matrix.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>

#include "PerlInput.h"

using namespace std;

class Lattice
{
public:
   int Echo_;
   int dbg_;

   enum LoadType {Temperature, Load};
   enum StateType {LHS, RHS, INTERMED};
   LoadType LoadParameter_;
   LoadType const& LoadParameter() const
   {
      return LoadParameter_;
   }
   int UseEigenValTFs_;  // 1 use eigenvalues of E2 as TFs, 0 don't
   int NumExtraTFs_;     // Number of TFs provided by inherited class
   int LSKAnalysis_;     // None,SecondOrder,ThridOrder - do LSK to appropriate order
   int FullPrint_;       // 1 Print out full vectors and matricies, 0 don't
   int GuessModes_;      // guess new RestrictToTranslatedSubSpace mode dofs
   string UseExtension_; // define an extension to use for input files

   Lattice(PerlInput const& Input, int const& Echo);
   virtual ~Lattice()
   {
   }

   virtual Vector const& DOF() const = 0;
   virtual void SetDOF(Vector const& dof) = 0;
   virtual double Entropy() const = 0;
   virtual double HeatCapacity() const = 0;
   virtual Vector const& StressDT() const = 0;
   virtual Matrix const& StiffnessDT() const = 0;
   virtual double Temp() const = 0;
   virtual void SetTemp(double const& temp) = 0;
   virtual Vector const& StressDL() const = 0;
   virtual Matrix const& StiffnessDL() const = 0;
   virtual double Lambda() const = 0;
   virtual void SetLambda(double const& lambda) = 0;
   void SetLoadParameter(double const& load);

   virtual double E0() const = 0;
   virtual Vector const& E1() const = 0;
   virtual Vector const& E1DLoad() const = 0;
   virtual Matrix const& E2() const = 0;
   virtual Matrix const& E3() const = 0;
   virtual Matrix const& E4() const = 0;
   int const UseEigenValTFs() const
   {
      return UseEigenValTFs_;
   }

   int const NumTestFunctions() const
   {
      return (UseEigenValTFs_) ? (DOF().Dim() + NumExtraTFs_) : NumExtraTFs_;
   }

   virtual int TestFunctions(Vector& TF1, StateType const& State = LHS, Vector* const EV2 = 0)
   const;
   virtual void ExtraTestFunctions(Vector& TF) const
   {
   }

   virtual Matrix const& LHSEigVect()
   {
      return EigVectLHS_static;
   }

   virtual Matrix const& RelativeEigVects() const
   {
      return EigVect_static;
   }

   virtual void DispersionCurves(Vector const& K, int const& NoPTS, char const* const prefix,
                                 ostream& out) const
   {
   }

   virtual int BlochWave(Vector& K) const
   {
      return -1;
   }

   virtual void LongWavelengthModuli(double const& dk, int const& gridsize,
                                     char const* const prefix, ostream& out) const
   {
   }

   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1) = 0;
   virtual void SetGridSize(int const& Grid) = 0;
   virtual void NeighborDistances(int const& cutoff, ostream& out) const
   {
   }

   virtual int CriticalPointInfo(int* const CPCrossingNum, int const& TFIndex,
                                 Vector const& DrDt, int const& CPorBif,
                                 int const& NumZeroEigenVals, double const& Tolerance,
                                 int const& Width, PerlInput const& Input, ostream& out);
   void ConsistencyCheck(double const& ConsistencyEpsilon, int const& Width, ostream& out);
   virtual void DebugMode()
   {
   }

   virtual char const* const Type() const = 0;
   enum PrintDetail {PrintLong, PrintShort};
   enum PrintPathSolutionType {NotSolutionPt, RegularPt, TurningPt, BifurcationPt, ExtraTFPt};
   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt) = 0;
   friend ostream& operator<<(ostream& out, Lattice& L)
   {
      L.Print(out, PrintShort); return out;
   }

protected:
   // "static" member variables
   // TestFunctions
   mutable int test_flag_static;
   mutable Matrix Stiffness_1_static;
   mutable Matrix Stiffness_2_static;
   mutable Matrix Stiffness_3_static;
   mutable Matrix Stiffness_temp_static;
   mutable Matrix Stiffness_diagonalized_static;
   mutable Matrix EigVect_static;
   mutable Matrix EigVectRHS_static;
   mutable Matrix EigVectLHS_static;
   mutable Matrix EV1_static;
   mutable Matrix EV2_static;
   mutable Matrix EV3_static;
   mutable Vector ExTF1_static;
   mutable Vector ExTF2_static;
   mutable Vector ExTF3_static;
};

#endif
