#ifndef RSE__QHQMultiChainTPP
#define RSE__QHQMultiChainTPP

#include "PerlInput.h"
#include "Lattice.h"
#include "ChainIterator.h"
#include "SCLDChainSum.h"
#include "KnownPairPotentials.h"
#include <CMatrix.h>
#include <Vector.h>

using namespace std;

class QHQMultiChainTPP : public Lattice
{
private:
   const static int DIM1;

   int INTERNAL_ATOMS;
   int DOFS;

   double InfluenceDist_;
   double NTemp_;
   // DOF[i] = [F S0 S1 S2 S3 ...]
   Vector DOF_;
   int LagrangeCB_;
   Matrix RefLattice_;
   int Density_;
   double NormModulus_;
   double Tref_;
   double kB_;
   double h_;
   double PhiRef_;
   double EntropyRef_;
   double HeatCapacityRef_;
   enum LDeriv {L0, DL};
   double Lambda_;
   Vector* BodyForce_;
   double* SpeciesMass_;
   double* AtomicMass_;

   mutable SCLDChainSum SCLDChainSum_;
   mutable ChainIterator ChainIter_;

   double RefNTemp_;
   double RefLambda_;

   mutable ChainIterator ChainIterNew_;
   mutable int FreqCached;
   mutable int FirstConverged;
   mutable int FirstPrintLong;

   int GridSize_;

   // Pair Potential data
   int NumberofSpecies_;
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials*** SpeciesPotential_;
   PairPotentials*** Potential_;

   // Misc
   Vector* AtomPositions_;

   double energy(PairPotentials::TDeriv const& dt = PairPotentials::T0) const;
   Vector const& stress(PairPotentials::TDeriv const& dt = PairPotentials::T0, LDeriv const& dl = L0) const;
   Matrix const& stiffness(PairPotentials::TDeriv const& dt = PairPotentials::T0,
                           LDeriv const& dl = L0) const;

   double FreeEnergy(PairPotentials::TDeriv const& dt = PairPotentials::T0) const;
   Vector const& Fstress(PairPotentials::TDeriv const& dt = PairPotentials::T0, LDeriv const& dl = L0) const;
   Matrix const& Fstiffness(PairPotentials::TDeriv const& dt = PairPotentials::T0, LDeriv const& dl = L0) const;

   void ReferenceDispersionCurves(Vector const& K, int const& NoPTS, char const* const prefix,
                                  ostream& out) const;

   void ReferenceHarmonic() const;
   void ReferenceV4() const;
   void ReferencePseudoHarmonic() const;
   CMatrix const& ReferenceDynamicalStiffness(Vector const& K, PairPotentials::TDeriv const& dt = PairPotentials::T0,
                                              int const& DOFderiv = 0, int const& x = 0, int const& y = 0) const;
   Vector const& ThermalExpansion() const;

   int ReferenceBlochWave(Vector& K) const;
   // Needed for DispersionCurves()
   //
   // find next eigval in position two based on previous two values
   // stored in zero and one position.
   static void interpolate(Matrix* const EigVals, int const& zero, int const& one, int const& two);
   // compair function for qsort
   static int comp(void const* const a, void const* const b);
   static int abscomp(void const* const a, void const* const b);

public:
   Vector const& AtomPositions(int const& i) const
   {
      return AtomPositions_[i];
   }

   // Virtual Functions required by Lattice
   Vector const& DOF() const
   {
      return DOF_;
   }

   Matrix const& RefLattice() const
   {
      return RefLattice_;
   }
   void SetDOF(Vector const& dof)
   {
      DOF_ = dof;
      FreqCached = 0; SCLDChainSum_.Recalc();
   }
   // Entropy is NEGATIVE dE/dT
   double Entropy() const
   {
      return -FreeEnergy(PairPotentials::DT);
   }
   double HeatCapacity() const
   {
      return -NTemp_* FreeEnergy(PairPotentials::D2T);
   }
   Vector const& StressDT() const
   {
      return Fstress(PairPotentials::DT);
   }
   Matrix const& StiffnessDT() const
   {
      return Fstiffness(PairPotentials::DT);
   }
   void SetTemp(double const& Ntemp)
   {
      NTemp_ = Ntemp;
      FreqCached = 0; SCLDChainSum_.Recalc();
   }
   Vector const& StressDL() const
   {
      return Fstress(PairPotentials::T0, DL);
   }
   Matrix const& StiffnessDL() const
   {
      return Fstiffness(PairPotentials::T0, DL);
   }
   virtual void SetGridSize(int const& Grid)
   {
      GridSize_ = Grid; ChainIter_(GridSize_, 0, 0);
   }
   // ThermalExpansion
   mutable Vector ThermalExp_static;

   double Temp() const
   {
      return NTemp_;
   }

   virtual Vector const& E1DLoad() const
   {
      return (LoadParameter_ == Temperature) ? StressDT() : StressDL();
   }
   double Lambda() const
   {
      return Lambda_;
   }
   void SetLambda(double const& lambda)
   {
      Lambda_ = lambda;
   }

   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Matrix const& E2() const;
   Matrix const& CondensedModuli() const;

   virtual Matrix const& E3() const;
   virtual Matrix const& E4() const;
   virtual void DispersionCurves(Vector const& K, int const& NoPTS, char const* const prefix,
                                 ostream& out) const
   {
      ReferenceDispersionCurves(K, NoPTS, prefix, out);
   }
   virtual int BlochWave(Vector& K) const
   {
      return ReferenceBlochWave(K);
   }
   virtual void LongWavelengthModuli(double const& dk, int const& gridsize,
                                     char const* const prefix, ostream& out) const;
   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1);
   virtual void NeighborDistances(int const& cutoff, ostream& out) const;
   virtual char const* const Type() const
   {
      return "QHQMultiChainTPP";
   }
   virtual void DebugMode();
   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   // Functions provided by QHQMultiChainTPP
   QHQMultiChainTPP(PerlInput const& Input, int const& Echo = 1, int const& Width = 20,
                    int const& Debug = 0);
   ~QHQMultiChainTPP();

   double InfluenceDist() const
   {
      return InfluenceDist_;
   }
   void SetInfluenceDist(double const& InfluenceDist)
   {
      InfluenceDist_ = InfluenceDist;
   }

   inline double Del(int const& i, int const& j) const
   {
      return i == j;
   }
   Vector const& BodyForce(int const& i) const
   {
      return BodyForce_[i];
   }
   double NormModulus() const
   {
      return NormModulus_;
   }
   friend ostream& operator<<(ostream& out, QHQMultiChainTPP& A);

private:
   double PI(double const* const Dx, double const* const DX) const;
   double PSI(double const* const DX) const;
   double OMEGA(double const* const Dx, int const& p, int const& q, int const& i) const;
   double SIGMA(int const& p, int const& q, int const& i, int const& j) const;
   double GAMMA(double const* const Dx, double const* const DX, int const& p, int const& q,
                int const& i) const;
   double THETA(double const* const DX, int const& p, int const& q, int const& i) const;
   double XI(int const& p, int const& q, int const& i, int const& j) const;
   double LAMDA(int const& p, int const& q, int const& i, int const& j) const;

   inline double DELTA(int const& s, int const& p, int const& q) const
   {
      return Del(s, q) - Del(s, p);
   }
   int FindLatticeSpacing(int const& iter);
   void RefineEqbm(double const& Tol, int const& MaxItr, ostream* const out);

   // "static" member variables
   // E1
   mutable Vector Phi1_static;
   // stress
   mutable Vector stress_static;
   mutable Vector Fstress_static;
   // E2
   mutable Matrix Phi2_static;
   mutable Matrix FPhi2_static;
   // stiffness
   mutable Matrix stiff_static;
   // CondensedModuli
   mutable Matrix CM_static;

   mutable Vector TE_static;

   // E3
   mutable Matrix Phi3_static;
   // E4
   mutable Matrix Phi4_static;
   // ReferenceDynamicalStiffness
   mutable CMatrix Dk_static;

   // ReferenceBlochWave
   mutable CMatrix A_static;
   mutable CMatrix A_DOF_static;
   mutable CMatrix A_DOF2_static;
   mutable CMatrix A_DOFDOF_static;

   mutable Matrix* HEigVals_static;
   mutable Matrix** HEigValsDOF_static;
   mutable Matrix*** HEigValsDOFDOF_static;

   mutable CMatrix** HEigVec_static;
   mutable CMatrix*** HEigVecDOF_static;
   mutable CMatrix**** HEigVecDOFDOF_static;

   mutable Matrix V4_static;
   mutable Matrix* V4DOF_static;
   mutable Matrix** V4DOFDOF_static;

   mutable Matrix* EigVals_previous_static;

   mutable Vector* Z_static;

   mutable Matrix* EigVals_static;
   mutable Matrix* EigValsT_static;
   mutable Matrix* EigValsTT_static;
   mutable Matrix** EigValsDOF_static;
   mutable Matrix** EigValsTDOF_static;
   mutable Matrix*** EigValsDOFDOF_static;

   mutable Matrix InitialEigVals_static;
   mutable int InitialEigVals_available;
   mutable Matrix InitialHEigVals_static;
   mutable int InitialHEigVals_available;

   // Print
   mutable Vector str_static;
   mutable Matrix pstiff_static;
   mutable Vector TestFunctVals_static;
};

#endif

