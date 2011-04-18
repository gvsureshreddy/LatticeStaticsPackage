#ifndef RSE__MultiLatticeTPP
#define RSE__MultiLatticeTPP

#include "PerlInput.h"
#include "Lattice.h"
#include "UnitCellIterator.h"
#include "CBKinematics.h"
#include "LagrangeCB.h"
#include "MixedCB.h"
#include "EulerCB.h"
#include "SymLagrangeWTransCB.h"
#include "SymLagrangeCB.h"
#include "PPSum.h"
#include "KnownPairPotentials.h"
#include <CMatrix.h>

using namespace std;

class MultiLatticeTPP : public Lattice
{
private:
   const static int DIM3;

   int InternalAtoms_;

   double InfluenceDist_;
   double NTemp_;
   // DOF[i] = [F00 F01 F02 F10 F11 F12 F20 F21 F22 S00 S01 S02 S11 S12 S13 S21... ...]
   // if using a FwithTransMapping CBkinematics
   // or
   // DOF[i] = [U11 U22 U33 U12 U13 U23 S11 S12 S13 S21... ...]
   // if using a UwithoutTransMapping CBkinematics
   CBKinematics* CBK_;
   int KillRotations_;
   int KillTranslations_;
   Vector KillOneRotation_;
   int Density_;
   double REFTemp_;
   double REFLambda_;
   double NormModulus_;
   double Tref_;
   double PhiRef_;
   double EntropyRef_;
   double HeatCapacityRef_;
   enum LDeriv {L0, DL};
   double Lambda_;
   Vector LoadingProportions_;
   double* EulerAng_;
   Matrix Rotation_;
   Matrix Loading_;
   Vector* BodyForce_;
   double* SpeciesMass_;
   double* AtomicMass_;

   mutable PPSum LatSum_;
   mutable UnitCellIterator UCIter_;
   int GridSize_;
   mutable Vector ExtraTestFunctions_;


   // Pair Potential data
   int NumberofSpecies_;
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials*** SpeciesPotential_;
   PairPotentials*** Potential_;

   // Misc
   double ConvexityDX_;

   double energy(PairPotentials::TDeriv const& dt = PairPotentials::T0) const;
   Vector const& stress(PairPotentials::TDeriv const& dt = PairPotentials::T0, LDeriv const& dl = L0)
   const;
   Matrix const& stiffness(PairPotentials::TDeriv const& dt = PairPotentials::T0,
                           LDeriv const& dl = L0) const;

   void ReferenceDispersionCurves(Vector const& K, int const& NoPTS, char const* const prefix,
                                  ostream& out) const;
   int ReferenceBlochWave(Vector& K) const;
   CMatrix const& ReferenceDynamicalStiffness(Vector const& K) const;
   // Needed for DispersionCurves()
   //
   // find next eigval in position two based on previous two values
   // stored in zero and one position.
   static void interpolate(Matrix* const EigVals, int const& zero, int const& one,
                           int const& two);
   // compair function for qsort
   static int comp(void const* const a, void const* const b);
   static int abscomp(void const* const a, void const* const b);

public:
   Vector const& AtomPositions(int const& i) const {return CBK_->AtomPositions(i);}


   // Virtual Functions required by Lattice
   Vector const& DOF() const {return CBK_->DOF();}
   void SetDOF(Vector const& dof) {CBK_->SetDOF(dof); LatSum_.Recalc();}
   // Entropy is NEGATIVE dE/dT
   double Entropy() const {return -energy(PairPotentials::DT);}
   double HeatCapacity() const {return -(NTemp_) * energy(PairPotentials::D2T);}
   Vector const& StressDT() const {return stress(PairPotentials::DT);}
   Matrix const& StiffnessDT() const {return stiffness(PairPotentials::DT);}
   double Temp() const {return NTemp_;}
   void SetTemp(double const& Ntemp) {NTemp_ = Ntemp; LatSum_.Recalc();}
   Vector const& StressDL() const {return stress(PairPotentials::T0, DL);}
   Matrix const& StiffnessDL() const {return stiffness(PairPotentials::T0, DL);}
   virtual Vector const& E1DLoad() const
   {return (LoadParameter_ == Temperature) ? StressDT() : StressDL();}
   double Lambda() const {return Lambda_;}
   double ConjugateToLambda() const;
   void SetLambda(double const& lambda) {Lambda_ = lambda;}

   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Matrix const& E2() const;
   virtual Matrix const& E3() const;
   virtual Matrix const& E4() const;
   virtual void ExtraTestFunctions(Vector& TF) const;
   virtual void DispersionCurves(Vector const& K, int const& NoPTS, char const* const prefix,
                                 ostream& out) const
   {ReferenceDispersionCurves(K, NoPTS, prefix, out);}
   virtual int BlochWave(Vector& K) const
   {return ReferenceBlochWave(K);}
   virtual void LongWavelengthModuli(double const& dk, int const& gridsize,
                                     char const* const prefix, ostream& out) const;
   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1);
   virtual void SetGridSize(int const& Grid) {GridSize_ = Grid; UCIter_(GridSize_);}
   void RefineEqbm(double const& Tol, int const& MaxItr, ostream* const out);
   virtual void NeighborDistances(int const& cutoff, ostream& out) const;
   virtual char const* const Type() const {return "MultiLatticeTPP";}
   virtual void DebugMode();
   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   void PrintCurrentCrystalParamaters(ostream& out) const;

   // Functions provided by MultiLatticeTPP
   MultiLatticeTPP(PerlInput const& Input, int const& Echo = 1, int const& Width = 20,
                   int const& Debug = 0);
   ~MultiLatticeTPP();
   double InfluenceDist() const {return InfluenceDist_;}
   void SetInfluenceDist(double const& InfluenceDist) {InfluenceDist_ = InfluenceDist;}
   inline double Del(int const& i, int const& j) const {return i == j;}
   Vector const& BodyForce(int const& i) const {return BodyForce_[i]; }
   double NormModulus() const {return NormModulus_;}
   Matrix const& RefLattice() const {return CBK_->RefLattice();}
   Matrix const& CondensedModuli() const;
   Vector const& ThermalExpansion() const;
   friend ostream& operator<<(ostream& out, MultiLatticeTPP& A);

private:
   int FindLatticeSpacing(int const& iter);

   // member variables used to avoid repeated memory allocation/deallocation
   // and thus, improve performance.

   // E0
   mutable double Phi0_static;
   mutable double Tsq_static[3];
   mutable double Rsq_static[3];
   // E1
   mutable Vector ME1_static;
   mutable double T_static[3];
   mutable double R_static[3];
   // stress
   mutable Vector S_static;
   // E2
   mutable Matrix ME2_static;
   // stiffness
   mutable Matrix Phi2_static;
   // CondensedModuli
   mutable Matrix CM_static;
   // ThermalExpansion
   mutable Vector ThermalExp_static;
   // E3
   mutable Matrix Phi3_static;
   // E4
   mutable Matrix Phi4_static;
   // ReferenceDynamicalStiffness
   mutable CMatrix Dk_static;
   // ReferenceBlochWave
   mutable CMatrix A_static;
   mutable Matrix EigVals_static;
   mutable Matrix InverseLat_static;
   mutable Vector Z_static;
   // Print
   mutable Vector str_static;
   mutable Matrix stiff_static;
   mutable Matrix CondEV_static;
   mutable Vector TE_static;
   mutable Matrix CondModuli_static;
   mutable Vector TestFunctVals_static;
   mutable Vector K_static;
};

#endif

