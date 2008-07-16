#ifndef RSE__MultiChainTPP
#define RSE__MultiChainTPP

#include "PerlInput.h"
#include "Lattice.h"
#include "ChainIterator.h"
#include "ChainSum.h"
#include "KnownPairPotentials.h"
#include <CMatrix.h>

using namespace std;

class MultiChainTPP : public Lattice
{
private:
   static int const DIM1;
   
   int INTERNAL_ATOMS;
   int DOFS;
   
   double InfluenceDist_;
   double NTemp_;
   // DOF[i] = [F S1 S2 S3 ...]
   Vector DOF_;
   int LagrangeCB_;
   Matrix RefLattice_;
   int Density_;
   double NormModulus_;
   double Tref_;
   double PhiRef_;
   double EntropyRef_;
   double HeatCapacityRef_;
   enum LDeriv {L0,DL};
   double Lambda_;
   Vector *BodyForce_;
   double *SpeciesMass_;
   double *AtomicMass_;
   
   mutable ChainSum ChainSum_;
   mutable ChainIterator ChainIter_;
   int GridSize_;
   
   //Pair Potential data
   int NumberofSpecies_;
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials ***SpeciesPotential_;
   PairPotentials ***Potential_;
   
   // Misc
   Vector *AtomPositions_;
   
   double energy(PairPotentials::TDeriv const& dt=PairPotentials::T0) const;
   Matrix const& stress(PairPotentials::TDeriv const& dt=PairPotentials::T0,LDeriv const& dl=L0)
      const;
   Matrix const& stiffness(PairPotentials::TDeriv const& dt=PairPotentials::T0,
                           LDeriv const& dl=L0) const;
   
   void ReferenceDispersionCurves(Vector const& K,int const& NoPTS,char const* const prefix,
                                  ostream& out) const;
   int ReferenceBlochWave(Vector& K) const;
   CMatrix const& ReferenceDynamicalStiffness(Vector const& K) const;
   // Needed for DispersionCurves()
   //
   // find next eigval in position two based on previous two values
   // stored in zero and one position.
   static void interpolate(Matrix* const EigVals,int const& zero,int const& one,int const& two);
   // compair function for qsort
   static int comp(void const* const a,void const* const b);
   static int abscomp(void const* const a,void const* const b);
   
public:
   Vector const& AtomPositions(int const& i) const {return AtomPositions_[i];}
   
   
   // Virtual Functions required by Lattice
   Vector const& DOF() const {return DOF_;}
   void SetDOF(Vector const& dof) {DOF_ = dof; ChainSum_.Recalc();}
   // Entropy is NEGATIVE dE/dT
   double Entropy() const {return -energy(PairPotentials::DT);}
   double HeatCapacity() const {return -NTemp_*energy(PairPotentials::D2T);}
   Matrix const& StressDT() const {return stress(PairPotentials::DT);}
   Matrix const& StiffnessDT() const {return stiffness(PairPotentials::DT);}
   double Temp() const {return NTemp_;}
   void SetTemp(double const& Ntemp) {NTemp_ = Ntemp; ChainSum_.Recalc();}
   Matrix const& StressDL() const {return stress(PairPotentials::T0,DL);}
   Matrix const& StiffnessDL() const {return stiffness(PairPotentials::T0,DL);}
   virtual Matrix const& E1DLoad() const
   {return (LoadParameter_==Temperature)?StressDT():StressDL();}
   double Lambda() const {return Lambda_;}
   void SetLambda(double const& lambda) {Lambda_ = lambda;}
   
   virtual double E0() const {return energy();}
   virtual Matrix const& E1() const {return stress();}
   virtual Matrix const& E2() const {return stiffness();}
   Matrix const& CondensedModuli() const;
   
   virtual Matrix const& E3() const;
   virtual Matrix const& E4() const;
   virtual void DispersionCurves(Vector const& K,int const& NoPTS,char const* const prefix,
                                 ostream& out) const
   {ReferenceDispersionCurves(K,NoPTS,prefix,out);}
   virtual int BlochWave(Vector& K) const {return ReferenceBlochWave(K);}
   virtual void LongWavelengthModuli(double const& dk,int const& gridsize,
                                     char const* const prefix,ostream& out) const;
   virtual void SetParameters(double const* const Vals,int const& ResetRef = 1);
   virtual void SetGridSize(int const& Grid) {GridSize_=Grid; ChainIter_(GridSize_);}
   virtual void NeighborDistances(int const& cutoff,ostream& out) const;
   virtual void DebugMode();
   virtual void Print(ostream& out,PrintDetail const& flag);
   
   // Functions provided by MultiChainTPP
   MultiChainTPP(PerlInput const& Input,int const& Echo=1,int const& Width=20,
                 int const& Debug=0);
   ~MultiChainTPP();
   
   double const& InfluenceDist() const {return InfluenceDist_;}
   void SetInfluenceDist(double const& InfluenceDist) {InfluenceDist_=InfluenceDist;}
   
   inline double Del(int const& i,int const& j) const {return i==j;}
   Vector const& BodyForce(int const& i) const {return BodyForce_[i];}
   double const& NormModulus() const {return NormModulus_;}
   friend ostream& operator<<(ostream& out,MultiChainTPP& A);
   
private:
   double PI(double const* const Dx,double const* const DX) const;
   double PSI(double const* const DX) const;
   double OMEGA(double const* const Dx,int const& p,int const& q,int const& i) const;
   double SIGMA(int const& p,int const& q,int const& i,int const& j) const;
   double GAMMA(double const* const Dx,double const* const DX,int const& p,int const& q,
                int const& i) const;
   double THETA(double const* const DX,int const& p,int const& q,int const& i) const;
   double XI(int const& p,int const& q,int const& i,int const& j) const;
   double LAMDA(int const& p,int const& q,int const& i,int const& j) const;
   
   inline double DELTA(int const& s,int const& p,int const& q) const
   {return Del(s,q) - Del(s,p);}
   int FindLatticeSpacing(int const& iter);
   void RefineEqbm(double const& Tol,int const& MaxItr,ostream* const out);
   
   // "static" member variables
   //stress
   mutable Matrix Phi1_static;
   // stiffness
   mutable Matrix Phi2_static;
   // CondenssedModuli
   mutable Matrix CM_static;
   // E3
   mutable Matrix Phi3_static;
   // E4
   mutable Matrix Phi4_static;
   // ReferenceDynamicalStiffness
   mutable CMatrix Dk_static;
   // ReferenceBlochWave
   mutable CMatrix A_static;
   mutable Matrix EigVals_static;
   // Print
   mutable Matrix str_static;
   mutable Matrix stiff_static;
   mutable Vector TestFunctVals_static;
};

#endif
