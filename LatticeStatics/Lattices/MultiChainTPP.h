#ifndef __MultiChainTPP
#define __MultiChainTPP

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
   const static int DIM1;
   
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
   
   ChainSum ChainSum_;
   
   ChainIterator ChainIter_;
   int GridSize_;
   
   //Pair Potential data
   int NumberofSpecies_;
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials ***SpeciesPotential_;
   PairPotentials ***Potential_;
   
   // Misc
   Vector *AtomPositions_;
   
   double energy(PairPotentials::TDeriv dt=PairPotentials::T0);
   Matrix stress(PairPotentials::TDeriv dt=PairPotentials::T0,LDeriv dl=L0);
   Matrix stiffness(PairPotentials::TDeriv dt=PairPotentials::T0,
                    LDeriv dl=L0);
   
   void ReferenceDispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out);
   int ReferenceBlochWave(Vector &K);
   CMatrix ReferenceDynamicalStiffness(Vector &K);
   // Needed for DispersionCurves()
   //
   // find next eigval in position two based on previous two values
   // stored in zero and one position.
   void interpolate(Matrix *EigVals,int zero,int one,int two);
   // compair function for qsort
   static int comp(const void *a,const void *b);
   static int abscomp(const void *a,const void *b);
   
public:
   Vector AtomPositions(int i) {return AtomPositions_[i];}
   
   
   // Virtual Functions required by Lattice
   Vector DOF() {return DOF_;}
   void SetDOF(const Vector &dof) {DOF_ = dof; ChainSum_.Recalc();}
   // Entropy is NEGATIVE dE/dT
   double Entropy() {return -energy(PairPotentials::DT);}
   double HeatCapacity() {return -(NTemp_*Tref_)*energy(PairPotentials::D2T);}
   Matrix StressDT() {return stress(PairPotentials::DT);}
   Matrix StiffnessDT() {return stiffness(PairPotentials::DT);}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp; ChainSum_.Recalc();}
   Matrix StressDL() {return stress(PairPotentials::T0,DL);}
   Matrix StiffnessDL() {return stiffness(PairPotentials::T0,DL);}
   virtual Matrix E1DLoad() {return (LoadParameter_==Temperature)?StressDT():StressDL();}
   double Lambda() {return Lambda_;}
   void SetLambda(const double &lambda) {Lambda_ = lambda;}
   
   virtual double E0() {return energy();}
   virtual Matrix E1() {return stress();}
   virtual Matrix E2() {return stiffness();}
   Matrix CondensedModuli();
   
   virtual Matrix E3();
   virtual Matrix E4();
   virtual void DispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out)
   {ReferenceDispersionCurves(K,NoPTS,prefix,out);}
   virtual int BlochWave(Vector &K)
   {return ReferenceBlochWave(K);}
   virtual void LongWavelengthModuli(double dk,int gridsize,const char *prefix,
                                     ostream &out);
   virtual void SetParameters(double *Vals,int ResetRef = 1);
   virtual void SetGridSize(int Grid) {GridSize_=Grid; ChainIter_(GridSize_);}
   virtual void NeighborDistances(int cutoff,ostream &out);
   virtual void DebugMode();
   virtual void Print(ostream &out,PrintDetail flag);
   
   // Functions provided by MultiChainTPP
   MultiChainTPP(PerlInput &Input,int Echo=1,int Width=20,int Debug=0);
   ~MultiChainTPP();
   
   double InfluenceDist() {return InfluenceDist_;}
   void SetInfluenceDist(double InfluenceDist) {InfluenceDist_=InfluenceDist;}
   
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) {return BodyForce_[i]; }
   double NormModulus() const {return NormModulus_;}
   friend ostream &operator<<(ostream &out,MultiChainTPP &A);
   
private:
   double PI(double *Dx,double *DX);
   double PSI(double *DX);
   double OMEGA(double *Dx,int p,int q,int i);
   double SIGMA(int p,int q,int i,int j);
   double GAMMA(double *Dx,double *DX,int p,int q,int i);
   double THETA(double *DX,int p,int q,int i);
   double XI(int p,int q,int i,int j);
   double LAMDA(int p,int q,int i,int j);
   
   inline double DELTA(int s,int p,int q) {return Del(s,q) - Del(s,p);}
   int FindLatticeSpacing(int iter);
   void RefineEqbm(double Tol,int MaxItr,ostream *out);
   
   // "static" member variables
   //stress
   Matrix Phi1_static;
   // stiffness
   Matrix Phi2_static;
   // E3
   Matrix Phi3_static;
   // E4
   Matrix Phi4_static;
   // ReferenceDynamicalStiffness
   CMatrix Dk_static;
   // ReferenceBlochWave
   CMatrix A_static;
   Matrix EigVals_static;
   // Print
   Matrix str_static;
   Matrix stiff_static;
   Vector TestFunctVals_static;
};

#endif
