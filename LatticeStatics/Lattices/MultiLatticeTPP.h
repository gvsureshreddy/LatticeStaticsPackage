#ifndef __MultiLatticeTPP
#define __MultiLatticeTPP

#include "Lattice.h"
#include "UnitCellIterator.h"
#include "CBKinematics.h"
#include "LagrangeCB.h"
#include "EulerCB.h"
#include "PPSum.h"
#include "KnownPairPotentials.h"
#include <CMatrix.h>

using namespace std;

class MultiLatticeTPP : public Lattice
{
private:
   const static int DIM3 = 3;
   
   int INTERNAL_ATOMS;
   int DOFS;

   unsigned InfluenceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U23 V11 V12 V13 V21... ...]
   Vector DOF_;
   Matrix RefLattice_;
   double NormModulus_;
   double Tref_;
   double PhiRef_;
   double EntropyRef_;
   double HeatCapacityRef_;
   enum LDeriv {L0,DL};
   double Lambda_;
   Vector LoadingProportions_;
   double EulerAng_[DIM3];
   Matrix Rotation_;
   Matrix Loading_;
   Vector *BodyForce_;
   double *AtomicMass_;

   CBKinematics *CBK_;
   PPSum LatSum_;

   UnitCellIterator UCIter_;
   int GridSize_;

   //Pair Potential data
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials ***Potential_;

   // Misc
   double ConvexityDX_;
   Vector *AtomPositions_;

   double energy(PairPotentials::TDeriv dt=PairPotentials::T0);
   Matrix stress(PairPotentials::TDeriv dt=PairPotentials::T0,LDeriv dl=L0);
   Matrix stiffness(PairPotentials::TDeriv dt=PairPotentials::T0,
		    LDeriv dl=L0);
   Matrix CondensedModuli();

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
   void SetDOF(const Vector &dof) {DOF_ = dof; LatSum_.Recalc();}
   // Entropy is NEGATIVE dE/dT
   double Entropy() {return -energy(PairPotentials::DT);}
   double HeatCapacity() {return -(NTemp_*Tref_)*energy(PairPotentials::D2T);}
   Matrix StressDT() {return stress(PairPotentials::DT);}
   Matrix StiffnessDT() {return stiffness(PairPotentials::DT);}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp; LatSum_.Recalc();}
   Matrix StressDL() {return stress(PairPotentials::T0,DL);}
   Matrix StiffnessDL() {return stiffness(PairPotentials::T0,DL);}
   double Lambda() {return Lambda_;}
   void SetLambda(const double &lambda) {Lambda_ = lambda;}

   virtual double Energy() {return energy();}
   virtual Matrix Stress() {return stress();}
   virtual Matrix Stiffness() {return stiffness();}

   virtual Matrix E3();
   virtual Matrix E4();
   virtual void DispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out)
   {ReferenceDispersionCurves(K,NoPTS,prefix,out);}
   virtual int BlochWave(Vector &K)
   {return ReferenceBlochWave(K);}
   virtual void LongWavelengthModuli(double dk,int gridsize,const char *prefix,
				     ostream &out);
   virtual void SetGridSize(int Grid) {GridSize_=Grid; UCIter_(GridSize_);}
   virtual void NeighborDistances(int cutoff,ostream &out);
   virtual void DebugMode();
   virtual void Print(ostream &out,PrintDetail flag);

   void PrintCurrentCrystalParamaters(ostream &out);
   
   // Functions provided by MultiLatticeTPP
   MultiLatticeTPP(char *datafile,const char *prefix,int Echo=1,int Width=20,int Debug=0);
   ~MultiLatticeTPP();
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) {return BodyForce_[i]; }
   double NormModulus() const {return NormModulus_;}
   friend ostream &operator<<(ostream &out,MultiLatticeTPP &A);

private:
   inline int INDU(int i,int j);
   inline int INDV(int i,int j);
   inline int INDUU(int k,int l,int m,int n);
   inline int INDVV(int k,int l,int m,int n);
   inline int INDUV(int i,int j,int m,int n);
   inline int INDVU(int m,int n,int i,int j);
   int FindLatticeSpacing(char *datafile,const char *prefix,int iter);
   void RefineEqbm(double Tol,int MaxItr,ostream *out);
   
};

#endif
