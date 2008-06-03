#ifndef __MultiLatticeTPP
#define __MultiLatticeTPP

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
   CBKinematics *CBK_;
   int KillRotations_;
   int KillTranslations_;
   Vector KillOneRotation_;
   int Density_;
   double NormModulus_;
   double Tref_;
   double PhiRef_;
   double EntropyRef_;
   double HeatCapacityRef_;
   enum LDeriv {L0,DL};
   double Lambda_;
   Vector LoadingProportions_;
   double *EulerAng_;
   Matrix Rotation_;
   Matrix Loading_;
   Vector *BodyForce_;
   double *SpeciesMass_;
   double *AtomicMass_;
   
   PPSum LatSum_;
   
   UnitCellIterator UCIter_;
   int GridSize_;
   
   //Pair Potential data
   int NumberofSpecies_;
   int AtomSpecies_[100]; // Max number of atoms in unit cell. might need to be changed...
   PairPotentials ***SpeciesPotential_;
   PairPotentials ***Potential_;
   
   // Misc
   double ConvexityDX_;
   
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
   const Vector& AtomPositions(int i) {return CBK_->AtomPositions(i);}
   
   
   // Virtual Functions required by Lattice
   Vector DOF() {return CBK_->DOF();}
   void SetDOF(const Vector &dof) {CBK_->SetDOF(dof); LatSum_.Recalc();}
   // Entropy is NEGATIVE dE/dT
   double Entropy() {return -energy(PairPotentials::DT);}
   double HeatCapacity() {return -(NTemp_)*energy(PairPotentials::D2T);}
   Matrix StressDT() {return stress(PairPotentials::DT);}
   Matrix StiffnessDT() {return stiffness(PairPotentials::DT);}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp; LatSum_.Recalc();}
   Matrix StressDL() {return stress(PairPotentials::T0,DL);}
   Matrix StiffnessDL() {return stiffness(PairPotentials::T0,DL);}
   virtual Matrix E1DLoad() {return (LoadParameter_==Temperature)?StressDT():StressDL();}
   double Lambda() {return Lambda_;}
   void SetLambda(const double &lambda) {Lambda_ = lambda;}
   
   virtual double E0();
   virtual Matrix E1();
   virtual Matrix E2();
   virtual Matrix E3();
   virtual Matrix E4();
   virtual void DispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out)
   {ReferenceDispersionCurves(K,NoPTS,prefix,out);}
   virtual int BlochWave(Vector &K)
   {return ReferenceBlochWave(K);}
   virtual void LongWavelengthModuli(double dk,int gridsize,const char *prefix,
                                     ostream &out);
   virtual void SetParameters(double *Vals,int ResetRef = 1);
   virtual void SetGridSize(int Grid) {GridSize_=Grid; UCIter_(GridSize_);}
   virtual void NeighborDistances(int cutoff,ostream &out);
   virtual void DebugMode();
   virtual void Print(ostream &out,PrintDetail flag);
   
   void PrintCurrentCrystalParamaters(ostream &out);
   
   // Functions provided by MultiLatticeTPP
   MultiLatticeTPP(PerlInput &Input,int Echo=1,int Width=20,int Debug=0);
   ~MultiLatticeTPP();
   double InfluenceDist() {return InfluenceDist_;}
   void SetInfluenceDist(double InfluenceDist) {InfluenceDist_=InfluenceDist;}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) {return BodyForce_[i]; }
   double NormModulus() const {return NormModulus_;}
   const Matrix& RefLattice() {return CBK_->RefLattice();}
   Matrix CondensedModuli();
   Matrix ThermalExpansion();
   friend ostream &operator<<(ostream &out,MultiLatticeTPP &A);
   
private:
   int FindLatticeSpacing(int iter);
   void RefineEqbm(double Tol,int MaxItr,ostream *out);

   // member variables used to avoid repeated memory allocation/deallocation
   // and thus, improve performance.

   // E0
   double Phi0, Tsq[3], Rsq[3];
   // E1
   Matrix ME1;
   double T[3], R[3];
   // stress
   Matrix S;
   // E2
   Matrix ME2;
   // stiffness
   Matrix Phi2;
   // E3
   Matrix Phi3;
   // E4
   Matrix Phi4;
   // ReferenceDynamicalStiffness
   CMatrix Dk;
   // ReferenceBlochWave
   CMatrix A;
   Matrix EigVals, InverseLat;
   Vector Z;
   // Print
   Matrix str, stiff, CondEV, TE, CondModuli;
   Vector TestFunctVals, K;
};

#endif
