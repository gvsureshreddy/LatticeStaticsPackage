#ifndef __NiTi9TPPLat
#define __NiTi9TPPLat

#include "GenericLat.h"
#include "UnitCellIterator.h"
#include "PPSum.h"
#include "RadiiMorse.h"
#include <CMatrix.h>

class NiTi9TPPLat : public GenericLat
{
private:
   const static int DIM3 = 3;
   const static int INTERNAL_ATOMS = 2;
   const static int DOFS = 9;

   double RefLen_;
   unsigned InfluanceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U23 V11 V12 V13]
   Vector DOF_;
   static const double LatticeBasis[DIM3][DIM3];
   Matrix RefLattice_;
   double ShearMod_;
   double Pressure_;
   Vector BodyForce_[INTERNAL_ATOMS];
   double AtomicMass_[INTERNAL_ATOMS];

   PPSum LatSum_;

   UnitCellIterator UCIter_;
   int GridSize_;

   //Pair Potential data
   RadiiMorse Potential_[INTERNAL_ATOMS][INTERNAL_ATOMS];

   // Misc
   double ConvexityDX_;
   static const double Alt[DIM3][DIM3][DIM3];
   Matrix A_;


   Matrix stress(PairPotentials::TDeriv dt=PairPotentials::T0);
   Matrix stiffness(int moduliflag=0,PairPotentials::TDeriv dt=PairPotentials::T0);
   Matrix CondensedModuli();

public:
   double RefLen() {return RefLen_;}

   // Virtual Functions required by GenericLat
   Vector DOF() {return DOF_;}
   void SetDOF(const Vector &dof) {DOF_ = dof; LatSum_.Recalc();}
   Matrix StressDT() {return stress(PairPotentials::DT);}
   Matrix StiffnessDT() {return stiffness(1,PairPotentials::DT);}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp;}

   // Virtual Functions required by Lattice
   virtual double Energy();
   virtual Matrix Stress() {return stress();}
   virtual Matrix Stiffness() {return stiffness();}
   virtual Matrix Moduli() {return stiffness(1);}
   virtual Matrix E3();
   virtual Matrix E4();
   virtual void Print(ostream &out,PrintDetail flag);
   
   // Functions provided by NiTi9TPPLat
   NiTi9TPPLat(char *datafile);
   ~NiTi9TPPLat() {}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) {return BodyForce_[i]; }
   double Pressure() const {return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const {return ShearMod_;}
   CMatrix DynamicalStiffness(Vector &Y);
   int BlochWave(Vector &Y);
   friend ostream &operator<<(ostream &out,NiTi9TPPLat &A);

private:
   double PI(double *Dx,double *DX,int r,int s);
   double PSI(double *DX,int r,int s,int t,int u);
   double OMEGA(double *Dx,int p,int q,int i, int j);
   double SIGMA(int p,int q,int i,int j,int k,int l);
   double GAMMA(double *Dx,double *DX,int p,int q,int i,int j,int k,int l);
   double THETA(double *DX,int p,int q,int i,int j,int k,int l,int m, int n);
   double XI(int p,int q,int i,int j,int k,int l,int m,int n);
   double LAMDA(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b);
   
   inline int INDU(int i,int j);
   inline int INDV(int i,int j);
   inline int INDUU(int k,int l,int m,int n);
   inline int INDVV(int k,int l,int m,int n);
   inline int INDUV(int i,int j,int m,int n);
   inline int INDVU(int m,int n,int i,int j);
   inline double DELTA(int s,int p,int q) {return Del(s,q) - Del(s,p);}
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
