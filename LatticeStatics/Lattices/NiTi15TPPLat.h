#ifndef __NiTi15TPPLat
#define __NiTi15TPPLat

#include "GenericLat.h"
#include "RadiiMorse.h"

#define DIM3 3
#define INTERNAL_ATOMS 4
#define DOFS 15

class NiTi15TPPLat : public GenericLat
{
private:
   double RefLen_;
   unsigned InfluanceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U12 V11 V12 V13 V21 V22 V23 V31 V32 V33]
   Vector DOF_;
   Matrix LatticeVec_;
   double ShearMod_;
   double Pressure_;
   Vector BodyForce_[INTERNAL_ATOMS];

   //Pair Potential data
   enum interaction {aa,bb,ab,NOINTERACTIONS};
   RadiiMorse Potential_[NOINTERACTIONS];

   // Misc
   double ConvexityDX_;
   static const double Alt[DIM3][DIM3][DIM3];
   static const double A[INTERNAL_ATOMS][DIM3];
   static const interaction INTER[INTERNAL_ATOMS][INTERNAL_ATOMS];

public:
   Vector DOF() {return DOF_;}
   Matrix StressDT();
   void SetDOF(const Vector &dof) { DOF_ = dof;}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp;}
   double RefLen() {return RefLen_;}

   // Virtual Functions required by Lattice
   virtual double Energy();
   virtual Matrix Stress();
   virtual Matrix Stiffness();
   virtual Matrix Moduli();
   virtual int StiffnessNulity(double *Min=NULL);
   virtual void Print(ostream &out,PrintDetail flag);
   virtual void CriticalPointInfo(double Tolerance,int Width,ostream &out);
   
   // Functions provided by NiTi15TPPLat
   NiTi15TPPLat(char *datafile);
   ~NiTi15TPPLat() {}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) { return BodyForce_[i]; }
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,NiTi15TPPLat &A);
private:
   inline double PI(const Vector &Dx,const Vector &DX,int r,int s);
   inline double PSI(const Vector &DX,int r,int s,int t,int u);
   inline double OMEGA(const Vector &Dx,int p,int q,int i, int j);
   inline double SIGMA(int p,int q,int i,int j,int k,int l);
   inline double GAMMA(const Vector &Dx,const Vector &DX,int p,int q,
		       int i,int j,int k,int l);
   inline double THETA(const Vector &DX,int p,int q,int i,int j,int k,int l,
		       int m, int n);
   inline double XI(int p,int q,int i,int j,int k,int l,int m,int n);
   inline double LAMDA(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b);
   
   double pwr(const double &x,const unsigned y);
   inline int INDU(int i,int j);
   inline int INDV(int i,int j);
   inline int INDUU(int k,int l,int m,int n);
   inline int INDVV(int k,int l,int m,int n);
   inline int INDUV(int i,int j,int m,int n);
   inline int INDVU(int m,int n,int i,int j);
   inline double DELTA(int s,int p,int q) {return Del(s,q) - Del(s,p);}
   Matrix Phi(unsigned moduliflag=0,PairPotentials::YDeriv dy=PairPotentials::Y0,
	      PairPotentials::TDeriv dt=PairPotentials::T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
