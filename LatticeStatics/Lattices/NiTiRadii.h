#ifndef __NiTiRadii
#define __NiTiRadii

#include "GenericLat.h"

#define DIM3 3
#define INTERNAL_ATOMS 4

class NiTiRadii : public GenericLat
{
private:
   double RefLen_;
   unsigned InfluanceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U12 S D A]
   Vector DOF_;
   Vector LatticeVec_[DIM3];
   double ShearMod_;
   double Pressure_;
   Vector BodyForce_[INTERNAL_ATOMS];

   //Pair Potential data
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   enum interaction {aa,bb,ab};
   double Tref_;
   double A0_aa, B0_aa, Alpha_aa, Rref_aa, Rtheta_aa, Tmelt_aa;
   double A0_bb, B0_bb, Alpha_bb, Rref_bb, Rtheta_bb, Tmelt_bb;
   double A0_ab, B0_ab, Alpha_ab, Rref_ab, Rtheta_ab, Tmelt_ab;

   // Misc
   double ConvexityDX_;
   static const double Alt[DIM3][DIM3][DIM3];

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
   virtual void CriticalPointInfo(int Width,ostream &out);
   
   // Functions provided by NiTiRadii
   NiTiRadii(char *datafile);
   ~NiTiRadii() {}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) { return BodyForce_[i]; }
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,NiTiRadii &A);
private:
   double PairPotential(interaction inter,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   inline double Beta(interaction inter,TDeriv dt=T0);
   inline double Rhat(interaction inter,TDeriv dt=T0);
   void GetLatticeVectorInfo(double *SX,double *DXPrimeS,double *DXPrimeD,
			     double *DXPrimeA,interaction &Inter,int p,int q);
   inline double PI(const Vector &Dx,const Vector &DX,int r,int s);
   inline double PSI(const Vector &DX,int r,int s,int t,int u);
   double pwr(const double &x,const unsigned y);
   inline int IND(int i,int j);
   inline int IND(int k,int l,int m,int n);
   Matrix Phi(unsigned moduliflag=0,YDeriv dy=Y0,TDeriv dt=T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
