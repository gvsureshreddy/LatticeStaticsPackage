#ifndef __NiTiShuffleTPPLat
#define __NiTiShuffleTPPLat

#include "Lattice.h"

#define DIM3 3

class NiTiShuffleTPPLat : public Lattice
{
private:
   double RefLen_;
   unsigned InfluanceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U12 S D]
   Vector DOF_;
   Vector LatticeVec_[DIM3];
   double ShearMod_;
   double Pressure_;
   Vector BodyForce_[4];

   //Pair Potential data
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   enum TDeriv {T0,DT};
   enum interaction {aa,bb,ab};
   double Tref_;
   double A0_aa, B0_aa, Alpha_aa, Rref_aa, Tmelt_aa;
   double A0_bb, B0_bb, Alpha_bb, Rref_bb, Tmelt_bb;
   double A0_ab, B0_ab, Alpha_ab, Rref_ab, Tmelt_ab;

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
   
   // Functions provided by NiTiShuffleTPPLat
   NiTiShuffleTPPLat(char *datafile);
   ~NiTiShuffleTPPLat() {}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) { return BodyForce_[i]; }
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,NiTiShuffleTPPLat &A);
private:
   double PairPotential(interaction inter,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   inline double Beta(interaction inter,TDeriv dt=T0);
   inline double Rhat(interaction inter,TDeriv dt=T0);
   void GetLatticeVectorInfo(double *SX,double *DXPrimeS,double *DXPrimeD
			     ,interaction &Inter,int p,int q);
   inline double PI(const Vector &Dx,const Vector &DX,int r,int s);
   inline double PSI(const Vector &DX,int r,int s,int t,int u);
   double pwr(const double &x,const unsigned y);
   inline int IND(int i,int j);
   inline int IND(int k,int l,int m,int n);
   Matrix Phi(unsigned moduliflag=0,YDeriv dy=Y0,TDeriv dt=T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
