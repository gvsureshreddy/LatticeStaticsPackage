#ifndef __NiTiPressTempPairPotLat
#define __NiTiPressTempPairPotLat

#include "UniDefTempLat.h"

#define DIM3 3

class NiTiPressTempPairPotLat : public UniDefTempLat
{
private:
   double RefLen_;
   unsigned InfluanceDist_;
   double Temp_;
   Matrix DefGrad_;
   double ShearMod_;
   double Pressure_;

   //Pair Potential data
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   enum interaction {aa,bb,ab};
   double A0_aa, B0_aa, Alpha_aa, Rref_aa, Tref_aa, Tmelt_aa;
   double A0_bb, B0_bb, Alpha_bb, Rref_bb, Tref_bb, Tmelt_bb;
   double A0_ab, B0_ab, Alpha_ab, Rref_ab, Tref_ab, Tmelt_ab;

   // Misc
   double ConvexityDX_;
   static const double Alt[DIM3][DIM3][DIM3];

public:
   // Virtual Functions required by UniDefLat
   virtual Matrix DefGrad() { return DefGrad_;}
   virtual Matrix StressDT();
   virtual void SetDefGrad(const Matrix &defgrad) { DefGrad_ = defgrad;}
   virtual double Temp() {return Temp_;}
   virtual void SetTemp(const double &temp) {Temp_ = temp;}
   double RefLen() {return RefLen_;}

   // Virtual Functions required by Lattice
   virtual double Energy();
   virtual Matrix Stress();
   virtual Matrix Stiffness();
   virtual Matrix Moduli();
   virtual int StiffnessNulity(double *Min=NULL);
   virtual void Print(ostream &out,PrintDetail flag);
   virtual void CriticalPointInfo(double Tolerance,int Width,ostream &out);
   
   // Functions provided by NiTiPressTempPairPotLat
   NiTiPressTempPairPotLat(char *datafile);
   ~NiTiPressTempPairPotLat() {}
   inline double Del(int i,int j) {return i==j;}
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,NiTiPressTempPairPotLat &A);
private:
   double PairPotential(interaction inter,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   inline double Beta(interaction inter,TDeriv dt=T0);
   inline double Rhat(interaction inter,TDeriv dt=T0);
   inline double PI(const Vector &Dx,const Vector &DX,int r,int s);
   inline double PSI(const Vector &DX,int r,int s,int t,int u);
   double pwr(const double &x,const unsigned y);
   inline int IND(int i,int j);
   inline int IND(int k,int l,int m,int n);
   Matrix Phi(unsigned moduliflag=0,YDeriv dy=Y0,TDeriv dt=T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
