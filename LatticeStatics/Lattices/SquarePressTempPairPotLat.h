#ifndef __SquarePressTempPairPotLat
#define __SquarePressTempPairPotLat

#include "UniDefTempLat.h"

class SquarePressTempPairPotLat : public UniDefTempLat
{
private:
   const static int DIM2 = 2;

   double RefLen_;
   unsigned InfluanceDist_;
   double Temp_;
   Matrix DefGrad_;
   double ShearMod_;
   double Pressure_;
   double ConvexityDX_;

   //Pair Potential data
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   double A0_, B0_, Alpha_, Rref_, Tref_, Tmelt_;

   // Misc
   static const double Alt[DIM2][DIM2];

public:
   // Virtual Functions required by UniDefLat
   virtual Matrix DefGrad() { return DefGrad_;}
   virtual Matrix StressDT() {return Phi(0,DY,DT);}
   virtual void SetDefGrad(const Matrix &defgrad) { DefGrad_ = defgrad;}
   virtual double Temp() {return Temp_;}
   virtual void SetTemp(const double &temp) {Temp_ = temp;}
   double RefLen() {return RefLen_;}

   // Virtual Functions required by Lattice
   virtual double Energy() {return Phi()[0][0];}
   virtual Matrix Stress() {return Phi(0,DY);}
   virtual Matrix Stiffness() {return Phi(0,D2Y);}
   virtual Matrix Moduli() {return Phi(1,D2Y);}
   virtual Matrix E3() {return Phi(0,D3Y);}
   virtual Matrix E4() {return Phi(0,D4Y);}
   virtual int StiffnessNulity(double *Min=NULL);
   virtual void Print(ostream &out,PrintDetail flag);
   virtual void CriticalPointInfo(const Vector &DrDt,double Tolerance,
				  char *datafile,int Width,ostream &out) {};
   
   // Functions provided by SquarePressTempPairPotLat
   SquarePressTempPairPotLat(char *datafile);
   ~SquarePressTempPairPotLat() {}
   inline double Del(int i,int j) {return i==j;}
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,SquarePressTempPairPotLat &A);
private:
   double PairPotential(double r2,YDeriv dy=Y0,TDeriv dt=T0);
   double Beta(TDeriv dt=T0);
   double Rhat(TDeriv dt=T0);
   double PI(const Vector &Dx,const Vector &DX,int r,int s);
   double PSI(const Vector &DX,int r,int s,int t,int u);
   double pwr(const double &x,const unsigned y);
   inline int IND(int i,int j);
   inline int IND(int k,int l,int m,int n);
   Matrix Phi(unsigned moduliflag=0,YDeriv dy=Y0,TDeriv dt=T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
