#ifndef __NiTiShuffle1TPPLat
#define __NiTiShuffle1TPPLat

#include "GenericLat.h"
#include "RadiiMorse.h"

class NiTiShuffle1TPPLat : public GenericLat
{
private:
   const static int DIM3 = 3;
   const static int INTERNAL_ATOMS = 4;

   double RefLen_;
   unsigned InfluanceDist_;
   double NTemp_;
   // DOF[i] = [U11 U22 U33 U12 U13 U12 S]
   Vector DOF_;
   Vector LatticeVec_[DIM3];
   double ShearMod_;
   double Pressure_;
   Vector BodyForce_[INTERNAL_ATOMS];

   //Pair Potential data
   enum interaction {aa,bb,ab,NOINTERACTIONS};
   RadiiMorse Potential_[NOINTERACTIONS];

   // Misc
   double ConvexityDX_;
   static const double Alt[DIM3][DIM3][DIM3];

public:
   double RefLen() {return RefLen_;}

   // Virtual Functions required by GenericLat
   Vector DOF() {return DOF_;}
   Matrix StressDT() {return Phi(0,PairPotentials::DY,PairPotentials::DT);}
   Matrix StiffnessDT() {return Phi(0,PairPotentials::D2Y,PairPotentials::DT);}
   void SetDOF(const Vector &dof) { DOF_ = dof;}
   double Temp() {return NTemp_;}
   void SetTemp(const double &Ntemp) {NTemp_ = Ntemp;}
   //
   virtual void CriticalPointInfo(const Vector &DrDt,double Tolerance,
				  char *datafile,int Width,ostream &out) {};

   // Virtual Functions required by Lattice
   virtual double Energy() {return Phi()[0][0];}
   virtual Matrix Stress() {return Phi(0,PairPotentials::DY);}
   virtual Matrix Stiffness() {return Phi(0,PairPotentials::D2Y);}
   virtual Matrix Moduli() {return Phi(1,PairPotentials::D2Y);}
   virtual Matrix E3() {return Phi(0,PairPotentials::D3Y);}
   virtual Matrix E4() {return Phi(0,PairPotentials::D4Y);}
   virtual void CriticalPointInfo(const Matrix &DrDu,
				  double Tolerance,int Width,ostream &out) {}
   virtual void Print(ostream &out,PrintDetail flag);
   
   // Functions provided by NiTiShuffle1TPPLat
   NiTiShuffle1TPPLat(char *datafile);
   ~NiTiShuffle1TPPLat() {}
   inline double Del(int i,int j) {return i==j;}
   Vector BodyForce(int i) { return BodyForce_[i]; }
   double Pressure() const { return Pressure_;}
   double SetPressure(double &p) { Pressure_ = p;}
   double ShearMod() const { return ShearMod_;}
   friend ostream &operator<<(ostream &out,NiTiShuffle1TPPLat &A);
private:
   void GetLatticeVectorInfo(double *SX,double *DXPrime,interaction &Inter,
			     int p,int q);
   inline double PI(const Vector &Dx,const Vector &DX,int r,int s);
   inline double PSI(const Vector &DX,int r,int s,int t,int u);
   double pwr(const double &x,const unsigned y);
   inline int IND(int i,int j);
   inline int IND(int k,int l,int m,int n);
   Matrix Phi(unsigned moduliflag=0,PairPotentials::YDeriv dy=PairPotentials::Y0,
	      PairPotentials::TDeriv dt=PairPotentials::T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
