#ifndef __NiTiShuffle1TPPLat
#define __NiTiShuffle1TPPLat

#include "GenericLat.h"
#include "RadiiMorse.h"

#define DIM3 3
#define INTERNAL_ATOMS 4

class NiTiShuffle1TPPLat : public GenericLat
{
private:
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
   Matrix Phi(unsigned moduliflag=0,RadiiMorse::YDeriv dy=RadiiMorse::Y0,
	      RadiiMorse::TDeriv dt=RadiiMorse::T0);
   int FindLatticeSpacing(int iter,double dx);
   
};

#endif
