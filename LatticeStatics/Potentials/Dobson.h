#ifndef __Dobson
#define __Dobson

#include "PairPotentials.h"

#define D3_jFACT -8000.0
#define D3_RCUTL 1.5
#define D3_RCUTU 1.7

using namespace std;

class Dobson: public PairPotentials
{
private:
   double Eps0_,Eps1_,Sigma0_,Sigma1_,rcut_;

public:

   Dobson() {};
   Dobson(double Eps0,double Eps1,double Sigma0,double Sigma1,double rcut);
   ~Dobson() {};
   friend ostream &operator<<(ostream &out,Dobson &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "Dobson";}

   double Eps0() {return Eps0_;}
   double Eps1() {return Eps1_;}
   double Sigma0() {return Sigma0_;}
   double Sigma1() {return Sigma1_;}
   double rcut() {return rcut_;}

   void SetEps0(double Eps0) {Eps0_=Eps0;}
   void SetEps1(double Eps1) {Eps1_=Eps1;}
   void SetSigma0(double Sigma0) {Sigma0_=Sigma0;}
   void SetSigma1(double Sigma1) {Sigma1_=Sigma1;}
   void Setrcut(double rcut) {rcut_ = rcut;}
   
private:
   double Eps(double NTemp,TDeriv dt=T0);
   double Sigma(double NTemp,TDeriv dt=T0);
   double j(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   double A(double NTemp,TDeriv dt=T0);
   double B(double NTemp,TDeriv dt=T0);

};

#endif
