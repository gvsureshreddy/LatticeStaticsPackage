#ifndef __LJ
#define __LJ

#include "PairPotentials.h"

using namespace std;

class LJ: public PairPotentials
{
protected:
   double Eps0_,Eps1_,Sigma0_,Sigma1_;

public:

   LJ() {};
   LJ(double Eps0,double Eps1,double Sigma0,double Sigma1);
   ~LJ() {};
   friend ostream &operator<<(ostream &out,LJ &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "LJ";}

   double Eps0() {return Eps0_;}
   double Eps1() {return Eps1_;}
   double Sigma0() {return Sigma0_;}
   double Sigma1() {return Sigma1_;}

   void SetEps0(double Eps0) {Eps0_=Eps0;}
   void SetEps1(double Eps1) {Eps1_=Eps1;}
   void SetSigma0(double Sigma0) {Sigma0_=Sigma0;}
   void SetSigma1(double Sigma1) {Sigma1_=Sigma1;}
private:
   double Eps(double NTemp,TDeriv dt=T0);
   double Sigma(double NTemp,TDeriv dt=T0);

};

#endif
