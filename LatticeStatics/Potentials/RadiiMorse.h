#ifndef __RadiiMorse
#define __RadiiMorse

#include "PairPotentials.h"

using namespace std;

class RadiiMorse: public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref1_, Rref2_, Rtheta1_, Rtheta2_;

public:

   RadiiMorse() {};
   RadiiMorse(double A0,double B0,double Alpha,double Rref1,double Rref2,
	      double Rtheta1,double Rtheta2);
   ~RadiiMorse() {};
   friend ostream &operator<<(ostream &out,RadiiMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorse";}

   double A0() {return A0_;}
   double B0() {return B0_;}
   double Alpha() {return Alpha_;}
   double Rref1() {return Rref1_;}
   double Rref2() {return Rref2_;}
   double Rtheta1() {return Rtheta1_;}
   double Rtheta2() {return Rtheta2_;}

   double SetA0(double A0) {A0_=A0;}
   double SetB0(double B0) {B0_=B0;}
   double SetAlpha(double Alpha) {Alpha_=Alpha;}
   double SetRref1(double Rref1) {Rref1_=Rref1;}
   double SetRref2(double Rref2) {Rref2_=Rref2;}
   double SetRtheta1(double Rtheta1) {Rtheta1_=Rtheta1;}
   double SetRtheta2(double Rtheta2) {Rtheta2_=Rtheta2;}
private:
   double A(double NTemp,TDeriv dt=T0);
   double Beta(double NTemp,TDeriv dt=T0);
   double Rhat(double NTemp,TDeriv dt=T0);

};

#endif
