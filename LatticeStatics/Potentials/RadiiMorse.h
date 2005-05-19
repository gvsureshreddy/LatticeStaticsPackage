#ifndef __RadiiMorse
#define __RadiiMorse

#include "PairPotentials.h"

using namespace std;

class RadiiMorse: public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref1_, Rref2_, Rtheta1_, Rtheta2_, Tref_;

public:

   RadiiMorse() {};
   RadiiMorse(double A0,double B0,double Alpha,double Rref1,double Rref2,
	      double Rtheta1,double Rtheta2,double Tref);
   ~RadiiMorse() {};
   friend ostream &operator<<(ostream &out,RadiiMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
private:
   double A(double NTemp,TDeriv dt=T0);
   double Beta(double NTemp,TDeriv dt=T0);
   double Rhat(double NTemp,TDeriv dt=T0);

};

#endif
