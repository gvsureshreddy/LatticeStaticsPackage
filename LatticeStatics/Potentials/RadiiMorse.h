#ifndef __RadiiMorse
#define __RadiiMorse

#include "PairPotentials.h"


class RadiiMorse: public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref_, Rtheta_, Tref_;

public:

   RadiiMorse() {};
   RadiiMorse(double A0,double B0,double Alpha,double Rref,double Rtheta,double Tref);
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
