#ifndef __RadiiMorse
#define __RadiiMorse

#include "PairPotentials.h"


class RadiiMorse: public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref_, Rtheta_, Tref_, Tmelt_;

public:

   RadiiMorse() {};
   RadiiMorse(double A0,double B0,double Alpha,double Rref,double Rtheta,
	      double Tref,double Tmelt);
   ~RadiiMorse() {};
   friend ostream &operator<<(ostream &out,RadiiMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
private:
   double Beta(double NTemp,TDeriv dt=T0);
   double Rhat(double NTemp,TDeriv dt=T0);   
};

#endif