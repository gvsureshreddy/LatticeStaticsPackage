#ifndef __TempMorse
#define __TempMorse

#include "PairPotentials.h"

using namespace std;


class TempMorse: public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref_, Tref_, Tmelt_;

public:
   TempMorse() {};
   TempMorse(double A0,double B0,double Alpha,double Rref,double Tref,double Tmelt);
   ~TempMorse() {};
   friend ostream &operator<<(ostream &out,TempMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
private:
   double Beta(double NTemp,TDeriv dt=T0);
   double Rhat(double NTemp,TDeriv dt=T0);
};

#endif
