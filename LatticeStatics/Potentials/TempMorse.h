#ifndef __TempMorse
#define __TempMorse

#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>


class TempMorse
{
private:
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   enum TDeriv {T0,DT};
   double A0_, B0_, Alpha_, Rref_, Tref_, Tmelt_;

public:
   TempMorse(double A0,double B0,double Alpha,double Rref,double Tref,double Tmelt);
   ~TempMorse() {};
   friend ostream &operator<<(ostream &out,TempMorse &A);
   double PairPotential(double Temp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
private:
   double Beta(double Temp,TDeriv dt=T0);
   double Rhat(double Temp,TDeriv dt=T0);   
};

#endif
