#ifndef __PairPotentials
#define __PairPotentials

#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>


class PairPotentials
{
public:
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y};
   enum TDeriv {T0,DT};

   ~PairPotentials() {};
   
   virtual double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0) = 0;
};

#endif
