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
   virtual void Print(ostream &out) = 0;
   friend ostream &operator<<(ostream &out,PairPotentials *PP)
   {PP->Print(out); return out;}
};

#endif
