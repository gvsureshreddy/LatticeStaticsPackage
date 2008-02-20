#ifndef __PairPotentials
#define __PairPotentials

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;


class PairPotentials
{
public:
   enum YDeriv {Y0,DY,D2Y,D3Y,D4Y,DYmax};
   enum TDeriv {T0,DT,D2T,DTmax};
   
   virtual ~PairPotentials() {};
   
   virtual double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0) = 0;
   virtual const char* Type() = 0;
   virtual void Print(ostream &out) = 0;
   friend ostream &operator<<(ostream &out,PairPotentials *PP)
   {PP->Print(out); return out;}
};

#endif
