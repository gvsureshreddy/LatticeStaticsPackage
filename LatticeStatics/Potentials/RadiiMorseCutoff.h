#ifndef __RadiiMorseCutoff
#define __RadiiMorseCutoff

#include "RadiiMorse.h"

using namespace std;

class RadiiMorseCutoff: public RadiiMorse
{
private:
   double Cutoff_;

public:

   RadiiMorseCutoff() {};
   RadiiMorseCutoff(double A0,double B0,double Alpha,double Rref1,double Rref2,
	      double Rtheta1,double Rtheta2,double Cutoff);
   ~RadiiMorseCutoff() {};
   friend ostream &operator<<(ostream &out,RadiiMorseCutoff &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorseCutoff";}

   double Cutoff() {return Cutoff_;}

   void SetCutoff(double Cutoff) {Cutoff_=Cutoff;}
private:
   double CutoffFunction(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);

};

#endif
