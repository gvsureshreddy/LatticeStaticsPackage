#ifndef __RadiiMorseCutoff2
#define __RadiiMorseCutoff2

#include "RadiiMorse.h"

using namespace std;

class RadiiMorseCutoff2: public RadiiMorse
{
private:
   double Cutoff_;
   
public:
   
   RadiiMorseCutoff2() {};
   RadiiMorseCutoff2(double A0,double AT,double B0,double BT,double Rref1,double Rtheta1,
                     double Rtheta1Pow,double Rref2,double Rtheta2,double Rtheta2Pow,
                     double Cutoff);
   ~RadiiMorseCutoff2() {};
   friend ostream &operator<<(ostream &out,RadiiMorseCutoff2 &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual int GetNoParameters() {return (1+RadiiMorse::GetNoParameters());}
   virtual void SetParameters(double *Vals);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorseCutoff2";}
   
   double Cutoff() {return Cutoff_;}
   
   void SetCutoff(double Cutoff) {Cutoff_=Cutoff;}
private:
   double CutoffFunction(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   
};

#endif
