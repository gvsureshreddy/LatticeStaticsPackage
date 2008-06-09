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
   RadiiMorseCutoff(double A0,double AT,double B0,double BT,double Rref1,double Rref2,
                    double Rtheta1,double Rtheta2,double Cutoff);
   ~RadiiMorseCutoff() {};
   friend ostream &operator<<(ostream &out,RadiiMorseCutoff &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual int GetNoParameters() {return (1+RadiiMorse::GetNoParameters());}
   virtual void SetParameters(double *Vals);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorseCutoff";}
   
   double Cutoff() {return Cutoff_;}
   
   void SetCutoff(double Cutoff) {Cutoff_=Cutoff;}
private:
   double CutoffFunction(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   
};

#endif
