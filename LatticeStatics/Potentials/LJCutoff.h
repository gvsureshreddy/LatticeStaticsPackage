#ifndef __LJCutoff
#define __LJCutoff

#include "LJ.h"

using namespace std;

class LJCutoff: public LJ
{
private:
   double Cutoff_;

public:

   LJCutoff() {};
   LJCutoff(double Eps0,double Eps1,double Sigma0,double Sigma1,double Cutoff);
   ~LJCutoff() {};
   friend ostream &operator<<(ostream &out,LJCutoff &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "LJCutoff";}

   double Cutoff() {return Cutoff_;}

   void SetCutoff(double Cutoff) {Cutoff_=Cutoff;}
private:

   double CutoffFunction(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);

};

#endif
