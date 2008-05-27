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
   virtual int GetNoParameters() {return 6;}
   virtual void SetParameters(double *Vals);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "TempMorse";}
   
   double A0() {return A0_;}
   double B0() {return B0_;}
   double Alpha() {return Alpha_;}
   double Rref() {return Rref_;}
   double Tref() {return Tref_;}
   double Tmelt() {return Tmelt_;}
   
   void SetA0(double A0) {A0_=A0;}
   void SetB0(double B0) {B0_=B0;}
   void SetAlpha(double Alpha) {Alpha_=Alpha;}
   void SetRref(double Rref) {Rref_=Rref;}
   void SetTref(double Tref) {Tref_=Tref;}
   void SetTmelt(double Tmelt) {Tmelt_=Tmelt;}
   
private:
   double Beta(double NTemp,TDeriv dt=T0);
   double Rhat(double NTemp,TDeriv dt=T0);
};

#endif
