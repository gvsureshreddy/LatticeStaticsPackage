#ifndef __RadiiMorse
#define __RadiiMorse

#include "PairPotentials.h"

using namespace std;

class RadiiMorse: public PairPotentials
{
protected:
   int Achk_[DTmax],
      Bchk_[DTmax],
      Rchk_[DTmax],
      Gchk_[DYmax][DTmax],
      Ichk_[DTmax];
   double Aval_[DTmax],
      Bval_[DTmax],
      Rval_[DTmax],
      Gval_[DYmax][DTmax],
      Ival_[DTmax];
   
   double A0_, B0_, Alpha_, Rref1_, Rref2_, Rtheta1_, Rtheta2_;
   
public:
   
   RadiiMorse() {};
   RadiiMorse(double A0,double B0,double Alpha,double Rref1,double Rref2,
              double Rtheta1,double Rtheta2);
   ~RadiiMorse() {};
   friend ostream &operator<<(ostream &out,RadiiMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual int GetNoParameters() {return 7;}
   virtual void SetParameters(double *Vals);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorse";}
   
   double A0() {return A0_;}
   double B0() {return B0_;}
   double Alpha() {return Alpha_;}
   double Rref1() {return Rref1_;}
   double Rref2() {return Rref2_;}
   double Rtheta1() {return Rtheta1_;}
   double Rtheta2() {return Rtheta2_;}
   
   void SetA0(double A0) {A0_=A0;}
   void SetB0(double B0) {B0_=B0;}
   void SetAlpha(double Alpha) {Alpha_=Alpha;}
   void SetRref1(double Rref1) {Rref1_=Rref1;}
   void SetRref2(double Rref2) {Rref2_=Rref2;}
   void SetRtheta1(double Rtheta1) {Rtheta1_=Rtheta1;}
   void SetRtheta2(double Rtheta2) {Rtheta2_=Rtheta2;}
private:
   double A(double NTemp,TDeriv dt=T0);
   inline double a(double NTemp,TDeriv dt=T0)
   {return (Achk_[dt])? Aval_[dt] : Aval_[dt]=A(NTemp,dt);}
   
   double Beta(double NTemp,TDeriv dt=T0);
   inline double b(double NTemp,TDeriv dt=T0)
   {return (Bchk_[dt])? Bval_[dt] : Bval_[dt]=Beta(NTemp,dt);}
   
   double Rhat(double NTemp,TDeriv dt=T0);
   inline double rhat(double NTemp,TDeriv dt=T0)
   {return (Rchk_[dt])? Rval_[dt] : Rval_[dt]=Rhat(NTemp,dt);}
   
   double G(double NTemp,double r2,YDeriv dy,TDeriv dt);
   inline double g(double NTemp,double r2,YDeriv dy,TDeriv dt)
   {return (Gchk_[dy][dt])? Gval_[dy][dt] : Gval_[dy][dt]=G(NTemp,r2,dy,dt);}
   
   double I(double NTemp,TDeriv dt=T0);
   inline double i(double NTemp,TDeriv dt=T0)
   {return (Ichk_[dt])? Ival_[dt] : Ival_[dt]=I(NTemp,dt);}
};

#endif
