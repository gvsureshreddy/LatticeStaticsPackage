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
   
   double A0_,
      AT_,
      B0_,
      BT_,
      Rref1_,
      Rtheta1_,
      Rtheta1Pow_,
      Rref2_,
      Rtheta2_,
      Rtheta2Pow_;
   
public:
   
   RadiiMorse() {};
   RadiiMorse(double A0,double AT,double B0,double BT,double Rref1,double Rtheta1,
              double Rtheta1Pow,double Rref2,double Rtheta2,double Rtheta2Pow);
   ~RadiiMorse() {};
   friend ostream &operator<<(ostream &out,RadiiMorse &A);
   double PairPotential(double NTemp,double r2,YDeriv dy=Y0,TDeriv dt=T0);
   virtual int GetNoParameters() {return 10;}
   virtual void SetParameters(double *Vals);
   virtual void Print(ostream &out);
   virtual const char* Type() {return "RadiiMorse";}
   
   double A0() {return A0_;}
   double AT() {return AT_;}
   double B0() {return B0_;}
   double BT() {return BT_;}
   double Rref1() {return Rref1_;}
   double Rtheta1() {return Rtheta1_;}
   double Rtheta1Pow() {return Rtheta1Pow_;}
   double Rref2() {return Rref2_;}
   double Rtheta2() {return Rtheta2_;}
   double Rtheta2Pow() {return Rtheta2Pow_;}
   
   void SetA0(double A0) {A0_=A0;}
   void SetAT(double AT) {AT_=AT;}
   void SetB0(double B0) {B0_=B0;}
   void SetBT(double BT) {BT_=BT;}
   void SetRref1(double Rref1) {Rref1_=Rref1;}
   void SetRtheta1(double Rtheta1) {Rtheta1_=Rtheta1;}
   void SetRtheta1Pow(double Rtheta1Pow) {Rtheta1Pow_=Rtheta1Pow;}
   void SetRref2(double Rref2) {Rref2_=Rref2;}
   void SetRtheta2(double Rtheta2) {Rtheta2_=Rtheta2;}
   void SetRtheta2Pow(double Rtheta2Pow) {Rtheta2Pow_=Rtheta2Pow;}
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
