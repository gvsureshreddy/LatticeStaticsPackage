#ifndef RSE__RadiiMorse
#define RSE__RadiiMorse

#include "PairPotentials.h"

using namespace std;

class RadiiMorse: public PairPotentials
{
protected:
   mutable int Achk_[DTmax];
   mutable int Bchk_[DTmax];
   mutable int Rchk_[DTmax];
   mutable int Gchk_[DYmax][DTmax];
   mutable int Ichk_[DTmax];
   mutable double Aval_[DTmax];
   mutable double Bval_[DTmax];
   mutable double Rval_[DTmax];
   mutable double Gval_[DYmax][DTmax];
   mutable double Ival_[DTmax];
   
   double A0_,
      AT_,
      B0_,
      BT_,
      Rref1_,
      Rtheta1_,
      Rref2_,
      Rtheta2_;
   
public:
   
   RadiiMorse() {};
   RadiiMorse(double const& A0,double const& AT,double const& B0,double const& BT,
              double const& Rref1,double const& Rtheta1,double const& Rref2,
              double const& Rtheta2);
   ~RadiiMorse() {};
   friend ostream& operator<<(ostream& out,RadiiMorse const& A);
   double PairPotential(double const& NTemp,double const& r2,YDeriv const& dy=Y0,
                        TDeriv const& dt=T0) const;
   virtual int GetNoParameters() const {return 8;}
   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const {return "RadiiMorse";}
   
   double const& A0() const {return A0_;}
   double const& AT() const {return AT_;}
   double const& B0() const {return B0_;}
   double const& BT() const {return BT_;}
   double const& Rref1() const {return Rref1_;}
   double const& Rtheta1() const {return Rtheta1_;}
   double const& Rref2() const {return Rref2_;}
   double const& Rtheta2() const {return Rtheta2_;}
   
   void SetA0(double const& A0) {A0_=A0;}
   void SetAT(double const& AT) {AT_=AT;}
   void SetB0(double const& B0) {B0_=B0;}
   void SetBT(double const& BT) {BT_=BT;}
   void SetRref1(double const& Rref1) {Rref1_=Rref1;}
   void SetRtheta1(double const& Rtheta1) {Rtheta1_=Rtheta1;}
   void SetRref2(double const& Rref2) {Rref2_=Rref2;}
   void SetRtheta2(double const& Rtheta2) {Rtheta2_=Rtheta2;}
private:
   virtual double A(double const& NTemp,TDeriv const& dt=T0) const;
   inline double a(double const& NTemp,TDeriv const& dt=T0) const
   {return (Achk_[dt])? Aval_[dt] : Aval_[dt]=A(NTemp,dt);}
   
   virtual double Beta(double const& NTemp,TDeriv const& dt=T0) const;
   inline double b(double const& NTemp,TDeriv const& dt=T0) const
   {return (Bchk_[dt])? Bval_[dt] : Bval_[dt]=Beta(NTemp,dt);}
   
   virtual double Rhat(double const& NTemp,TDeriv const& dt=T0) const;
   inline double rhat(double const& NTemp,TDeriv const& dt=T0) const
   {return (Rchk_[dt])? Rval_[dt] : Rval_[dt]=Rhat(NTemp,dt);}
   
   double G(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const;
   inline double g(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const
   {return (Gchk_[dy][dt])? Gval_[dy][dt] : Gval_[dy][dt]=G(NTemp,r2,dy,dt);}
   
   double I(double const& NTemp,TDeriv const& dt=T0) const;
   inline double i(double const& NTemp,TDeriv const& dt=T0) const
   {return (Ichk_[dt])? Ival_[dt] : Ival_[dt]=I(NTemp,dt);}
};

#endif
