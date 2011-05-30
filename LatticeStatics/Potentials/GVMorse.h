#ifndef __GVMorse
#define __GVMorse

#include "RadiiMorse.h"

using namespace std;

class GVMorse : public RadiiMorse
{
private:
   double
      ATPow_,
      BTPow_,
      Rtheta1Pow_,
      Rtheta2Pow_;

public:
   GVMorse()
   {
   }

   GVMorse(double const& A0, double const& AT, double const& ATPow, double const& B0,
           double const& BT, double const& BTPow, double const& Rref1, double const& Rtheta1,
           double const& Rtheta1Pow, double const& Rref2, double const& Rtheta2,
           double const& Rtheta2Pow);
   ~GVMorse()
   {
   }

   friend ostream& operator<<(ostream& out, GVMorse const& A);
   virtual int GetNoParameters() const
   {
      return 12;
   }

   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const
   {
      return "GVMorse";
   }

   double const& ATPow() const
   {
      return ATPow_;
   }

   double const& BTPow() const
   {
      return BTPow_;
   }

   double const& Rtheta1Pow() const
   {
      return Rtheta1Pow_;
   }

   double const& Rtheta2Pow() const
   {
      return Rtheta2Pow_;
   }

   void SetATPow(double const& ATPow)
   {
      ATPow_ = ATPow;
   }

   void SetBTPow(double const& BTPow)
   {
      BTPow_ = BTPow;
   }

   void SetRtheta1Pow(double const& Rtheta1Pow)
   {
      Rtheta1Pow_ = Rtheta1Pow;
   }

   void SetRtheta2Pow(double const& Rtheta2Pow)
   {
      Rtheta2Pow_ = Rtheta2Pow;
   }

private:
   virtual double A(double const& NTemp, TDeriv const& dt = T0) const;
   virtual double Beta(double const& NTemp, TDeriv const& dt = T0) const;
   virtual double Rhat(double const& NTemp, TDeriv const& dt = T0) const;
};

#endif
