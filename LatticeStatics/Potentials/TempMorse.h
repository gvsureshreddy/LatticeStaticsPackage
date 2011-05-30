#ifndef RSE__TempMorse
#define RSE__TempMorse

#include "PairPotentials.h"

using namespace std;


class TempMorse : public PairPotentials
{
private:
   double A0_, B0_, Alpha_, Rref_, Tref_, Tmelt_;

public:
   TempMorse()
   {
   }

   TempMorse(double const& A0, double const& B0, double const& Alpha, double const& Rref,
             double const& Tref, double const& Tmelt);
   ~TempMorse()
   {
   }

   friend ostream& operator<<(ostream& out, TempMorse const& A);
   double PairPotential(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                        TDeriv const& dt = T0) const;
   virtual int GetNoParameters() const
   {
      return 6;
   }

   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const
   {
      return "TempMorse";
   }

   double const& A0() const
   {
      return A0_;
   }

   double const& B0() const
   {
      return B0_;
   }

   double const& Alpha() const
   {
      return Alpha_;
   }

   double const& Rref() const
   {
      return Rref_;
   }

   double const& Tref() const
   {
      return Tref_;
   }

   double const& Tmelt() const
   {
      return Tmelt_;
   }

   void SetA0(double const& A0)
   {
      A0_ = A0;
   }

   void SetB0(double const& B0)
   {
      B0_ = B0;
   }

   void SetAlpha(double const& Alpha)
   {
      Alpha_ = Alpha;
   }

   void SetRref(double const& Rref)
   {
      Rref_ = Rref;
   }

   void SetTref(double const& Tref)
   {
      Tref_ = Tref;
   }

   void SetTmelt(double const& Tmelt)
   {
      Tmelt_ = Tmelt;
   }

private:
   double Beta(double const& NTemp, TDeriv const& dt = T0) const;
   double Rhat(double const& NTemp, TDeriv const& dt = T0) const;
};

#endif
