#ifndef RSE__RadiiMorseConstCutoff
#define RSE__RadiiMorseConstCutoff

#include "RadiiMorse.h"

using namespace std;

class RadiiMorseConstCutoff : public RadiiMorse
{
private:
   double Cutoff_;

public:
   RadiiMorseConstCutoff()
   {
   }

   RadiiMorseConstCutoff(double const& A0, double const& AT, double const& B0, double const& BT,
                    double const& Rref1, double const& Rtheta1, double const& Rref2,
                    double const& Rtheta2, double const& Cutoff);
   ~RadiiMorseConstCutoff()
   {
   }

   friend ostream& operator<<(ostream& out, RadiiMorseConstCutoff const& A);
   double PairPotential(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                        TDeriv const& dt = T0) const;
   virtual int GetNoParameters() const
   {
      return (1 + RadiiMorse::GetNoParameters());
   }

   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const
   {
      return "RadiiMorseConstCutoff";
   }

   double const& Cutoff() const
   {
      return Cutoff_;
   }

   void SetCutoff(double const& Cutoff)
   {
      Cutoff_ = Cutoff;
   }

private:
   double CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                         TDeriv const& dt = T0) const;
};

#endif
