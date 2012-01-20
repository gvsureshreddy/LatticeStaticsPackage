#include "RadiiMorseConstCutoff.h"
#include <cstdlib>

RadiiMorseConstCutoff::RadiiMorseConstCutoff(double const& A0, double const& AT, double const& B0,
                                   double const& BT, double const& Rref1, double const& Rtheta1,
                                   double const& Rref2, double const& Rtheta2,
                                   double const& Cutoff) :
   RadiiMorse(A0, AT, B0, BT, Rref1, Rtheta1, Rref2, Rtheta2), Cutoff_(Cutoff)
{
}

void RadiiMorseConstCutoff::SetParameters(double const* const Vals)
{
   SetCutoff(Vals[0]);
   RadiiMorse::SetParameters(&(Vals[1]));
}

double RadiiMorseConstCutoff::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                             TDeriv const& dt) const
{
   double c2 = Cutoff_ * Cutoff_;
   double val = 0;

   switch (dy)
   {
      case Y0:
         val = -(RadiiMorse::PairPotential(NTemp, c2, dy, dt));
         break;
      case DY:
         val = 0.0;
         break;
      case D2Y:
         val = 0.0;
         break;
      case D3Y:
         val = 0.0;
         break;
      case D4Y:
         val = 0.0;
         break;
      case D5Y:
         val = 0.0;
         break;
      case D6Y:
         val = 0.0;
         break;
      case DYmax:
      default:
         cerr << "Error in RadiiMorseConstCutoff::CutoffFunction()\n";
         exit(-1);
   }
   return val;
}

double RadiiMorseConstCutoff::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
                                       TDeriv const& dt) const
{
   if (r2 >= Cutoff_ * Cutoff_)
   {
      return 0.0;
   }
   else
   {
      return (RadiiMorse::PairPotential(NTemp, r2, dy, dt)
              + CutoffFunction(NTemp, r2, dy, dt));
   }
}

void RadiiMorseConstCutoff::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "A0=" << setw(W) << A0_
       << "; AT=" << setw(W) << AT_
       << "; B0=" << setw(W) << B0_
       << "; BT=" << setw(W) << BT_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream& operator<<(ostream& out, RadiiMorseConstCutoff const& A)
{
   A.Print(out);
   return out;
}
