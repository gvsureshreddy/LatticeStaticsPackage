#include "LJConstCutoff.h"
#include <cstdlib>

LJConstCutoff::LJConstCutoff(double const& Eps0, double const& Eps1, double const& Sigma0,
                             double const& Sigma1, double const& Cutoff) :
   LJ(Eps0, Eps1, Sigma0, Sigma1), Cutoff_(Cutoff)
{
}

void LJConstCutoff::SetParameters(double const* const Vals)
{
   SetCutoff(Vals[0]);
   LJ::SetParameters(&(Vals[1]));
}

double LJConstCutoff::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                     TDeriv const& dt) const
{
   double val = 0;

   switch (dy)
   {
      case Y0:
         val = -LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, Y0, dt);
         break;
      case DY:
         val = 0.0;
         break;
      case D2Y:
         val = 0.0;
         break;
      case D3Y:
      case D4Y:
         val = 0.0;
         break;
      case DYmax:
      default:
         cerr << "Error in LJConstCutoff::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJConstCutoff::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
                                    TDeriv const& dt) const
{
   if (r2 >= Cutoff_ * Cutoff_)
   {
      return 0.0;
   }
   else
   {
      return (LJ::PairPotential(NTemp, r2, dy, dt)
              + CutoffFunction(NTemp, r2, dy, dt));
   }
}

void LJConstCutoff::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream& operator<<(ostream& out, LJConstCutoff const& A)
{
   A.Print(out);
   return out;
}

