#include "LJLinearCutoff.h"
#include <cstdlib>

LJLinearCutoff::LJLinearCutoff(double const& Eps0, double const& Eps1, double const& Sigma0,
                               double const& Sigma1, double const& Cutoff) :
   LJ(Eps0, Eps1, Sigma0, Sigma1), Cutoff_(Cutoff)
{
}

void LJLinearCutoff::SetParameters(double const* const Vals)
{
   SetCutoff(Vals[0]);
   LJ::SetParameters(&(Vals[1]));
}

double LJLinearCutoff::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                      TDeriv const& dt) const
{
   double val = 0;

   double A = -LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt) * 2.0 * Cutoff_;

   switch (dy)
   {
      case Y0:
         val = A * sqrt(r2) - (LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, Y0, dt)
                               + A * Cutoff_);
         break;
      case DY:
         val = A / (2.0 * sqrt(r2));
         break;
      case D2Y:
         val = -A / (4.0 * r2 * sqrt(r2));
         break;
      case D3Y:
         val = 3.0 * A / (8.0 * r2 * r2 * sqrt(r2));
         break;
      case D4Y:
         val = -15.0 * A / (16.0 * r2 * r2 * r2 * sqrt(r2));
         break;
      case DYmax:
      default:
         cerr << "Error in LJLinearCutoff::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJLinearCutoff::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
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

void LJLinearCutoff::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream& operator<<(ostream& out, LJLinearCutoff const& A)
{
   A.Print(out);
   return out;
}
