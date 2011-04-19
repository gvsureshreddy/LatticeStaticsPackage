#include "LJQuadraticCutoff.h"
#include <cstdlib>

LJQuadraticCutoff::LJQuadraticCutoff(double const& Eps0, double const& Eps1, double const& Sigma0,
                                     double const& Sigma1, double const& Cutoff) :
   LJ(Eps0, Eps1, Sigma0, Sigma1), Cutoff_(Cutoff)
{
}

void LJQuadraticCutoff::SetParameters(double const* const Vals)
{
   SetCutoff(Vals[0]);
   LJ::SetParameters(&(Vals[1]));
}

double LJQuadraticCutoff::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                         TDeriv const& dt) const
{
   double val = 0;

   double A = -(LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, D2Y, dt) * 4.0 * Cutoff_ * Cutoff_
                + LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt) * 2.0) / 2.0;
   double B = -(LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt) * 2.0 * Cutoff_) - 2 * A * Cutoff_;


   switch (dy)
   {
      case Y0:
         val = A * r2 + B* sqrt(r2) - (A * Cutoff_ * Cutoff_ + B * Cutoff_
                                       + LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, Y0, dt));
         break;
      case DY:
         val = A + B / (2.0 * sqrt(r2));
         break;
      case D2Y:
         val = -B / (4.0 * r2 * sqrt(r2));
         break;
      case D3Y:
         val = 3.0 * B / (8.0 * r2 * r2 * sqrt(r2));
         break;
      case D4Y:
         val = -15.0 * B / (16.0 * r2 * r2 * r2 * sqrt(r2));
         break;
      case DYmax:
      default:
         cerr << "Error in LJQuadraticCutoff::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJQuadraticCutoff::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
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

void LJQuadraticCutoff::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream& operator<<(ostream& out, LJQuadraticCutoff const& A)
{
   A.Print(out);
   return out;
}

