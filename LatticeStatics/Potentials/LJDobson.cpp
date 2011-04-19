#include "LJDobson.h"
#include <cstdlib>

LJDobson::LJDobson(double const& Eps0, double const& Eps1, double const& Sigma0,
                   double const& Sigma1, double const& Cutoff) :
   LJ(Eps0, Eps1, Sigma0, Sigma1), Cutoff_(Cutoff)
{
}

void LJDobson::SetParameters(double const* const Vals)
{
   SetCutoff(Vals[0]);
   LJ::SetParameters(&(Vals[1]));
}

double LJDobson::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                TDeriv const& dt) const
{
   double val = 0;

   double A, B;

   switch (dy)
   {
      case Y0:
         A = -(1.0 / (2.0 * Cutoff_)) * (LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt) * 2.0 * Cutoff_);
         B = (Cutoff_ / 2.0) * (LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt)) -
             LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, Y0, dt);
         val = A * r2 + B;
         break;
      case DY:
         A = -(1.0 / (2.0 * Cutoff_)) * (LJ::PairPotential(NTemp, Cutoff_ * Cutoff_, DY, dt) * 2.0 * Cutoff_);
         val = A;
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
         cerr << "Error in LJDobson::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJDobson::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
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

void LJDobson::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream& operator<<(ostream& out, LJDobson const& A)
{
   A.Print(out);
   return out;
}

