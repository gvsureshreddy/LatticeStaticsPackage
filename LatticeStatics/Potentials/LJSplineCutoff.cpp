#include "LJSplineCutoff.h"
#include <cstdlib>

LJSplineCutoff::LJSplineCutoff(double const& Eps0, double const& Eps1, double const& Sigma0,
                               double const& Sigma1, double const& CutoffStart,
                               double const& CutoffEnd) :
   LJ(Eps0, Eps1, Sigma0, Sigma1), CutoffStart_(CutoffStart), CutoffEnd_(CutoffEnd),
   CoeffsMat_(3, 3), CoeffsVec_(3), Coeffs_(3)
{
   double diff = CutoffStart_ - CutoffEnd_;

   CoeffsMat_[0][0] = diff * diff * diff * diff * diff;
   CoeffsMat_[0][1] = diff * diff * diff * diff;
   CoeffsMat_[0][2] = diff * diff * diff;
   CoeffsMat_[1][0] = 5.0 * CoeffsMat_[0][1];
   CoeffsMat_[1][1] = 4.0 * CoeffsMat_[0][2];
   CoeffsMat_[1][2] = 3.0 * diff * diff;
   CoeffsMat_[2][0] = 20.0 * CoeffsMat_[0][2];
   CoeffsMat_[2][1] = 12.0 * diff * diff;
   CoeffsMat_[2][2] = 6.0 * diff;

   CoeffsMat_ = CoeffsMat_.Inverse();
}

void LJSplineCutoff::SetParameters(double const* const Vals)
{
   SetCutoffStart(Vals[0]);
   SetCutoffEnd(Vals[1]);
   LJ::SetParameters(&(Vals[2]));
}

double LJSplineCutoff::CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy,
                                      TDeriv const& dt) const
{
   double val = 0;

   CoeffsVec_[0] = LJ::PairPotential(NTemp, CutoffStart_ * CutoffStart_, Y0, dt);
   CoeffsVec_[1] = LJ::PairPotential(NTemp, CutoffStart_ * CutoffStart_, DY, dt) * 2.0 * CutoffStart_;
   CoeffsVec_[2] =
      LJ::PairPotential(NTemp, CutoffStart_ * CutoffStart_, D2Y, dt) * 4.0 * CutoffStart_ * CutoffStart_
      + LJ::PairPotential(NTemp, CutoffStart_ * CutoffStart_, DY, dt) * 2.0;

   for (int i = 0; i < 3; ++i)
   {
      Coeffs_[i] = 0.0;
      for (int j = 0; j < 3; ++j)
      {
         Coeffs_[i] += CoeffsMat_[i][j] * CoeffsVec_[j];
      }
   }

   double x = sqrt(r2) - CutoffEnd_;

   switch (dy)
   {
      case Y0:
         val = Coeffs_[0] * x * x * x * x * x + Coeffs_[1] * x * x * x * x + Coeffs_[2] * x * x * x;
         break;
      case DY:
         val = (5.0 * Coeffs_[0] * x * x * x * x + 4.0 * Coeffs_[1] * x * x * x + 3.0 * Coeffs_[2] * x * x)
               / (2.0 * sqrt(r2));
         break;
      case D2Y:
         val = (20.0 * Coeffs_[0] * x * x * x + 12.0 * Coeffs_[1] * x * x + 6.0 * Coeffs_[2] * x)
               / (4.0 * r2)
               - (5.0 * Coeffs_[0] * x * x * x * x + 4.0 * Coeffs_[1] * x * x * x + 3.0 * Coeffs_[2] * x * x)
               / (4.0 * r2 * sqrt(r2));
         break;
      case D3Y:
         val = (60.0 * Coeffs_[0] * x * x + 24.0 * Coeffs_[1] * x + 6.0 * Coeffs_[2])
               / (8.0 * r2 * sqrt(r2))
               - (20.0 * Coeffs_[0] * x * x * x + 12.0 * Coeffs_[1] * x * x + 6.0 * Coeffs_[2] * x)
               * (1.0 / (4.0 * r2 * r2) + 1.0 / (4.0 * r2 * sqrt(r2)))
               + 3.0 * (5.0 * Coeffs_[0] * x * x * x * x + 4.0 * Coeffs_[1] * x * x * x + 3.0 * Coeffs_[2] * x * x)
               / (8.0 * r2 * r2 * sqrt(r2));
         break;
      case D4Y:
         val = (120.0 * Coeffs_[0] * x + 24.0 * Coeffs_[1])
               / (16.0 * r2 * r2)
               - (60.0 * Coeffs_[0] * x * x + 24.0 * Coeffs_[1] * x + 6.0 * Coeffs_[2])
               * (3.0 / (16.0 * r2 * r2 * sqrt(r2)) + 1.0 / (8.0 * r2 * sqrt(r2)) + 1.0 / (8.0 * r2 * r2))
               + (20.0 * Coeffs_[0] * x * x * x + 12.0 * Coeffs_[1] * x * x + 6.0 * Coeffs_[2] * x)
               * (11.0 / (16.0 * r2 * r2 * r2) + 3.0 / (8.0 * r2 * r2 * sqrt(r2)))
               - 15.0 * (5.0 * Coeffs_[0] * x * x * x * x + 4.0 * Coeffs_[1] * x * x * x + 3.0 * Coeffs_[2] * x * x)
               / (16.0 * r2 * r2 * r2 * sqrt(r2));
         break;
      case DYmax:
      default:
         cerr << "Error in LJSplineCutoff::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJSplineCutoff::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
                                     TDeriv const& dt) const
{
   if (r2 < CutoffStart_ * CutoffStart_)
   {
      return LJ::PairPotential(NTemp, r2, dy, dt);
   }
   else if (r2 < CutoffEnd_ * CutoffEnd_)
   {
      return CutoffFunction(NTemp, r2, dy, dt);
   }
   else
   {
      return 0.0;
   }
}

void LJSplineCutoff::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; CutoffStart=" << setw(W) << CutoffStart_
       << "; CutoffEnd=" << setw(W) << CutoffEnd_;
}

ostream& operator<<(ostream& out, LJSplineCutoff const& A)
{
   A.Print(out);
   return out;
}
