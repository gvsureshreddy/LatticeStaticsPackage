#include "GVMorse.h"
#include <cstdlib>

GVMorse::GVMorse(double const& A0, double const& AT, double const& ATPow, double const& B0,
                 double const& BT, double const& BTPow, double const& Rref1, double const& Rtheta1,
                 double const& Rtheta1Pow, double const& Rref2, double const& Rtheta2,
                 double const& Rtheta2Pow) :
   RadiiMorse(A0, AT, B0, BT, Rref1, Rtheta1, Rref2, Rtheta2), ATPow_(ATPow), BTPow_(BTPow),
   Rtheta1Pow_(Rtheta1Pow), Rtheta2Pow_(Rtheta2Pow)
{
}

void GVMorse::SetParameters(double const* const Vals)
{
   SetA0(Vals[0]);
   SetAT(Vals[1]);
   SetATPow(Vals[2]);
   SetB0(Vals[3]);
   SetBT(Vals[4]);
   SetBTPow(Vals[5]);
   SetRref1(Vals[6]);
   SetRtheta1(Vals[7]);
   SetRtheta1Pow(Vals[8]);
   SetRref2(Vals[9]);
   SetRtheta2(Vals[10]);
   SetRtheta2Pow(Vals[11]);
}

double GVMorse::A(double const& NTemp, TDeriv const& dt) const
{
   double retval;
   switch (dt)
   {
      case T0:
         retval = A0_ + AT_ * (pow(NTemp, ATPow_) - 1.0);
         break;
      case DT:
         retval = AT_ * ATPow_ * pow(NTemp, ATPow_ - 1.0);
         break;
      case D2T:
         retval = AT_ * ATPow_ * (ATPow_ - 1.0) * pow(NTemp, ATPow_ - 2.0);
         break;
      default:
         cerr << "Error in GVMorse::A" << "\n";
         exit(-1);
   }

   return retval;
}

double GVMorse::Beta(double const& NTemp, TDeriv const& dt) const
{
   double retval;

   switch (dt)
   {
      case T0:
         retval = B0_ + BT_ * (pow(NTemp, BTPow_) - 1.0);
         break;
      case DT:
         retval = BT_ * BTPow_ * pow(NTemp, BTPow_ - 1.0);
         break;
      case D2T:
         retval = BT_ * BTPow_ * (BTPow_ - 1.0) * pow(NTemp, BTPow_ - 2.0);
         break;
      default:
         cerr << "Error in GVMorse::Beta" << "\n";
         exit(-1);
   }

   return retval;
}

double GVMorse::Rhat(double const& NTemp, TDeriv const& dt) const
{
   double rhat;

   switch (dt)
   {
      case T0:
         rhat = Rref1_ + Rtheta1_ * (exp(Rtheta1Pow_ * (NTemp - 1.0)) - 1.0);
         break;
      case DT:
         rhat = Rtheta1_ * exp(Rtheta1Pow_ * (NTemp - 1.0)) * Rtheta1Pow_;
         break;
      case D2T:
         rhat = Rtheta1_ * exp(Rtheta1Pow_ * (NTemp - 1.0)) * Rtheta1Pow_ * Rtheta1Pow_;
         break;
      default:
         cerr << "Error in GVMorse::Rhat" << "\n";
         exit(-1);
   }

   return rhat;
}

void GVMorse::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "A0=" << setw(W) << A0_
       << "; AT=" << setw(W) << AT_
       << "; ATPow=" << setw(W) << ATPow_
       << "; B0=" << setw(W) << B0_
       << "; BT=" << setw(W) << BT_
       << "; BTPow=" << setw(W) << BTPow_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rtheta1Pow=" << setw(W) << Rtheta1Pow_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_
       << "; Rtheta2Pow=" << setw(W) << Rtheta2Pow_;
}

ostream& operator<<(ostream& out, GVMorse const& A)
{
   A.Print(out);
   return out;
}

