#include "RadiiMorseCutoff2.h"

RadiiMorseCutoff2::RadiiMorseCutoff2(double A0,double B0,double Alpha,double Rref1,double Rref2,
                                   double Rtheta1,double Rtheta2,double Cutoff):
   RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2),Cutoff_(Cutoff)
{
}

double RadiiMorseCutoff2::CutoffFunction(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double c2=Cutoff_*Cutoff_;
   double val=0;

   switch (dy)
   {
      case Y0:
         val = - (RadiiMorse::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(2.0*c2))
            *r2*r2
            + ((RadiiMorse::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/2.0)*c2 -
               RadiiMorse::PairPotential(NTemp,Cutoff_*Cutoff_,Y0,dt));
         break;
      case DY:
         val = -(RadiiMorse::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(c2))*r2;
         break;
      case D2Y:
         val = -(RadiiMorse::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(c2));
         break;
      case D3Y:
      case D4Y:
      case DYmax:
      default:
         cerr << "Error in RadiiMorseCutoff2::CutoffFunction()\n";
         exit(-1);
   }
   return val;
}

double RadiiMorseCutoff2::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   if (r2 >= Cutoff_*Cutoff_)
      return 0.0;
   else
      return (RadiiMorse::PairPotential(NTemp,r2,dy,dt)
              + CutoffFunction(NTemp,r2,dy,dt));
}

void RadiiMorseCutoff2::Print(ostream &out)
{
   int W=out.width();
   
   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; B0=" << setw(W) << B0_
       << "; Alpha=" << setw(W) << Alpha_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream &operator<<(ostream &out,RadiiMorseCutoff2 &A)
{
   A.Print(out);
   return out;
}
