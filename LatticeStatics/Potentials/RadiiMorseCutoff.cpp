#include "RadiiMorseCutoff.h"

RadiiMorseCutoff::RadiiMorseCutoff(double A0,double AT,double B0,double BT,double Rref1,
                                   double Rtheta1,double Rtheta1Pow,double Rref2,
                                   double Rtheta2,double Rtheta2Pow,double Cutoff):
   RadiiMorse(A0,AT,B0,BT,Rref1,Rtheta1,Rtheta1Pow,Rref2,Rtheta2,Rtheta2Pow),Cutoff_(Cutoff)
{
}

void RadiiMorseCutoff::SetParameters(double *Vals)
{
   SetCutoff(Vals[0]);
   RadiiMorse::SetParameters(&(Vals[1]));
}

double RadiiMorseCutoff::CutoffFunction(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double c2=Cutoff_*Cutoff_;
   double val=0;
   
   switch (dy)
   {
      case Y0:
         val = -(RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/(6*c2))*r2*r2*r2
            + ((RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/2.0)*c2
               -RadiiMorse::PairPotential(NTemp,c2,DY,dt))*r2
            + (-(RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/3.0)*c2*c2
               +RadiiMorse::PairPotential(NTemp,c2,DY,dt)*c2
               -RadiiMorse::PairPotential(NTemp,c2,Y0,dt));
         break;
      case DY:
         val = -(RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/(2*c2))*r2*r2
            + ((RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/2.0)*c2
               -RadiiMorse::PairPotential(NTemp,c2,DY,dt));
         break;
      case D2Y:
         val = -(RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/c2)*r2;
         break;
      case D3Y:
         val = -(RadiiMorse::PairPotential(NTemp,c2,D2Y,dt)/c2);
         break;
      case D4Y:
         val = 0.0;
         break;
      case DYmax:
      default:
         cerr << "Error in RadiiMorseCutoff::CutoffFunction()\n";
         exit(-1);
   }
   return val;
}

double RadiiMorseCutoff::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   if (r2 >= Cutoff_*Cutoff_)
      return 0.0;
   else
      return (RadiiMorse::PairPotential(NTemp,r2,dy,dt)
              + CutoffFunction(NTemp,r2,dy,dt));
}

void RadiiMorseCutoff::Print(ostream &out)
{
   int W=out.width();
   
   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; AT=" << setw(W) << AT_
       << "; B0=" << setw(W) << B0_
       << "; BT=" << setw(W) << BT_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rtheta1Pow=" << setw(W) << Rtheta1Pow_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_
       << "; Rtheta2Pow=" << setw(W) << Rtheta2Pow_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream &operator<<(ostream &out,RadiiMorseCutoff &A)
{
   A.Print(out);
   return out;
}
