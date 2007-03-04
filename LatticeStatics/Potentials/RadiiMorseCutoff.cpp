#include "RadiiMorseCutoff.h"

RadiiMorseCutoff::RadiiMorseCutoff(double A0,double B0,double Alpha,double Rref1,double Rref2,
		       double Rtheta1,double Rtheta2,double Cutoff):
   RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2),Cutoff_(Cutoff)
{
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
       << "; B0=" << setw(W) << B0_
       << "; Alpha=" << setw(W) << Alpha_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream &operator<<(ostream &out,RadiiMorseCutoff &A)
{
   A.Print(out);
   return out;
}
