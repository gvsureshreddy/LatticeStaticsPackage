#include "LJCutoff.h"

LJCutoff::LJCutoff(double Eps0,double Eps1,double Sigma0,double Sigma1,double Cutoff):
   LJ(Eps0,Eps1,Sigma0,Sigma1),Cutoff_(Cutoff)
{
}

double LJCutoff::CutoffFunction(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double val=0;
   
   switch (dy)
   {
      case Y0:
         val = - (LJ::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(2.0*Cutoff_*Cutoff_))*r2*r2
            + ((LJ::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/2.0)*Cutoff_*Cutoff_ -
               LJ::PairPotential(NTemp,Cutoff_*Cutoff_,Y0,dt));
         break;
      case DY:
         val = -(LJ::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(Cutoff_*Cutoff_))*r2;
         break;
      case D2Y:
         val = -(LJ::PairPotential(NTemp,Cutoff_*Cutoff_,DY,dt)/(Cutoff_*Cutoff_));
         break;
      case D3Y:
      case D4Y:
         val = 0.0;
         break;
      case DYmax:
      default:
         cerr << "Error in LJCutoff::CutoffFunction()\n";
         exit(-1);
         break;
   }
   return val;
}

double LJCutoff::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   if (r2 >= Cutoff_*Cutoff_)
      return 0.0;
   else
      return (LJ::PairPotential(NTemp,r2,dy,dt)
              + CutoffFunction(NTemp,r2,dy,dt));
}

void LJCutoff::Print(ostream &out)
{
   int W=out.width();
   
   out.width(0);
   
   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; Cutoff=" << setw(W) << Cutoff_;
}

ostream &operator<<(ostream &out,LJCutoff &A)
{
   A.Print(out);
   return out;
}
