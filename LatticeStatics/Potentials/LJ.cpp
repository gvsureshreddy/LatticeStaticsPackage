#include "LJ.h"

LJ::LJ(double Eps0,double Eps1,double Sigma0,double Sigma1):
   Eps0_(Eps0),Eps1_(Eps1),Sigma0_(Sigma0),Sigma1_(Sigma1)
{
}

void LJ::SetParameters(double *Vals)
{
   SetEps0(Vals[0]);
   SetEps1(Vals[1]);
   SetSigma0(Vals[2]);
   SetSigma1(Vals[3]);
}

double LJ::Eps(double NTemp,TDeriv dt)
{
   double retval;
   switch (dt)
   {
      case T0:
         retval = 4.0*(Eps0_ + Eps1_*(NTemp - 1.0));
         break;
      case DT:
         retval = 4.0*Eps1_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in LJ::Eps" << "\n";
         exit(-1);
   }
   
   return retval;
}

double LJ::Sigma(double NTemp,TDeriv dt)
{
   double retval;
   
   switch (dt)
   {
      case T0:
         retval = Sigma0_ + Sigma1_*(NTemp-1.0);
         break;
      case DT:
         retval = Sigma1_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in LJ::Sigma" << "\n";
         exit(-1);
   }
   
   return retval;
}

double LJ::G(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double val;
   
   switch (dt)
   {
      case T0:
         val = pow(s(NTemp,T0),6.0);
         break;
      case DT:
         val = 6.0*pow(s(NTemp,T0),5.0)*s(NTemp,DT);
         break;
      case D2T:
         val = 30.0*pow(s(NTemp,T0),4.0)*pow(s(NTemp,DT),2.0)
            + 6.0*pow(s(NTemp,T0),5.0)*s(NTemp,D2T);
         break;
      default:
         cerr << "Error in LJ::G dt" << "\n";
         exit(-1);
         break;
   }
   
   switch (dy)
   {
      case Y0:
         val *= pow(r2,-3.0);
         break;
      case DY:
         val *= -3.0*pow(r2,-4.0);
         break;
      case D2Y:
         val *= 12.0*pow(r2,-5.0);
         break;
      case D3Y:
         val *= -60.0*pow(r2,-6.0);
         break;
      case D4Y:
         val *= 360.0*pow(r2,-7.0);
         break;
      default:
         cerr << "Error in LJ::G dy" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}

double LJ::H(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double val;
   
   switch (dt)
   {
      case T0:
         val = pow(s(NTemp,T0),12.0);
         break;
      case DT:
         val = 12.0*pow(s(NTemp,T0),11.0)*s(NTemp,DT);
         break;
      case D2T:
         val = 132.0*pow(s(NTemp,T0),10.0)*pow(s(NTemp,DT),2.0)
            + 12.0*pow(s(NTemp,T0),11.0)*s(NTemp,D2T);
         break;
      default:
         cerr << "Error in LJ::H dt" << "\n";
         exit(-1);
         break;
   }
   
   switch (dy)
   {
      case Y0:
         val *= pow(r2,-6.0);
         break;
      case DY:
         val *= -6.0*pow(r2,-7.0);
         break;
      case D2Y:
         val *= 42.0*pow(r2,-8.0);
         break;
      case D3Y:
         val *= -336.0*pow(r2,-9.0);
         break;
      case D4Y:
         val *= 3024.0*pow(r2,-10.0);
         break;
      default:
         cerr << "Error in LJ::H dy" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}

double LJ::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   for (int i=0;i<DTmax;++i)
   {
      EpsChk_[i]=
         SigmaChk_[i]=0;
      for (int j=0;j<DYmax;++j)
      {
         Gchk_[j][i] =
            Hchk_[j][i] = 0;
      }
   }
   
   double val=0.0;
   switch (dt)
   {
      case T0:
         val = e(NTemp)*(h(NTemp,r2,dy,T0) - g(NTemp,r2,dy,T0));
         break;
      case DT:
         val = e(NTemp,DT)*(h(NTemp,r2,dy,T0) - g(NTemp,r2,dy,T0))
            + e(NTemp)*(h(NTemp,r2,dy,DT) - g(NTemp,r2,dy,DT));
         break;
      case D2T:
         val = e(NTemp,D2T)*(h(NTemp,r2,dy,T0) - g(NTemp,r2,dy,T0))
            + 2.0*e(NTemp,DT)*(h(NTemp,r2,dy,DT) - g(NTemp,r2,dy,DT))
            + e(NTemp)*(h(NTemp,r2,dy,D2T) - g(NTemp,r2,dy,D2T));
         break;
      default:
         cerr << "Error in LJ::PairPotential" << "\n";
         exit(-1);
   }
   
   return val;
}

void LJ::Print(ostream &out)
{
   int W=out.width();
   
   out.width(0);
   
   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_;
}

ostream &operator<<(ostream &out,LJ &A)
{
   A.Print(out);
   return out;
}
