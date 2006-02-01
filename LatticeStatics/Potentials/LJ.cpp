#include "LJ.h"

LJ::LJ(double Eps0,double Eps1,double Sigma0,double Sigma1):
   Eps0_(Eps0),Eps1_(Eps1),Sigma0_(Sigma0),Sigma1_(Sigma1)
{
}

double LJ::Eps(double NTemp,TDeriv dt)
{
   double retval;
   switch (dt)
   {
      case T0:
	 retval = Eps0_ + Eps1_*(NTemp - 1.0);
	 break;
      case DT:
	 retval = Eps1_;
	 break;
      case D2T:
	 retval = 0.0;
	 break;
      default:
	 cerr << "Error in LJ::A" << endl;
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
	 cerr << "Error in LJ::Beta" << endl;
	 exit(-1);
   }

   return retval;
}

double LJ::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double S=Sigma(NTemp),
      E=Eps(NTemp),
      r = sqrt(r2);

   double val=0;

   switch (dy)
   {
      case Y0:
	 switch (dt)
	 {
	    case T0:
	       val = 4.0*E*(pow(S/r,12.0) - pow(S/r,6.0));
	       break;
	    case DT:
	       val = 4.0*Eps(NTemp,DT)*(pow(S/r,12.0) - pow(S/r,6.0))
		  + (24.0*E*Sigma(NTemp,DT)/r)*(2.0*pow(S/r,11.0) - pow(S/r,5.0));
	       break;
	    case D2T:
	       val = 4.0*Eps(NTemp,D2T)*(pow(S/r,12.0) - pow(S/r,6.0))
		  + (24.0/r)*(2.0*Eps(NTemp,DT)*Sigma(NTemp,DT) + E*Sigma(NTemp,D2T))
		  *(2.0*pow(S/r,11.0) - pow(S/r,5.0))
		  + (24*E*pow(Sigma(NTemp,DT),2.0)/r2)*(22.0*pow(S/r,10) - 5.0*pow(S/r,4));
	       break;
	    default:
	       cerr << "Error in LJ::PairPotential" << endl;
	       exit(-1);}
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -(12.0*E*S/pow(r2,1.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0));
	       break;
	    case DT:
	       val = -(12.0/pow(r2,1.5))*(Eps(NTemp,DT)*S + E*Sigma(NTemp,DT))
		  *(2.0*pow(S/r,11) - pow(S/r,5))
		  - (12.0*E*S*Sigma(NTemp,DT)/pow(r2,2.0))
		  *(22.0*pow(S/r,10) - 5.0*pow(S/r,4.0));
	       break;
	    default:
	       cerr << "Error in LJ::PairPotential -- DY,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = (18.0*E*S/pow(r2,2.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
		  +(12.0*E*S*S/pow(r2,3.0))*(22.0*pow(S/r,10.0) - pow(S/r,4.0));
	       break;
	    case DT:
	       val = (18.0/pow(r2,2.5))*(Eps(NTemp,DT)*S + E*Sigma(NTemp,DT))
		  *(2.0*pow(S/r,11.0) - pow(S/r,5.0))
		  +(6.0/pow(r2,3.0))*(2.0*Eps(NTemp,DT)*S*S + 7.0*E*S*Sigma(NTemp,DT))
		  *(22.0*pow(S/r,10) - 5.0*pow(S/r,4))
		  +(240.0*E*S*S*Sigma(NTemp,DT)/pow(r2,3.5))
		  *(11.0*pow(S/r,9.0) - pow(S/r,3.0));
	       break;
	    default:
	       cerr << "Error in LJ::PairPotential -- D2Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D3Y:
	 switch (dt)
	 {
	    case T0:
	       val = (-45.0*E*S/pow(r2,3.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
		  -(45.0*E*S*S/pow(r2,4.0))*(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
		  -(120.0*E*S*S*S/pow(r2,4.5))*(11.0*pow(S/r,9.0) - pow(S/r,3.0));
	       break;
	    case DT:
	       cerr << "D4phi/Dy3DT Not Coded... " << endl;
	       exit(-1);
	       break;
	    default:
	       cerr << "Error in LJ::PairPotential -- D3Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D4Y:
	 switch (dt)
	 {
	    case T0:
	       val = (315.0*E*S/(2.0*pow(r2,4.5)))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
		  +(405.0*E*S*S/(2.0*pow(r2,5.0)))*(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
		  +(990.0*E*S*S*S/pow(r2,5.5))*(11.0*pow(S/r,9.0) - pow(S/r,3.0))
		  +(180.0*E*S*S*S*S/pow(r2,6.0))*(33.0*pow(S/r,8.0) - pow(S/r,2.0));
	       break;
	    case DT:
	       cerr << "D5phi/Dy4DT Not Coded... " << endl;
	       exit(-1);
	       break;
	    default:
	       cerr << "Error in LJ::PairPotential -- D4Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
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
