#include "RadiiMorse.h"

RadiiMorse::RadiiMorse(double A0,double B0,double Alpha,double Rref1,double Rref2,
		       double Rtheta1,double Rtheta2):
   A0_(A0),B0_(B0),Alpha_(Alpha),Rref1_(Rref1),Rref2_(Rref2),
   Rtheta1_(Rtheta1),Rtheta2_(Rtheta2)
{
}

double RadiiMorse::A(double NTemp,TDeriv dt)
{
   double retval;
   switch (dt)
   {
      case T0:
	 retval = A0_;
	 break;
      case DT:
	 retval = 0.0;
	 break;
      case D2T:
	 retval = 0.0;
	 break;
      default:
	 cerr << "Error in RadiiMorse::A" << endl;
	 exit(-1);
   }

   return retval;
}

double RadiiMorse::Beta(double NTemp,TDeriv dt)
{
   double retval;
   
   switch (dt)
   {
      case T0:
	 retval = B0_ + Alpha_*(NTemp-1.0);
	 break;
      case DT:
	 retval = Alpha_;
	 break;
      case D2T:
	 retval = 0.0;
	 break;
      default:
	 cerr << "Error in RadiiMorse::Beta" << endl;
	 exit(-1);
   }

   return retval;
}

double RadiiMorse::Rhat(double NTemp,TDeriv dt)
{
   double rhat;

   switch (dt)
   {
      case T0:
	 rhat = (Rref1_ + Rtheta1_*(NTemp - 1.0))*(Rref2_ + Rtheta2_*(NTemp - 1.0));
	 break;
      case DT:
	 rhat = Rtheta1_*(Rref2_ + Rtheta2_*(NTemp - 1.0))
	    + Rtheta2_*(Rref1_ + Rtheta1_*(NTemp - 1.0));
	 break;
      case D2T:
	 rhat = Rtheta1_*Rtheta2_ + Rtheta2_*Rtheta1_;
	 break;
      default:
	 cerr << "Error in RadiiMorse::Rhat" << endl;
	 exit(-1);
   }

   return rhat;
}

double RadiiMorse::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double At=A(NTemp),
      //double At(NTemp),
      beta=Beta(NTemp),
      rhat=Rhat(NTemp),
      r = sqrt(r2),
      Exp_temp=exp(-beta*(r/rhat - 1.0));

   double val=0;

   switch (dy)
   {
      case Y0:
	 switch (dt)
	 {
	    case T0:
	       val = At*Exp_temp*(Exp_temp - 2.0);
	       break;
	    case DT:
	       val = A(NTemp,DT)*Exp_temp*(Exp_temp - 2.0)
		  - 2.0*At*(Beta(NTemp,DT)*(r/rhat - 1.0)-beta*Rhat(NTemp,DT)*(r/(rhat*rhat)))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	    case D2T:
	       val = A(NTemp,D2T)*Exp_temp*(Exp_temp-2.0)
		  - At*(Beta(NTemp,DT)*(r/rhat - 1.0) - (beta*Rhat(NTemp,DT)*r)/(rhat*rhat))
		  *(Beta(NTemp,DT)*(r/rhat - 1.0) - (beta*Rhat(NTemp,DT)*r)/(rhat*rhat))
		  *Exp_temp*(4.0*Exp_temp - 2.0)
		  - At*(Beta(NTemp,D2T)*(r/rhat - 1.0)
			- Beta(NTemp,DT)*(2.0*r*Rhat(NTemp,DT)/(rhat*rhat))
			- beta*(r*Rhat(NTemp,D2T)/(rhat*rhat)
				- 2.0*r*Rhat(NTemp,DT)*Rhat(NTemp,DT)/(rhat*rhat*rhat)))
		  *Exp_temp*(2.0*Exp_temp - 2.0);
	       break;
	    default:
	       cerr << "Error in RadiiMorese::PairPotential" << endl;
	       exit(-1);}
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -At*(beta/(r*rhat))*Exp_temp*(Exp_temp - 1.0);
	       break;
	    case DT:
	       val = -A(NTemp,DT)*
		  (beta/(rhat*r))*Exp_temp*(Exp_temp - 1.0)
		  + At*Exp_temp*
		  (((beta*Rhat(NTemp,DT) - Beta(NTemp,DT)*rhat)/(rhat*rhat*r))*
		   (Exp_temp - 1.0)
		   +(beta/(rhat*r))*(Beta(NTemp,DT)*(r/rhat - 1.0)
				     - beta*r*Rhat(NTemp,DT)/(rhat*rhat))
		   *(2.0*Exp_temp - 1.0));
	       break;
	    default:
	       cerr << "Error in RadiiMorse::PairPotential -- DY,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = At*(beta/(2.0*rhat*r2))
		  *((beta/rhat + 1.0/r)*Exp_temp*(Exp_temp - 1.0)
		    + (beta/rhat)*Exp_temp*Exp_temp);
	       break;
	    case DT:
	       val = (A(NTemp,DT))*(beta/(2*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +At*((Beta(NTemp,DT)*rhat
			- beta*Rhat(NTemp,DT))/
		       (2.0*rhat*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +At*(beta/(2.0*rhat*r2))
		  *(((Beta(NTemp,DT)*rhat - beta*Rhat(NTemp,DT))/(rhat*rhat))*
		    (2.0*Exp_temp*Exp_temp - Exp_temp)
		    +(beta*Rhat(NTemp,DT)*r/(rhat*rhat)
		      - Beta(NTemp,DT)*(r/rhat - 1.0))
		    *(2.0*beta/rhat*Exp_temp*Exp_temp +
		      (beta/rhat + 1.0/r)
		      *(2.0*Exp_temp*Exp_temp - Exp_temp)));
	       break;
	    default:
	       cerr << "Error in RadiiMorse::PairPotential -- D2Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D3Y:
	 switch (dt)
	 {
	    case T0:
	       val = -At*(beta/(2.0*rhat*r2))
		  *((3.0*beta/(2.0*rhat*r2) + 3.0/(2.0*r2*r) + beta*beta
		     /(2.0*rhat*rhat*r))*Exp_temp*(Exp_temp - 1.0)
		    +(3.0*beta/(2.0*rhat*r2) + 3.0*beta*beta/(2.0*rhat*rhat*r))
		    *Exp_temp*Exp_temp);
	       break;
	    case DT:
	       cerr << "D4phi/Dy3DT Not Coded... " << endl;
	       exit(-1);
	       break;
	    default:
	       cerr << "Error in RadiiMorse::PairPotential -- D3Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
      case D4Y:
	 switch (dt)
	 {
	    case T0:
	       val = At*(beta/(2.0*rhat*r2))
		  *Exp_temp*((15.0/(4.0*r2*r2*r))*(Exp_temp - 1.0)
			     +(15.0*beta/(4.0*rhat*r2*r2))*(2.0*Exp_temp - 1.0)
			     +(6.0*beta*beta/(4.0*rhat*rhat*r2*r))
			     *(4.0*Exp_temp - 1.0)
			     +(beta*beta*beta/(4.0*rhat*rhat*rhat*r2))
			     *(8.0*Exp_temp - 1.0));
	       break;
	    case DT:
	       cerr << "D5phi/Dy4DT Not Coded... " << endl;
	       exit(-1);
	       break;
	    default:
	       cerr << "Error in RadiiMorse::PairPotential -- D4Y,D2T not programmed" << endl;
	       exit(-1);
	 }
	 break;
   }
   return val;
}

void RadiiMorse::Print(ostream &out)
{
   int W=out.width();

   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; B0=" << setw(W) << B0_
       << "; Alpha=" << setw(W) << Alpha_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_;
}

ostream &operator<<(ostream &out,RadiiMorse &A)
{
   A.Print(out);
   return out;
}
