#include "RadiiMorse.h"

RadiiMorse::RadiiMorse(double A0,double B0,double Alpha,double Rref,double Rtheta,
		       double Tref,double Tmelt):
   A0_(A0),B0_(B0),Alpha_(Alpha),Rref_(Rref),Rtheta_(Rtheta),
   Tref_(Tref),Tmelt_(Tmelt)
{
}

double RadiiMorse::Beta(double NTemp,TDeriv dt)
{
   switch (dt)
   {
      case T0:
	 return B0_*(1.0+Alpha_*(NTemp-1.0));
	 break;
      case DT:
	 return B0_*Alpha_;
	 break;
   }

   cerr << "Error in RadiiMorse::Beta" << endl;
   exit(-1);
}

double RadiiMorse::Rhat(double NTemp,TDeriv dt)
{
   double rhat;

   switch (dt)
   {
      case T0:
	 rhat = Rref_*(1.0 + Rtheta_*(NTemp - 1.0));
	 break;
      case DT:
	 rhat = Rref_*Rtheta_;
	 break;
   }

   return rhat;
}

double RadiiMorse::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double beta=Beta(NTemp),
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
	       val = A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0);
	       break;
	    case DT:
	       val = -A0_*(Tref_/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0)
		  - 2.0*A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*(Beta(NTemp,DT)
							      *(r/rhat - 1.0)
							      -beta*Rhat(NTemp,DT)
							      *(r/(rhat*rhat)))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	 }
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))
		  *(beta/(r*rhat))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	    case DT:
	       val = (A0_*Tref_/(4.0*Tmelt_))*
		  (beta/(rhat*r))*Exp_temp*(Exp_temp - 1.0)
		  + A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*Exp_temp*
		  (((beta*Rhat(NTemp,DT) - Beta(NTemp,DT)*rhat)/(rhat*rhat*r))*
		   (Exp_temp - 1.0)
		   +(beta/(rhat*r))*(Beta(NTemp,DT)*(r/rhat - 1.0)
				     - beta*r*Rhat(NTemp,DT)/(rhat*rhat))
		   *(2.0*Exp_temp - 1.0));
	       break;
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *((beta/rhat + 1.0/r)*Exp_temp*(Exp_temp - 1.0)
		    + (beta/rhat)*Exp_temp*Exp_temp);
	       break;
	    case DT:
	       val = (-A0_*Tref_/(4.0*Tmelt_))*(beta/(2*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*((Beta(NTemp,DT)*rhat
							  - beta*Rhat(NTemp,DT))/
							 (2.0*rhat*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *(((Beta(NTemp,DT)*rhat - beta*Rhat(NTemp,DT))/(rhat*rhat))*
		    (2.0*Exp_temp*Exp_temp - Exp_temp)
		    +(beta*Rhat(NTemp,DT)*r/(rhat*rhat)
		      - Beta(NTemp,DT)*(r/rhat - 1.0))
		    *(2.0*beta/rhat*Exp_temp*Exp_temp +
		      (beta/rhat + 1.0/r)
		      *(2.0*Exp_temp*Exp_temp - Exp_temp)));
	       break;
	 }
	 break;
      case D3Y:
	 switch (dt)
	 {
	    case T0:
	       val = -A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *((3.0*beta/(2.0*rhat*r2) + 3.0/(2.0*r2*r) + beta*beta
		     /(2.0*rhat*rhat*r))*Exp_temp*(Exp_temp - 1.0)
		    +(3.0*beta/(2.0*rhat*r2) + 3.0*beta*beta/(2.0*rhat*rhat*r))
		    *Exp_temp*Exp_temp);
	       break;
	    case DT:
	       cerr << "D4phi/Dy3DT Not Coded... " << endl;
	       exit(-1);
	       break;
	 }
	 break;
      case D4Y:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 - NTemp*Tref_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
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
	 }
	 break;
   }
   return val;
}

ostream &operator<<(ostream &out,RadiiMorse &A)
{
   int W=out.width();

   out.width(0);
   
   out << "A0=" << setw(W) << A.A0_
       << "; B0=" << setw(W) << A.B0_
       << "; Alpha=" << setw(W) << A.Alpha_
       << "; Rref=" << setw(W) << A.Rref_
       << "; Rtheta=" << setw(W) << A.Rtheta_
       << "; Tref=" << setw(W) << A.Tref_
       << "; Tmelt=" << setw(W) << A.Tmelt_;

   return out;
}

