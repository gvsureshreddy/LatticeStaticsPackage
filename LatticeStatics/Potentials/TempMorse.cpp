#include "TempMorse.h"

TempMorse::TempMorse(double A0,double B0,double Alpha,double Rref,double Tref,
		     double Tmelt):
   A0_(A0),B0_(B0),Alpha_(Alpha),Rref_(Rref),Tref_(Tref),Tmelt_(Tmelt)
{
}

double TempMorse::Beta(double Temp,TDeriv dt)
{
   switch (dt)
   {
      case T0:
	 return B0_*(1.0+Alpha_*(Temp-Tref_)/(Tref_));
	 break;
      case DT:
	 return B0_*Alpha_/Tref_;
	 break;
   }

   cerr << "Error in TempMorse::Beta" << endl;
   exit(-1);
}

double TempMorse::Rhat(double Temp,TDeriv dt)
{
   double num,den,rhat;

   switch (dt)
   {
      case T0:
	 num = 1.0 - (1.0/(2.0*B0_))*
	    log(1.0 - Temp/(4.0*Tmelt_));
	 break;
      case DT:
	 num = -(1.0/(2.0*B0_))*
	    (1.0/(1.0 - Temp/(4.0*Tmelt_)))*(-1.0/(4.0*Tmelt_));
	 break;
   }

   den = 1.0 - (1.0/(2.0*B0_))*
      log(1.0 - Tref_/(4.0*Tmelt_));

   rhat = Rref_*(num/den);

   return rhat;
}

double TempMorse::PairPotential(double Temp,double r2,YDeriv dy,TDeriv dt)
{
   double beta=Beta(Temp),
      rhat=Rhat(Temp),
      r = sqrt(r2),
      Exp_temp=exp(-beta*(r/rhat - 1.0));

   double val=0;

   switch (dy)
   {
      case Y0:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 - Temp/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0);
	       break;
	    case DT:
	       val = -A0_*(1.0/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0)
		  - 2.0*A0_*(1.0 - Temp/(4.0*Tmelt_))*(Beta(Temp,DT)
						       *(r/rhat -1.0)
							-beta*Rhat(Temp,DT)
							*(r/(rhat*rhat)))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	 }
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -A0_*(1.0 - Temp/(4.0*Tmelt_))
		  *(beta/(r*rhat))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	    case DT:
	       val = (A0_/(4.0*Tmelt_))*
		  (beta/(rhat*r))*Exp_temp*(Exp_temp - 1.0)
		  + A0_*(1.0 - Temp/(4.0*Tmelt_))*Exp_temp*
		  (((beta*Rhat(Temp,DT) - Beta(Temp,DT)*rhat)/(rhat*rhat*r))*
		   (Exp_temp - 1.0)
		   +(beta/(rhat*r))*(Beta(Temp,DT)*(r/rhat - 1.0)
				     - beta*r*Rhat(Temp,DT)/(rhat*rhat))
		   *(2.0*Exp_temp - 1.0));
	       break;
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 -Temp/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *((beta/rhat + 1.0/r)*Exp_temp*(Exp_temp - 1.0)
		    + (beta/rhat)*Exp_temp*Exp_temp);
	       break;
	    case DT:
	       val = (A0_/(4.0*Tmelt_))*(beta/(2*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - Temp/(4.0*Tmelt_))*((Beta(Temp,DT)*rhat
						    - beta*Rhat(Temp,DT))/
						   (2.0*rhat*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - Temp/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *(((Beta(Temp,DT)*rhat - beta*Rhat(Temp,DT))/(rhat*rhat))*
		    (2.0*Exp_temp*Exp_temp - Exp_temp)
		    +(beta*Rhat(Temp,DT)*r/(rhat*rhat)
		      - Beta(Temp,DT)*(r/rhat - 1.0))
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
	       val = -A0_*(1.0 - Temp/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
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
	       val = A0_*(1.0 - Temp/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
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

ostream &operator<<(ostream &out,TempMorse &A)
{
   int W=out.width();

   out.width(0);
   
   out << "Potential Parameters : " << "\n\t"
       << "A0=" << setw(W) << A.A0_
       << "; B0=" << setw(W) << A.B0_
       << "; Alpha=" << setw(W) << A.Alpha_
       << "; Rref=" << setw(W) << A.Rref_
       << "; Tref=" << setw(W) << A.Tref_
       << "; Tmelt=" << setw(W) << A.Tmelt_ << endl;

   return out;
}

