#include "TempMorse.h"

TempMorse::TempMorse(double A0,double B0,double Alpha,double Rref,double Tref,
		     double Tmelt):
   A0_(A0),B0_(B0),Alpha_(Alpha),Rref_(Rref),Tref_(Tref),Tmelt_(Tmelt)
{
}

double TempMorse::Beta(double NTemp,TDeriv dt)
{
   double retval;
   switch (dt)
   {
      case T0:
	 retval = B0_*(1.0+Alpha_*(NTemp-1.0));
	 break;
      case DT:
	 retval = B0_*Alpha_;
	 break;
      default:
	 cerr << "Error in TempMorse::Beta" << endl;
	 exit(-1);
   }

   return retval;
}

double TempMorse::Rhat(double NTemp,TDeriv dt)
{
   double num,den,rhat;

   switch (dt)
   {
      case T0:
	 num = 1.0 - (1.0/(2.0*B0_))*
	    log(1.0 - NTemp*Tref_/(4.0*Tmelt_));
	 break;
      case DT:
	 num = -(1.0/(2.0*B0_))*
	    (1.0/(1.0 - NTemp*Tref_/(4.0*Tmelt_)))*(-Tref_/(4.0*Tmelt_));
	 break;
      default:
	 cerr << "Error in TempMorse::Rhat" << endl;
	 exit(-1);
   }

   den = 1.0 - (1.0/(2.0*B0_))*
      log(1.0 - Tref_/(4.0*Tmelt_));

   rhat = Rref_*(num/den);

   return rhat;
}

double TempMorse::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
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
	    default:
	       cerr << "Error in TempMorse::PairPotential -- Y0,D2T not programmed" << endl;
	       exit(-1);
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
	    default:
	       cerr << "Error in TempMorse::PairPotential -- DY,D2T not programmed" << endl;
	       exit(-1);
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
	    default:
	       cerr << "Error in TempMorse::PairPotential -- D2Y,D2T not programmed" << endl;
	       exit(-1);
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
	    default:
	       cerr << "Error in TempMorse::PairPotential -- D3Y,D2T not programmed" << endl;
	       exit(-1);
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
	    default:
	       cerr << "Error in TempMorse::PairPotential -- Dy4,D2T Not Coded... " << endl;
	       exit(-1);
	       break;
	 }
	 break;
      default:
	 cerr << "Error in TempMorse::PairPotential -- Unknown dy" << endl;
	 exit(-1);
   }
   return val;
}

void TempMorse::Print(ostream &out)
{
   int W=out.width();

   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; B0=" << setw(W) << B0_
       << "; Alpha=" << setw(W) << Alpha_
       << "; Rref=" << setw(W) << Rref_
       << "; Tref=" << setw(W) << Tref_
       << "; Tmelt=" << setw(W) << Tmelt_;
}

ostream &operator<<(ostream &out,TempMorse &A)
{
   A.Print(out);
   return out;
}
