#include "SquarePressTempPairPotLat.h"
#include <math.h>

#include "UtilityFunctions.h"

SquarePressTempPairPotLat::SquarePressTempPairPotLat(char *datafile)
{
   // First Size Defgrad
   DefGrad_.SetIdentity(DIM2);

   // Get Potential Parameters
   if(!GetParameter("^A0",datafile,"%lf",&A0_)) exit(-1);
   if(!GetParameter("^B0",datafile,"%lf",&B0_)) exit(-1);
   if(!GetParameter("^Alpha",datafile,"%lf",&Alpha_)) exit(-1);
   if(!GetParameter("^Rref",datafile,"%lf",&Rref_)) exit(-1);
   if(!GetParameter("^Tref",datafile,"%lf",&Rref_)) exit(-1);
   if(!GetParameter("^Tmelt",datafile,"%lf",&Tmelt_)) exit(-1);

   // Get Lattice parameters
   if(!GetParameter("^RefLen",datafile,"%lf",&RefLen_)) exit(-1);
   if(!GetParameter("^InfluanceDist",datafile,"%lf",&InfluanceDist_)) exit(-1);
   if(!GetParameter("^Temp",datafile,"%lf",&Temp_)) exit(-1);
   if(!GetParameter("^Pressure",datafile,"%lf",&Pressure_)) exit(-1);
   if(!GetParameter("^ConvexityDX",datafile,"%lf",&ConvexityDX_)) exit(-1);

   // needed to initialize reference length
   unsigned iter;
   double DX;
   if(!GetParameter("^MaxIterations",datafile,"%u",&iter)) exit(-1);
   if(!GetParameter("^InitializeStepSize",datafile,"%lf",&DX)) exit(-1);

   
   int err=0;
   err=FindLatticeSpacing(iter,DX);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << endl;
      exit(-1);
   }
}

int SquarePressTempPairPotLat::FindLatticeSpacing(int iter,double dx)
{
   double oldPressure=Pressure_,
      oldTemp=Temp_;
   
   Pressure_=0.0;
   Temp_=Tref_;
   ShearMod_=1.0;
   DefGrad_.SetIdentity(DIM2);

   double s=Stress()[0][0];

   cout << setw(20) << Stress() << endl;

   double
      val     = fabs(s),
      sign    = s/val,
      newsign = 0;

   int i=0;

   while ((sign*newsign >= 0) && (i < iter))
   {
      cout << setw(20) << RefLen_ << setw(20) << s << endl;
      RefLen_ -= dx*sign;
      s = Stress()[0][0];
      
      newsign = s/fabs(s);
      i++;
   }

   if (i >= iter)
      return 1;
   else
   {
      i=0;
      while ((fabs(s) > (1.0e-14)*val) && (i < iter))
      {
	 cout << setw(20) << RefLen_ << setw(20) << s << endl;
	 i++;
	 RefLen_ -= dx*newsign/(pwr(2,i));
	 s = Stress()[0][0];
	 newsign=s/fabs(s);
      }
      if (i > iter)
	 return 1;
   }

   ShearMod_ = fabs(Stiffness()[2][2]);
   Temp_=oldTemp;
   Pressure_=oldPressure;
   
   return 0;
}

   
// Pair Potential Routines

double SquarePressTempPairPotLat::Beta(TDeriv dt)
{
   switch (dt)
   {
      case T0:
	 return B0_*(1.0+Alpha_*(Temp_-Tref_)/(Tref_));
	 break;
      case DT:
	 return B0_*Alpha_/Tref_;
	 break;
   }

   cerr << "Error in SquarePressTempPairPotLat::Beta" << endl;
   exit(-1);
}

double SquarePressTempPairPotLat::Rhat(TDeriv dt)
{
   double num,den,rhat;

   switch (dt)
   {
      case T0:
	 num = 1.0 - (1.0/(2.0*B0_))*
	    log(1.0 - Temp_/(4.0*Tmelt_));
	 break;
      case DT:
	 num = -(1.0/(2.0*B0_))*
	    (1.0/(1.0 - Temp_/(4.0*Tmelt_)))*(-1.0/(4.0*Tmelt_));
	 break;
   }

   den = 1.0 - (1.0/(2.0*B0_))*
      log(1.0 - Tref_/(4.0*Tmelt_));

   rhat = Rref_*(num/den);

   return rhat;
}

double SquarePressTempPairPotLat::PairPotential(double r2,YDeriv dy,TDeriv dt)
{
   double beta=Beta(),
      rhat=Rhat(),
      r = sqrt(r2),
      Exp_temp=exp(-beta*(r/rhat - 1.0));

   double val=0;

   switch (dy)
   {
      case Y0:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 - Temp_/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0);
	       break;
	    case DT:
	       val = -A0_*(1.0/(4.0*Tmelt_))*Exp_temp*(Exp_temp - 2.0)
		  - 2.0*A0_*(1.0 - Temp_/(4.0*Tmelt_))*(Beta(DT)*(r/rhat -1.0)
							-beta*Rhat(DT)
							*(r/(rhat*rhat)))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	 }
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -A0_*(1.0 - Temp_/(4.0*Tmelt_))
		  *(beta/(r*rhat))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	    case DT:
	       val = (A0_/(4.0*Tmelt_))*
		  (beta/(rhat*r))*Exp_temp*(Exp_temp - 1.0)
		  + A0_*(1.0 - Temp_/(4.0*Tmelt_))*Exp_temp*
		  (((beta*Rhat(DT) - Beta(DT)*rhat)/(rhat*rhat*r))*
		   (Exp_temp - 1.0)
		   +(beta/(rhat*r))*(Beta(DT)*(r/rhat - 1.0)
				     - beta*r*Rhat(DT)/(rhat*rhat))
		   *(2.0*Exp_temp - 1.0));
	       break;
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = A0_*(1.0 -Temp_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *((beta/rhat + 1.0/r)*Exp_temp*(Exp_temp - 1.0)
		    + (beta/rhat)*Exp_temp*Exp_temp);
	       break;
	    case DT:
	       val = (-A0_/(4.0*Tmelt_))*(beta/(2*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - Temp_/(4.0*Tmelt_))*((Beta(DT)*rhat
						    - beta*Rhat(DT))/
						   (2.0*rhat*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0_*(1.0 - Temp_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
		  *(((Beta(DT)*rhat - beta*Rhat(DT))/(rhat*rhat))*
		    (2.0*Exp_temp*Exp_temp - Exp_temp)
		    +(beta*Rhat(DT)*r/(rhat*rhat) - Beta(DT)*(r/rhat - 1.0))
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
	       val = -A0_*(1.0 - Temp_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
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
	       val = A0_*(1.0 - Temp_/(4.0*Tmelt_))*(beta/(2.0*rhat*r2))
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
	       
// End Pair Potential Routines

// Lattice Routines

double SquarePressTempPairPotLat::PI(const Vector &Dx,const Vector &DX,
				 int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double SquarePressTempPairPotLat::PSI(const Vector &DX,
				  int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}


double SquarePressTempPairPotLat::pwr(const double &x,const unsigned y)
{
   if (y==1)
   {
      return x;
   }
   else
   {
      return x*pwr(x,y-1);
   }
}

int SquarePressTempPairPotLat::IND(int i,int j)
{
   if (i==j)
      return i;
   else
      return 1+i+j;
}

int SquarePressTempPairPotLat::IND(int k,int l,int m,int n)
{
   if (k==l)
   {
      if (m==n)
      {
	 return 3*k+m;
      }
      else
      {
	 return 3*k+1+m+n;
      }
   }
   else
   {
      if (m==n)
      {
	 return 3*(1+k+l) + m;
      }
      else
      {
	 return 3*(1+k+l) + 1+m+n;
      }
   }
}

Matrix SquarePressTempPairPotLat::Phi(unsigned moduliflag,YDeriv dy,TDeriv dt)
{
   Matrix Phi;

   switch (dy)
   {
      case Y0:
	 Phi.Resize(1,1,0.0);
	 break;
      case DY:
	 Phi.Resize(DIM2,DIM2,0.0);
	 break;
      case D2Y:
	 Phi.Resize(3,3,0.0);
	 break;
      case D3Y:
	 Phi.Resize(9,3,0.0);
	 break;
      case D4Y:
	 Phi.Resize(9,9,0.0);
	 break;
   }

   Vector X(DIM2),dummy(DIM2,0.0);
   Vector DX(DIM2),Dx(DIM2);
   Matrix Uinv = DefGrad_.Inverse();
   double J=DefGrad_.Det();
   double r2,phi,phi1,phi2,Influancedist[DIM2];
   int k,l,Top[DIM2],Bottom[DIM2],CurrentInfluanceDist;

   for (int i=0;i<DIM2;i++)
      Influancedist[i]=(fabs(Uinv[0][i]) < fabs(Uinv[1][i]) ?
			fabs(Uinv[1][i]) : fabs(Uinv[0][i]))*InfluanceDist_;

   X=dummy;
   for (int p=0;p<DIM2;p++)
   {
      // set influance distance based on cube size
      //
      // alos setup to be large enough to encompass Eulerian sphere
      CurrentInfluanceDist = int(ceil(Influancedist[p]));

      Top[p] = CurrentInfluanceDist;
      Bottom[p] = -CurrentInfluanceDist;
   }

   for (X[0] = Bottom[0];X[0] <= Top[0];X[0]++)
   {
      for (X[1] = Bottom[1];X[1] <= Top[1];X[1]++)
      {
	 DX = X*RefLen_;
	 Dx = DefGrad_ * DX;

	 r2=Dx*Dx;
	 // Only use Sphere of Influance (current)
	 if (r2==0 || r2 > InfluanceDist_*InfluanceDist_)
	 {
	    continue;
	 }

	 phi = PairPotential(r2,dy,dt);
	 switch (dy)
	 {
	    case Y0:
	       Phi[0][0]+=phi;
	       break;
	    case DY:
	       for (int i=0;i<DIM2;i++)
	       {
		  for (int j=i;j<DIM2;j++)
		  {
		     Phi[i][j] = Phi[j][i] +=
			phi*PI(Dx,DX,i,j);
		  }
	       }
	       break;
	    case D2Y:
	       phi1=PairPotential(r2,DY,dt);

	       for (int i=0;i<DIM2;i++)
		  for (int j=i;j<DIM2;j++)
		     for (int k=0;k<DIM2;k++)
			for (int l=k;l<DIM2;l++)
			{
			   Phi[IND(i,j)][IND(k,l)]+=
			      phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)
			      +phi1*(0.5)*PSI(DX,i,j,k,l);
			}
	       break;
	    case D3Y:
	       phi1=PairPotential(r2,D2Y,dt);

	       for (int i=0;i<DIM2;i++)
		  for (int j=i;j<DIM2;j++)
		     for (int k=0;k<DIM2;k++)
			for (int l=k;l<DIM2;l++)
			   for (int m=0;m<DIM2;m++)
			      for (int n=m;n<DIM2;n++)
			      {
				 Phi[IND(i,j,k,l)][IND(m,n)] +=
				    phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PI(Dx,DX,m,n)
				    +0.5*phi1*(PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
					       PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
					       PI(Dx,DX,m,n)*PSI(DX,i,j,k,l));
			      }
	       break;
	    case D4Y:
	       phi1=PairPotential(r2,D3Y,dt);
	       phi2=PairPotential(r2,D2Y,dt);

	       for (int i=0;i<DIM2;i++)
		  for (int j=i;j<DIM2;j++)
		     for (int k=0;k<DIM2;k++)
			for (int l=k;l<DIM2;l++)
			   for (int m=0;m<DIM2;m++)
			      for (int n=m;n<DIM2;n++)
				 for (int p=0;p<DIM2;p++)
				    for (int q=p;q<DIM2;q++)
				    {
				       Phi[IND(i,j,k,l)][IND(m,n,p,q)]+=
					  phi*(PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*
					       PI(Dx,DX,m,n)*PI(Dx,DX,p,q)) +
					  0.5*phi1*PI(Dx,DX,p,q)*(
					     PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
					     PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
					     PI(Dx,DX,m,n)*PSI(DX,i,j,k,l)) +
					  0.25*phi2*(
					     PSI(DX,i,j,m,n)*PSI(DX,k,l,p,q) +
					     PSI(DX,k,l,m,n)*PSI(DX,i,j,p,q) +
					     PSI(DX,i,j,k,l)*PSI(DX,m,n,p,q));
				    }
	       break;
	 }
      }
   }

   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*RefLen_*RefLen_*ShearMod_);

   if (moduliflag)
   {
      return Phi;
   }
   
   // Add in External Work Terms
   // Recall that Pressure_ is applied stress
   // Thus E = W - pJ

   if (dt == T0)
   {
      switch (dy)
      {
	 case Y0:
	    Phi[0][0] = Phi[0][0] - Pressure_*J;
	    break;
	 case DY:
	    Phi = Phi - Pressure_*J*Uinv;
	    break;
	 case D2Y:
	    for (int i=0;i<DIM2;i++)
	       for (int j=i;j<DIM2;j++)
		  for (int k=0;k<DIM2;k++)
		     for (int l=k;l<DIM2;l++)
		     {
			Phi[IND(i,j)][IND(k,l)] -=
			   (Pressure_/8.0)*(
			      Alt[i][k]*Alt[j][l] +
			      Alt[i][l]*Alt[j][k] +
			      Alt[j][k]*Alt[i][l] +
			      Alt[j][l]*Alt[i][k] +
			      Alt[k][i]*Alt[l][j] +
			      Alt[k][j]*Alt[l][i] +
			      Alt[l][i]*Alt[k][j] +
			      Alt[l][j]*Alt[k][i]);
		     }
	    break;
	 case D3Y:
	    cerr << "External work terms not programmed for D3Y" << endl;
	    exit(-1);
	    break;
	 case D4Y:
	    cerr << "External work terms not programmed for D4Y" << endl;
	    exit(-2);
	    break;
      }
   }

   return Phi;
}

int SquarePressTempPairPotLat::StiffnessNulity(double *Min)
{
   int NoNegEigVal = 0;
   int index = 0;

   Matrix EigenValues(1,3);

   EigenValues=SymEigVal(Stiffness());
   if (Min != NULL) *Min = fabs(EigenValues[0][0]);
   for (int i=0;i<3;i++)
   {
      if (EigenValues[0][i] < 0.0) NoNegEigVal++;
      if ((Min != NULL)
	  && (fabs(EigenValues[0][i]) < *Min))
      {
	 *Min = fabs(EigenValues[0][i]);
	 index = i;
      }
   }

   if (Min != NULL) *Min = EigenValues[0][index];
   return NoNegEigVal;
}

void SquarePressTempPairPotLat::Print(ostream &out,PrintDetail flag)
{
   int W=out.width();

   out.width(0);
   cout.width(0);

   double MinEigVal;
   int NoNegEigVal=0;

   double energy = Energy();
   Matrix
      stress = Stress(),
      stiffness = Stiffness(),
      moduli = Moduli(),
      EigenValues(1,3);

   EigenValues=SymEigVal(stiffness);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<3;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   int Rank1Convex = FullScanRank1Convex2D(moduli,ConvexityDX_);

   switch (flag)
   {
      case PrintLong:
	 out << "SquarePressTempPairPotLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	     << "Potential Parameters : "
	     << "A0=  " << setw(W) << A0_
	     << "; B0=  " << setw(W) << B0_
	     << "; Alpha=" << setw(W) << Alpha_ << endl
	     << "                       "
	     << "Rref=" << setw(W) << Rref_
	     << "; Tref=" << setw(W) << Tref_
	     << "; Tmelt=" << setw(W) << Tmelt_ << endl
	     << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 cout << "SquarePressTempPairPotLat:" << endl << endl
	      << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	      << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	      << "Potential Parameters : "
	      << "A0=  " << setw(W) << A0_
	      << "; B0=  " << setw(W) << B0_
	      << "; Alpha=" << setw(W) << Alpha_ << endl
	      << "                       "
	      << "Rref=" << setw(W) << Rref_
	      << "; Tref=" << setw(W) << Tref_
	      << "; Tmelt=" << setw(W) << Tmelt_ << endl
	      << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 // passthrough to short
      case PrintShort:
	 out << "Temperature : " << setw(W) << Temp_ << endl
	     << "Pressure (Normalized): " << setw(W) << Pressure_ << endl
	     << "Deformation Gradient:" << setw(W) << DefGrad_
	     << "Potential Value (Normalized):" << setw(W) << energy << endl
	     << "Stress (Normalized):" << setw(W) << stress
	     << "Stiffness (Normalized):" << setw(W) << stiffness
	     << "Rank 1 Convex:" << setw(W) << Rank1Convex << endl
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl;
	 cout << "Temperature : " << setw(W) << Temp_ << endl
	      << "Pressure (Normalized): " << setw(W) << Pressure_ << endl
	      << "Deformation Gradient:" << setw(W) << DefGrad_
	      << "Potential Value (Normalized):" << setw(W) << energy << endl
	      << "Stress (Normalized):" << setw(W) << stress
	      << "Stiffness (Normalized):" << setw(W) << stiffness
	      << "Rank 1 Convex:" << setw(W) << Rank1Convex << endl
	      << "Eigenvalue Info:"  << setw(W) << EigenValues
	      << "Bifurcation Info:" << setw(W) << MinEigVal
	      << setw(W) << NoNegEigVal << endl;
	 break;
   }
}

ostream &operator<<(ostream &out,SquarePressTempPairPotLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double SquarePressTempPairPotLat::Alt[DIM2][DIM2] = {0.0, 1.0, -1.0, 0.0};
