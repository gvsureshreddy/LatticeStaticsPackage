#include "NiTiPressTempPairPotLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTiPressTempPairPotLat::NiTiPressTempPairPotLat(char *datafile)
{
   // First Size Defgrad
   DefGrad_.SetIdentity(DIM3);

   // Get Potential Parameters
   GetParameter("^A0_aa",datafile,"%lf",&A0_aa);
   GetParameter("^B0_aa",datafile,"%lf",&B0_aa);
   GetParameter("^Alpha_aa",datafile,"%lf",&Alpha_aa);
   GetParameter("^Rref_aa",datafile,"%lf",&Rref_aa);
   GetParameter("^Tref_aa",datafile,"%lf",&Tref_aa);
   GetParameter("Tmelt_aa",datafile,"%lf",&Tmelt_aa);

   GetParameter("^A0_bb",datafile,"%lf",&A0_bb);
   GetParameter("^B0_bb",datafile,"%lf",&B0_bb);
   GetParameter("^Alpha_bb",datafile,"%lf",&Alpha_bb);
   GetParameter("^Rref_bb",datafile,"%lf",&Rref_bb);
   GetParameter("^Tref_bb",datafile,"%lf",&Tref_bb);
   GetParameter("Tmelt_bb",datafile,"%lf",&Tmelt_bb);

   GetParameter("^A0_ab",datafile,"%lf",&A0_ab);
   GetParameter("^B0_ab",datafile,"%lf",&B0_ab);
   GetParameter("^Alpha_ab",datafile,"%lf",&Alpha_ab);
   GetParameter("^Rref_ab",datafile,"%lf",&Rref_ab);
   GetParameter("^Tref_ab",datafile,"%lf",&Tref_ab);
   GetParameter("Tmelt_ab",datafile,"%lf",&Tmelt_ab);
   
   // Get Lattice parameters
   GetParameter("^RefLen",datafile,"%lf",&RefLen_);
   GetParameter("^InfluanceDist",datafile,"%u",&InfluanceDist_);
   GetParameter("^Temp",datafile,"%lf",&Temp_);
   GetParameter("^Pressure",datafile,"%lf",&Pressure_);
   GetParameter("^ConvexityDX",datafile,"%lf",&ConvexityDX_);
   
   // needed to initialize reference length
   int iter;
   double DX;
   GetParameter("^MaxIterations",datafile,"%u",&iter);
   GetParameter("^InitializeStepSize",datafile,"%lf",&DX);

   
   int err=0;
   err=FindLatticeSpacing(iter,DX);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << endl;
      exit(-1);
   }
}

int NiTiPressTempPairPotLat::FindLatticeSpacing(int iter,double dx)
{
   double oldPressure=Pressure_,
      oldTemp=Temp_;
   
   Pressure_=0.0;
   Temp_=Tref_aa;
   ShearMod_=1.0;
   DefGrad_.SetIdentity(DIM3);

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

   ShearMod_ = fabs(Stiffness()[5][5]);
   Temp_=oldTemp;
   Pressure_=oldPressure;
   
   return 0;
}

   
// Pair Potential Routines

double NiTiPressTempPairPotLat::Beta(interaction inter,TDeriv dt)
{
   double B0,Alpha,Tref;
   
   switch (inter)
   {
      case aa:
	 B0 = B0_aa;
	 Alpha = Alpha_aa;
	 Tref = Tref_aa;
	 break;
      case bb:
	 B0 = B0_bb;
	 Alpha = Alpha_bb;
	 Tref = Tref_bb;
	 break;
      case ab:
	 B0 = B0_ab;
	 Alpha = Alpha_ab;
	 Tref = Tref_ab;
	 break;
   }
   
   switch (dt)
   {
      case T0:
	 return B0*(1.0+Alpha*(Temp_-Tref)/(Tref));
	 break;
      case DT:
	 return B0*Alpha/Tref;
	 break;
   }

   cerr << "Error in NiTiPressTempPairPotLat::Beta" << endl;
   exit(-1);
}

double NiTiPressTempPairPotLat::Rhat(interaction inter,TDeriv dt)
{
   double num,den,rhat;
   double B0,Tref,Tmelt,Rref;

   switch (inter)
   {
      case aa:
	 B0 = B0_aa;
	 Tref = Tref_aa;
	 Tmelt = Tmelt_aa;
	 Rref = Rref_aa;
	 break;
      case bb:
	 B0 = B0_bb;
	 Tref = Tref_bb;
	 Tmelt = Tmelt_bb;
	 Rref = Rref_bb;
	 break;
      case ab:
	 return (Rhat(aa,dt) + Rhat(bb,dt)) / 2.0;
	 break;
   }

   switch (dt)
   {
      case T0:
	 num = 1.0 - (1.0/(2.0*B0))*
	    log(1.0 - Temp_/(4.0*Tmelt));
	 break;
      case DT:
	 num = -(1.0/(2.0*B0))*
	    (1.0/(1.0 - Temp_/(4.0*Tmelt)))*(-1.0/(4.0*Tmelt));
	 break;
   }

   den = 1.0 - (1.0/(2.0*B0))*
      log(1.0 - Tref/(4.0*Tmelt));

   rhat = Rref*(num/den);

   return rhat;
}

double NiTiPressTempPairPotLat::PairPotential(interaction inter,double r2,
					      YDeriv dy,TDeriv dt)
{
   double beta=Beta(inter),
      rhat=Rhat(inter),
      r = sqrt(r2),
      Exp_temp=exp(-beta*(r/rhat - 1.0));

   double val=0;
   double A0,Tmelt;

   switch (inter)
   {
      case aa:
	 A0 = A0_aa;
	 Tmelt = Tmelt_aa;
	 break;
      case bb:
	 A0 = A0_bb;
	 Tmelt = Tmelt_bb;
	 break;
      case ab:
	 A0 = A0_ab;
	 Tmelt = Tmelt_ab;
	 break;
   }

   switch (dy)
   {
      case Y0:
	 switch (dt)
	 {
	    case T0:
	       val = A0*(1.0 - Temp_/(4.0*Tmelt))*Exp_temp*(Exp_temp - 2.0);
	       break;
	    case DT:
	       val = -A0*(1.0/(4.0*Tmelt))*Exp_temp*(Exp_temp - 2.0)
		  - 2.0*A0*(1.0 - Temp_/(4.0*Tmelt))*(Beta(inter,DT)
						       *(r/rhat -1.0)
							-beta*Rhat(inter,DT)
							*(r/(rhat*rhat)))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	 }
	 break;
      case DY:
	 switch (dt)
	 {
	    case T0:
	       val = -A0*(1.0 - Temp_/(4.0*Tmelt))
		  *(beta/(r*rhat))
		  *Exp_temp*(Exp_temp - 1.0);
	       break;
	    case DT:
	       val = (A0/(4.0*Tmelt))*
		  (beta/(rhat*r))*Exp_temp*(Exp_temp - 1.0)
		  + A0*(1.0 - Temp_/(4.0*Tmelt))*Exp_temp*
		  (((beta*Rhat(inter,DT) - Beta(inter,DT)*rhat)/(rhat*rhat*r))*
		   (Exp_temp - 1.0)
		   +(beta/(rhat*r))*(Beta(inter,DT)*(r/rhat - 1.0)
				     - beta*r*Rhat(inter,DT)/(rhat*rhat))
		   *(2.0*Exp_temp - 1.0));
	       break;
	 }
	 break;
      case D2Y:
	 switch (dt)
	 {
	    case T0:
	       val = A0*(1.0 -Temp_/(4.0*Tmelt))*(beta/(2.0*rhat*r2))
		  *((beta/rhat + 1.0/r)*Exp_temp*(Exp_temp - 1.0)
		    + (beta/rhat)*Exp_temp*Exp_temp);
	       break;
	    case DT:
	       val = (A0/(4.0*Tmelt))*(beta/(2*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0*(1.0 - Temp_/(4.0*Tmelt))*((Beta(inter,DT)*rhat
						    - beta*Rhat(inter,DT))/
						   (2.0*rhat*rhat*r2))
		  *((beta/rhat + 1.0/r)*(Exp_temp*Exp_temp - Exp_temp)
		    +(beta/rhat)*Exp_temp*Exp_temp)
		  +A0*(1.0 - Temp_/(4.0*Tmelt))*(beta/(2.0*rhat*r2))
		  *(((Beta(inter,DT)*rhat - beta*Rhat(inter,DT))/(rhat*rhat))*
		    (2.0*Exp_temp*Exp_temp - Exp_temp)
		    +(beta*Rhat(inter,DT)*r/(rhat*rhat)
		      - Beta(inter,DT)*(r/rhat - 1.0))
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
	       val = -A0*(1.0 - Temp_/(4.0*Tmelt))*(beta/(2.0*rhat*r2))
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
	       val = A0*(1.0 - Temp_/(4.0*Tmelt))*(beta/(2.0*rhat*r2))
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

double NiTiPressTempPairPotLat::PI(const Vector &Dx,const Vector &DX,
				 int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double NiTiPressTempPairPotLat::PSI(const Vector &DX,
				  int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}


double NiTiPressTempPairPotLat::pwr(const double &x,const unsigned y)
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

int NiTiPressTempPairPotLat::IND(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

int NiTiPressTempPairPotLat::IND(int k,int l,int m,int n)
{
   if (k==l)
   {
      if (m==n)
      {
	 return 6*k+m;
      }
      else
      {
	 return 6*k+2+m+n;
      }
   }
   else
   {
      if (m==n)
      {
	 return 6*(2+k+l) + m;
      }
      else
      {
	 return 6*(2+k+l) + 2+m+n;
      }
   }
}

Matrix NiTiPressTempPairPotLat::Phi(unsigned moduliflag,YDeriv dy,TDeriv dt)
{
   Matrix Phi;

   switch (dy)
   {
      case Y0:
	 Phi.Resize(1,1,0.0);
	 break;
      case DY:
	 Phi.Resize(DIM3,DIM3,0.0);
	 break;
      case D2Y:
	 Phi.Resize(6,6,0.0);
	 break;
      case D3Y:
	 Phi.Resize(36,6,0.0);
	 break;
      case D4Y:
	 Phi.Resize(36,36,0.0);
	 break;
   }

   Vector X(DIM3),dummy(DIM3,0.0);
   Vector DX(DIM3),Dx(DIM3);
   Vector ABOffset(DIM3,0.5);
   Matrix Uinv = DefGrad_.Inverse();
   double J=DefGrad_.Det();
   double r2,phi,phi1,phi2,Influancedist[DIM3];
   int k,l,Top[DIM3],Bottom[DIM3],CurrentInfluanceDist;
   interaction Inter;

   for (int i=0;i<DIM3;i++)
      Influancedist[i]=(fabs(Uinv[0][i]) < fabs(Uinv[1][i]) ?
			(fabs(Uinv[1][i]) < fabs(Uinv[2][i]) ?
			fabs(Uinv[2][i]) : fabs(Uinv[1][i])) :
			(fabs(Uinv[0][i]) < fabs(Uinv[2][i]) ?
			 fabs(Uinv[2][i]) : fabs(Uinv[0][i])))*InfluanceDist_;

   X=dummy;
   for (int p=0;p<DIM3;p++)
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
	 for (X[2] = Bottom[2];X[2] <= Top[2];X[2]++)
	 {
	    for (int inter=0;inter<4;inter++)
	    {
	       switch (inter)
	       {
		  case 0:
		     Inter = aa;
		     DX = X*RefLen_;
		     break;
		  case 1:
		     Inter = bb;
		     DX = X*RefLen_;
		     break;
		  case 2:
		     Inter = ab;
		     DX = (X + ABOffset)*RefLen_;
		     break;
		  case 3:
		     Inter = ab;
		     DX = (X - ABOffset)*RefLen_;
		     break;
	       }
		     
	       Dx = DefGrad_ * DX;
	    
	       r2=Dx*Dx;
	       // Only use Sphere of Influance (current)
	       if (r2==0 || r2 > InfluanceDist_*InfluanceDist_)
	       {
		  if (inter == 0) inter++;
		  continue;
	       }

	       phi = PairPotential(Inter,r2,dy,dt);
	       switch (dy)
	       {
		  case Y0:
		     Phi[0][0]+=phi;
		     break;
		  case DY:
		     for (int i=0;i<DIM3;i++)
		     {
			for (int j=i;j<DIM3;j++)
			{
			   Phi[i][j] = Phi[j][i] +=
			      phi*PI(Dx,DX,i,j);
			}
		     }
		     break;
		  case D2Y:
		     phi1=PairPotential(Inter,r2,DY,dt);
		     
		     for (int i=0;i<DIM3;i++)
			for (int j=i;j<DIM3;j++)
			   for (int k=0;k<DIM3;k++)
			      for (int l=k;l<DIM3;l++)
			      {
				 Phi[IND(i,j)][IND(k,l)]+=
				    phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)
				    +phi1*(0.5)*PSI(DX,i,j,k,l);
			      }
		     break;
		  case D3Y:
		     phi1=PairPotential(Inter,r2,D2Y,dt);

		     for (int i=0;i<DIM3;i++)
			for (int j=i;j<DIM3;j++)
			   for (int k=0;k<DIM3;k++)
			      for (int l=k;l<DIM3;l++)
				 for (int m=0;m<DIM3;m++)
				    for (int n=m;n<DIM3;n++)
				    {
				       Phi[IND(i,j,k,l)][IND(m,n)] +=
					  phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PI(Dx,DX,m,n)
					  +0.5*phi1*(PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
						     PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
						     PI(Dx,DX,m,n)*PSI(DX,i,j,k,l));
				    }
		     break;
		  case D4Y:
		     phi1=PairPotential(Inter,r2,D3Y,dt);
		     phi2=PairPotential(Inter,r2,D2Y,dt);

		     for (int i=0;i<DIM3;i++)
			for (int j=i;j<DIM3;j++)
			   for (int k=0;k<DIM3;k++)
			      for (int l=k;l<DIM3;l++)
				 for (int m=0;m<DIM3;m++)
				    for (int n=m;n<DIM3;n++)
				       for (int p=0;p<DIM3;p++)
					  for (int q=p;q<DIM3;q++)
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
      }
   }
   
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*RefLen_*RefLen_*RefLen_*ShearMod_);

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
	    for (int i=0;i<DIM3;i++)
	       for (int j=i;j<DIM3;j++)
		  for (int k=0;k<DIM3;k++)
		     for (int l=k;l<DIM3;l++)
		     {
			for (int q=0;q<DIM3;q++)
			   for (int s=0;s<DIM3;s++)
			   {
			      Phi[IND(i,j)][IND(k,l)] -=
				 (Pressure_/8.0)*(
				    Alt[i][k][q]*Alt[j][l][s] + Alt[j][k][q]*Alt[i][l][s] +
				    Alt[i][l][q]*Alt[j][k][s] + Alt[j][l][q]*Alt[i][k][s] +
				    Alt[i][q][k]*Alt[j][s][l] + Alt[j][q][k]*Alt[i][s][l] +
				    Alt[i][q][l]*Alt[j][s][k] + Alt[j][q][l]*Alt[i][s][k]
				    )*DefGrad_[q][s];
			   }
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


double NiTiPressTempPairPotLat::Energy()
{
   return Phi()[0][0];
}

Matrix NiTiPressTempPairPotLat::Stress()
{
   return Phi(0,DY);
}

Matrix NiTiPressTempPairPotLat::StressDT()
{
   return Phi(0,DY,DT);
}

Matrix NiTiPressTempPairPotLat::Stiffness()
{
   return Phi(0,D2Y);
}

Matrix NiTiPressTempPairPotLat::Moduli()
{
   return Phi(1,D2Y);
}

int NiTiPressTempPairPotLat::StiffnessNulity(double *Min)
{
   int NoNegEigVal = 0;
   int index = 0;

   Matrix EigenValues(1,6);

   EigenValues=SymEigVal(Stiffness());
   if (Min != NULL) *Min = fabs(EigenValues[0][0]);
   for (int i=0;i<6;i++)
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

void NiTiPressTempPairPotLat::Print(ostream &out,PrintDetail flag)
{
   int W=out.width();

   out.width(0);

   double MinEigVal;
   int NoNegEigVal=0;
   
   Matrix
      stiffness = Stiffness(),
      moduli = Moduli(),
      EigenValues(1,6);

   EigenValues=SymEigVal(stiffness);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<6;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   switch (flag)
   {
      case PrintLong:
	 out << "NiTiPressTempPairPotLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	     << "Potential Parameters : "
	     << "A0_aa=  " << setw(W) << A0_aa
	     << "; B0_aa=  " << setw(W) << B0_aa
	     << "; Alpha_aa=" << setw(W) << Alpha_aa << endl
	     << "                       "
	     << "Rref_aa=" << setw(W) << Rref_aa
	     << "; Tref_aa=" << setw(W) << Tref_aa
	     << "; Tmelt_aa=" << setw(W) << Tmelt_aa << endl
	     << "                       "
	     << "A0_bb=  " << setw(W) << A0_bb
	     << "; B0_bb=  " << setw(W) << B0_bb
	     << "; Alpha_bb=" << setw(W) << Alpha_bb << endl
	     << "                       "
	     << "Rref_bb=" << setw(W) << Rref_bb
	     << "; Tref_bb=" << setw(W) << Tref_bb
	     << "; Tmelt_bb=" << setw(W) << Tmelt_bb << endl
	     << "                       "
	     << "A0_ab=  " << setw(W) << A0_ab
	     << "; B0_ab=  " << setw(W) << B0_ab
	     << "; Alpha_ab=" << setw(W) << Alpha_ab << endl
	     << "                       "
	     << "Rref_ab=" << setw(W) << Rref_ab
	     << "; Tref_ab=" << setw(W) << Tref_ab
	     << "; Tmelt_ab=" << setw(W) << Tmelt_ab << endl
	     << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 // passthrough to short
      case PrintShort:
	 out << "Temperature : " << setw(W) << Temp_ << endl
	     << "Pressure (Normalized): " << setw(W) << Pressure_ << endl
	     << "Deformation Gradient:" << setw(W) << DefGrad_
	     << "Potential Value (Normalized):" << setw(W) << Energy() << endl
	     << "Stress (Normalized):" << setw(W) << Stress()
	     << "Stiffness (Normalized):" << setw(W) << stiffness
	     << "Rank 1 Convex:"<< setw(W)
	     << Rank1Convex3D(moduli,ConvexityDX_) << endl
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl;
	 break;
   }
}

ostream &operator<<(ostream &out,NiTiPressTempPairPotLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double NiTiPressTempPairPotLat::Alt[DIM3][DIM3][DIM3]= {0.0, 0.0, 0.0,
							      0.0, 0.0, 1.0,
							      0.0, -1.0,0.0,
							      0.0, 0.0, -1.0,
							      0.0, 0.0, 0.0,
							      1.0, 0.0, 0.0,
							      0.0, 1.0, 0.0,
							      -1.0,0.0, 0.0,
							      0.0, 0.0, 0.0};
