#include "NiTi9TPPLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTi9TPPLat::NiTi9TPPLat(char *datafile)
{
   // First Size DOF
   DOF_.Resize(DOFS);
   // Set LatticeVec_
   LatticeVec_.Resize(DIM3,DIM3,0.0);
   LatticeVec_[0][0] = 1.0;
   LatticeVec_[1][1] = 1.0;
   LatticeVec_[2][2] = 1.0;
   
   // Setup Bodyforce_
   BodyForce_[0].Resize(DIM3,0.0);
   BodyForce_[1].Resize(DIM3,0.0);

   // Get Potential Parameters
   double Tref,A0,B0,Alpha,Rref,Rtheta,Tmelt;
   GetParameter("^Tref",datafile,"%lf",&Tref);
   GetParameter("^A0_aa",datafile,"%lf",&A0);
   GetParameter("^B0_aa",datafile,"%lf",&B0);
   GetParameter("^Alpha_aa",datafile,"%lf",&Alpha);
   GetParameter("^Rref_aa",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_aa",datafile,"%lf",&Rtheta);
   GetParameter("^Tmelt_aa",datafile,"%lf",&Tmelt);
   Potential_[aa]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);

   GetParameter("^A0_bb",datafile,"%lf",&A0);
   GetParameter("^B0_bb",datafile,"%lf",&B0);
   GetParameter("^Alpha_bb",datafile,"%lf",&Alpha);
   GetParameter("^Rref_bb",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_bb",datafile,"%lf",&Rtheta);
   GetParameter("^Tmelt_bb",datafile,"%lf",&Tmelt);
   Potential_[bb]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);

   GetParameter("^A0_ab",datafile,"%lf",&A0);
   GetParameter("^B0_ab",datafile,"%lf",&B0);
   GetParameter("^Alpha_ab",datafile,"%lf",&Alpha);
   GetParameter("^Rref_ab",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_ab",datafile,"%lf",&Rtheta);
   GetParameter("^Tmelt_ab",datafile,"%lf",&Tmelt);
   Potential_[ab]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);
   
   // Get Lattice parameters
   GetParameter("^RefLen",datafile,"%lf",&RefLen_);
   GetParameter("^InfluanceDist",datafile,"%u",&InfluanceDist_);
   GetParameter("^NTemp",datafile,"%lf",&NTemp_);
   GetParameter("^Pressure",datafile,"%lf",&Pressure_);
   GetParameter("^ConvexityDX",datafile,"%lf",&ConvexityDX_);
   char temp[30];
   for (int i=0;i<INTERNAL_ATOMS;i++)
   {
      sprintf(temp,"^AtomicMass_%u",i);
      GetParameter(temp,datafile,"%lf",&(AtomicMass_[i]));
   }
   
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

int NiTi9TPPLat::FindLatticeSpacing(int iter,double dx)
{
   double oldPressure=Pressure_,
      oldTemp=NTemp_;

   Pressure_=0.0;
   NTemp_=1.0;
   ShearMod_=1.0;
   DOF_[0] = DOF_[1] = DOF_[2] = 1.0;
   for (int i=3;i<DOFS;i++)
   {
      DOF_[i] = 0.0;
   }

   double s=Stress()[0][0];

   cout << setw(20) << Stress() << endl;

   double
      val     = fabs(s),
      sign    = s/val,
      newsign = 0;

   int i=0;

   if (!(val < 1.0e-14))
   {
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
	 while ((fabs(s) > 1.0e-13) && (i < iter))
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
   }

   ShearMod_ = 0.25*fabs(Stiffness()[5][5]);
   NTemp_=oldTemp;
   Pressure_=oldPressure;

   return 0;
}   

// Lattice Routines

double NiTi9TPPLat::PI(const Vector &Dx,const Vector &DX,
			int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double NiTi9TPPLat::PSI(const Vector &DX,
			 int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}

double NiTi9TPPLat::OMEGA(const Vector &Dx,int p,int q,int i, int j)
{
   double ret=0;
   
   ret=0;
   for (int s=0;s<DIM3;s++)
   {
      for (int t=0;t<DIM3;t++)
      {
	 ret += (LatticeVec_[j][s]*DOF_[INDU(s,t)]*Dx[t] +
		 Dx[t]*DOF_[INDU(t,s)]*LatticeVec_[j][s]);
      }
   }
   ret *= DELTA(i,p,q);

   return ret;
}

double NiTi9TPPLat::SIGMA(int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   for (int s=0;s<DIM3;s++)
   {
      for (int t=0;t<DIM3;t++)
      {
	 for (int r=0;r<DIM3;r++)
	 {
	    tmp += (LatticeVec_[j][s]*DOF_[INDU(s,t)]*DOF_[INDU(t,r)]*LatticeVec_[l][r] +
		    LatticeVec_[l][s]*DOF_[INDU(s,t)]*DOF_[INDU(t,r)]*LatticeVec_[j][r]);
	 }
      }
   }
   
   return DELTA(i,p,q)*DELTA(k,p,q)*tmp;
}

double NiTi9TPPLat::GAMMA(const Vector &Dx,const Vector &DX,
			   int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   
   for (int s=0;s<DIM3;s++)
   {
      tmp += (LatticeVec_[l][s]*DOF_[INDU(s,i)]*DX[j] +
	      LatticeVec_[l][s]*DOF_[INDU(s,j)]*DX[i] +
	      DX[i]*DOF_[INDU(j,s)]*LatticeVec_[l][s] +
	      DX[j]*DOF_[INDU(i,s)]*LatticeVec_[l][s]);
   }
   
   
   return (0.5*DELTA(k,p,q)*(LatticeVec_[l][i]*Dx[j] + LatticeVec_[l][j]*Dx[i] + tmp));
}

double NiTi9TPPLat::THETA(const Vector &DX,int p,int q,int i,int j,int k,int l,
			   int m, int n)
{
   return (0.5*DELTA(m,p,q)*(Del(i,k)*LatticeVec_[n][j]*DX[l]
			     + Del(i,k)*DX[j]*LatticeVec_[n][l]
			     + Del(i,l)*LatticeVec_[n][j]*DX[k]
			     + Del(i,l)*DX[j]*LatticeVec_[n][k]
			     + Del(j,k)*LatticeVec_[n][i]*DX[l]
			     + Del(j,k)*DX[i]*LatticeVec_[n][l]
			     + Del(j,l)*LatticeVec_[n][i]*DX[k]
			     + Del(j,l)*DX[i]*LatticeVec_[n][k]));
}

double NiTi9TPPLat::XI(int p,int q,int i,int j,int k,int l,int m,int n)
{
   double tmp=0;

   for (int s=0;s<DIM3;s++)
   {
      tmp += (LatticeVec_[j][m]*DOF_[INDU(n,s)]*LatticeVec_[l][s]
	      + LatticeVec_[j][n]*DOF_[INDU(m,s)]*LatticeVec_[l][s]
	      + LatticeVec_[j][s]*DOF_[INDU(s,m)]*LatticeVec_[l][n]
	      + LatticeVec_[j][s]*DOF_[INDU(s,n)]*LatticeVec_[l][m]
	      + LatticeVec_[l][m]*DOF_[INDU(n,s)]*LatticeVec_[j][s]
	      + LatticeVec_[l][n]*DOF_[INDU(m,s)]*LatticeVec_[j][s]
	      + LatticeVec_[l][s]*DOF_[INDU(s,m)]*LatticeVec_[j][n]
	      + LatticeVec_[l][s]*DOF_[INDU(s,n)]*LatticeVec_[j][m]);
   }

   return (0.5*DELTA(i,p,q)*DELTA(k,p,q)*tmp);
}

double NiTi9TPPLat::LAMDA(int p,int q,int i,int j,int k,int l,int m,int n,int a,
			   int b)
{
   return (0.5*DELTA(m,p,q)*DELTA(a,p,q)*(Del(i,k)*LatticeVec_[n][j]*LatticeVec_[b][l]
					  + Del(i,k)*LatticeVec_[b][j]*LatticeVec_[n][l]
					  + Del(i,l)*LatticeVec_[n][j]*LatticeVec_[b][k]
					  + Del(i,l)*LatticeVec_[b][j]*LatticeVec_[n][k]
					  + Del(j,k)*LatticeVec_[n][i]*LatticeVec_[b][l]
					  + Del(j,k)*LatticeVec_[b][i]*LatticeVec_[n][l]
					  + Del(j,l)*LatticeVec_[n][i]*LatticeVec_[b][k]
					  + Del(j,l)*LatticeVec_[b][i]*LatticeVec_[n][k]));
}

double NiTi9TPPLat::pwr(const double &x,const unsigned y)
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

inline int NiTi9TPPLat::INDU(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

inline int NiTi9TPPLat::INDUU(int k,int l,int m,int n)
{
   if (k==l)
   {
      if (m==n)
      {
	 return DOFS*k+m;
      }
      else
      {
	 return DOFS*k + 2+m+n;
      }
   }
   else
   {
      if (m==n)
      {
	 return DOFS*(2+k+l) + m;
      }
      else
      {
	 return DOFS*(2+k+l) + 2+m+n;
      }
   }
}

inline int NiTi9TPPLat::INDV(int i,int j)
{
   if (!i)
   {
      cerr << "Error: INDV(i,j) i==0!!!!!" << endl;
      exit(-1);
   }

   return 6 + (i-1)*3 + j;
}

inline int NiTi9TPPLat::INDVV(int k,int l,int m,int n)
{
   if (!k || !m)
   {
      cerr << "Error : INDVV(k,l,m,n) i==0 OR m==0!!!!!!" << endl;
      exit(-1);
   }

   return 6*DOFS + 6
      + DOFS*( (k-1)*3 + l )
      + ( (m-1)*3 + n );
}

inline int NiTi9TPPLat::INDUV(int i,int j,int m,int n)
{
   if (!m)
   {
      cerr << "Error : INDUV(i,j,m,n) m==0!!!!!!" << endl;
      exit(-1);
   }

   if (i==j)
   {
      return DOFS*i + 6 + (m-1)*3+n;
   }
   else
   {
      return DOFS*(2+i+j) + 6 + (m-1)*3+n;
   }
}

inline int NiTi9TPPLat::INDVU(int m,int n,int i,int j)
{
   if (!m)
   {
      cerr << "Error : INDVU(m,n,i,j) m==0!!!!!!" << endl;
      exit(-1);
   }

   if (i==j)
   {
      return 6*DOFS + DOFS*( (m-1)*3+n ) + i;
   }
   else
   {
      return 6*DOFS + DOFS*( (m-1)*3+n ) + 2+i+j;
   }
}

Matrix NiTi9TPPLat::Phi(unsigned moduliflag,PairPotentials::YDeriv dy,
			 PairPotentials::TDeriv dt)
{
   static Matrix U(DIM3,DIM3);
   static Matrix V(INTERNAL_ATOMS,3);
   static Matrix Eigvals(1,DIM3);
   static double X[DIM3];
   static Vector DX(DIM3),Dx(DIM3);
   static Vector Direction(DIM3);
   static double ForceNorm;
   static double J;
   static int p,q;
   static int i,j,k,l,m,n,s,t;
   static double r2,phi,phi1,phi2,phi3,Influancedist[DIM3],tmp;
   static int Top[DIM3],Bottom[DIM3],CurrentInfluanceDist;
   static interaction Inter;
   static Matrix Phi;

   switch (dy)
   {
      case PairPotentials::Y0:
	 Phi.Resize(1,1,0.0);
	 break;
      case PairPotentials::DY:
	 Phi.Resize(1,DOFS,0.0);
	 break;
      case PairPotentials::D2Y:
	 Phi.Resize(DOFS,DOFS,0.0);
	 break;
      case PairPotentials::D3Y:
	 Phi.Resize(DOFS*DOFS,DOFS,0.0);
	 break;
      case PairPotentials::D4Y:
	 Phi.Resize(DOFS*DOFS,DOFS*DOFS,0.0);
	 break;
   }


   U[0][0] = DOF_[0];
   U[1][1] = DOF_[1];
   U[2][2] = DOF_[2];
   U[0][1] = U[1][0] = DOF_[3];
   U[0][2] = U[2][0] = DOF_[4];
   U[1][2] = U[2][1] = DOF_[5];

   V[0][0] = 0.0;
   V[0][1] = 0.0;
   V[0][2] = 0.0;
   V[1][0] = DOF_[6];
   V[1][1] = DOF_[7];
   V[1][2] = DOF_[8];

   // find largest eigenvalue of the inverse transformation
   // (i.e. from current to ref) and use influence cube of
   // that size...
   //
   // Use the fact that eigs of Uinv = 1/ eigs of U.
   J = U.Det();
   Eigvals = SymEigVal(U);
   tmp = Eigvals[0][0];
   for (i=0;i<DIM3;i++)
      if (Eigvals[0][i] < tmp) tmp = Eigvals[0][i];
   
   // Set to inverse eigenvalue
   tmp = 1.0/tmp;
   for (i=0;i<DIM3;i++)
   {
      Influancedist[i]=tmp*InfluanceDist_;
   }
   
   for (p=0;p<DIM3;p++)
   {
      // set influance distance based on cube size
      //
      // also setup to be large enough to encompass Eulerian sphere
      CurrentInfluanceDist = int(ceil(Influancedist[p]));

      Top[p] = CurrentInfluanceDist;
      Bottom[p] = -CurrentInfluanceDist;

      // misc initialization
      for (t=0;t<INTERNAL_ATOMS;t++)
      {
	 BodyForce_[t][p] = 0.0;
      }
   }
   // misc initialization
   ForceNorm = 0.0;

   for (p=0;p<INTERNAL_ATOMS;p++)
   {
      for (q=0;q<INTERNAL_ATOMS;q++)
      {
	 Inter = INTER[p][q];
	 
	 for (X[0] = Bottom[0];X[0] <= Top[0];X[0]++)
	 {
	    for (X[1] = Bottom[1];X[1] <= Top[1];X[1]++)
	    {
	       for (X[2] = Bottom[2];X[2] <= Top[2];X[2]++)
	       {
		  for (i=0;i<DIM3;i++)
		  {
		     DX[i] = 0.0;
				     
		     for (j=0;j<DIM3;j++)
		     {
			DX[i] += (X[j] + A[q][j] - A[p][j] + V[q][j] - V[p][j])*
			   LatticeVec_[j][i]*RefLen_;

		     }
		  }

		  r2 = 0.0;
		  for (i=0;i<DIM3;i++)
		  {
		     Dx[i] = 0.0;

		     for (j=0;j<DIM3;j++)
		     {
			Dx[i] += U[i][j] * DX[j];
		     }
		     r2 += Dx[i]*Dx[i];
		  }
		  // Only use Sphere of Influance (current)
		  if (r2==0 || r2 > InfluanceDist_*InfluanceDist_)
		  {
		     continue;
		  }

		  // Calculate bodyforce
		  // NOTE: phi1 = d(phi)/d(r2)
		  // We need d(phi)/dr = 2*r*d(phi)/d(r2)
		  phi1 = 2.0*sqrt(r2)*
		     Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,
						     PairPotentials::T0);
		  if (ForceNorm < fabs(-phi1/2.0))
		  {
		     ForceNorm = fabs(-phi1/2.0);
		  }
		  for (i=0;i<DIM3;i++)
		  {
		     BodyForce_[p][i] += -phi1*Dx[i]/(2.0*sqrt(r2));
		  }
		  
		  // Calculate Phi
		  phi = Potential_[Inter].PairPotential(NTemp_,r2,dy,dt);
		  switch (dy)
		  {
		     case PairPotentials::Y0:
			Phi[0][0]+=phi;
			break;
		     case PairPotentials::DY:
			for (i=0;i<DIM3;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      Phi[0][INDU(i,j)] += phi*PI(Dx,DX,i,j);
			   }
			}
			for (i=1;i<INTERNAL_ATOMS;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      Phi[0][INDV(i,j)] += phi*OMEGA(Dx,p,q,i,j);
			   }
			}
			break;
		     case PairPotentials::D2Y:
			phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,
							     dt);
			
			//Upper Diag Block (6,6)
			for (i=0;i<DIM3;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      for (k=0;k<DIM3;k++)
			      {
				 for (l=0;l<DIM3;l++)
				 {
				    Phi[INDU(i,j)][INDU(k,l)]+=
				       phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)
				       +phi1*(0.5)*PSI(DX,i,j,k,l);
				 }
			      }
			   }
			}
			//Lower Diag Block (3,3)
			for (i=1;i<INTERNAL_ATOMS;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      for (k=1;k<INTERNAL_ATOMS;k++)
			      {
				 for (l=0;l<DIM3;l++)
				 {
				    Phi[INDV(i,j)][INDV(k,l)]+=
				       phi*OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,k,l)
				       +phi1*SIGMA(p,q,i,j,k,l);
				 }
			      }
			   }
			}
			//Off Diag Blocks
			for (i=0;i<DIM3;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      for (k=1;k<INTERNAL_ATOMS;k++)
			      {
				 for (l=0;l<DIM3;l++)
				 {
				    Phi[INDU(i,j)][INDV(k,l)] =
				       Phi[INDV(k,l)][INDU(i,j)] +=
				       phi*PI(Dx,DX,i,j)*OMEGA(Dx,p,q,k,l)
				       +phi1*GAMMA(Dx,DX,p,q,i,j,k,l);
				 }
			      }
			   }
			}		
			break;
		     case PairPotentials::D3Y:
			phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y,dt);
			phi2=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,dt);

			// DU^3 block
			for (i=0;i<DIM3;i++)
			   for (j=0;j<DIM3;j++)
			      for (k=0;k<DIM3;k++)
				 for (l=0;l<DIM3;l++)
				    for (m=0;m<DIM3;m++)
				       for (n=0;n<DIM3;n++)
				       {
					  Phi[INDUU(i,j,k,l)][INDU(m,n)] +=
					     phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PI(Dx,DX,m,n)
					     +0.5*phi1*(PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
							PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
							PI(Dx,DX,m,n)*PSI(DX,i,j,k,l));
				       }
			// DV^3 block
			for (i=1;i<INTERNAL_ATOMS;i++)
			   for (j=0;j<DIM3;j++)
			      for (k=1;k<INTERNAL_ATOMS;k++)
				 for (l=0;l<DIM3;l++)
				    for (m=1;m<INTERNAL_ATOMS;m++)
				       for (n=0;n<DIM3;n++)
				       {
					  Phi[INDVV(i,j,k,l)][INDV(m,n)] +=
					     phi*OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,k,l)*OMEGA(Dx,p,q,m,n)
					     +phi1*(OMEGA(Dx,p,q,k,l)*SIGMA(p,q,i,j,m,n)
						    + OMEGA(Dx,p,q,i,j)*SIGMA(p,q,k,l,m,n)
						    + OMEGA(Dx,p,q,m,n)*SIGMA(p,q,i,j,k,l));
				       }
			// DU^2DV blocks
			for (i=0;i<DIM3;i++)
			   for (j=0;j<DIM3;j++)
			      for (k=0;k<DIM3;k++)
				 for (l=0;l<DIM3;l++)
				    for (m=1;m<INTERNAL_ATOMS;m++)
				       for (n=0;n<DIM3;n++)
				       {
					  Phi[INDUU(i,j,k,l)][INDV(m,n)] =
					     Phi[INDUV(i,j,m,n)][INDU(k,l)] =
					     Phi[INDVU(m,n,i,j)][INDU(k,l)] += (
						phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*OMEGA(Dx,p,q,m,n)
						+phi1*(PI(Dx,DX,k,l)*GAMMA(Dx,DX,p,q,i,j,m,n)
						       + PI(Dx,DX,i,j)*GAMMA(Dx,DX,p,q,k,l,m,n)
						       + 0.5*OMEGA(Dx,p,q,m,n)*PSI(DX,i,j,k,l))
						+phi2*THETA(DX,p,q,i,j,k,l,m,n));
				       }
			// DV^2DU blocks
			for (i=1;i<INTERNAL_ATOMS;i++)
			   for (j=0;j<DIM3;j++)
			      for (k=1;k<INTERNAL_ATOMS;k++)
				 for (l=0;l<DIM3;l++)
				    for (m=0;m<DIM3;m++)
				       for (n=0;n<DIM3;n++)
				       {
					  Phi[INDVV(i,j,k,l)][INDU(m,n)] =
					     Phi[INDVU(i,j,m,n)][INDV(k,l)] =
					     Phi[INDUV(m,n,i,j)][INDV(k,l)] += (
						phi*OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,k,l)*PI(Dx,DX,m,n)
						+phi1*(OMEGA(Dx,p,q,k,l)*GAMMA(Dx,DX,p,q,m,n,i,j)
						       + OMEGA(Dx,p,q,i,j)*GAMMA(Dx,DX,p,q,m,n,k,l)
						       + PI(Dx,DX,m,n)*SIGMA(p,q,i,j,k,l))
						+phi2*XI(p,q,i,j,k,l,m,n));
				       }
			break;
		     case PairPotentials::D4Y:
			 phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D3Y,dt);
			 phi2=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y,dt);
			 phi3=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,dt);

			 // DU^4 block 
			 for (i=0;i<DIM3;i++)
			    for (j=0;j<DIM3;j++)
			       for (k=0;k<DIM3;k++)
				  for (l=0;l<DIM3;l++)
				     for (m=0;m<DIM3;m++)
					for (n=0;n<DIM3;n++)
					   for (s=0;s<DIM3;s++)
					      for (t=0;t<DIM3;t++)
					      {
						 Phi[INDUU(i,j,k,l)][INDUU(m,n,s,t)]+=
						    phi*(PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*
							 PI(Dx,DX,m,n)*PI(Dx,DX,s,t)) +
						    0.5*phi1*(
						       PI(Dx,DX,k,l)*PI(Dx,DX,m,n)*PSI(DX,i,j,s,t)
						       + PI(Dx,DX,i,j)*PI(Dx,DX,m,n)*PSI(DX,k,l,s,t)
						       + PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PSI(DX,m,n,s,t)
						       + PI(Dx,DX,k,l)*PI(Dx,DX,s,t)*PSI(DX,i,j,m,n)
						       + PI(Dx,DX,i,j)*PI(Dx,DX,s,t)*PSI(DX,k,l,m,n)
						       + PI(Dx,DX,m,n)*PI(Dx,DX,s,t)*PSI(DX,i,j,k,l))
						    +0.25*phi2*(
						       PSI(DX,i,j,m,n)*PSI(DX,k,l,s,t)
						       + PSI(DX,i,j,s,t)*PSI(DX,k,l,m,n)
						       + PSI(DX,i,j,k,l)*PSI(DX,m,n,s,t));
					      }
			 // DV^4 block
			 for (i=1;i<INTERNAL_ATOMS;i++)
			    for (j=0;j<DIM3;j++)
			       for (k=1;k<INTERNAL_ATOMS;k++)
				  for (l=0;l<DIM3;l++)
				     for (m=1;m<INTERNAL_ATOMS;m++)
					for (n=0;n<DIM3;n++)
					   for (s=1;s<INTERNAL_ATOMS;s++)
					      for (t=0;t<DIM3;t++)
					      {
						 Phi[INDVV(i,j,k,l)][INDVV(m,n,s,t)] +=
						    phi*(OMEGA(Dx,p,q,i,j)
							 *OMEGA(Dx,p,q,k,l)
							 *OMEGA(Dx,p,q,m,n)
							 *OMEGA(Dx,p,q,s,t))
						    +phi1*(
						       OMEGA(Dx,p,q,k,l)*OMEGA(Dx,p,q,m,n)*SIGMA(p,q,i,j,s,t)
						       + OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,m,n)*SIGMA(p,q,k,l,s,t)
						       + OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,k,l)*SIGMA(p,q,m,n,s,t)
						       + OMEGA(Dx,p,q,k,l)*OMEGA(Dx,p,q,s,t)*SIGMA(p,q,i,j,m,n)
						       + OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,s,t)*SIGMA(p,q,k,l,m,n)
						       + OMEGA(Dx,p,q,m,n)*OMEGA(Dx,p,q,s,t)*SIGMA(p,q,i,j,k,l))
						    +phi2*(
						       SIGMA(p,q,i,j,m,n)*SIGMA(p,q,k,l,s,t)
						       + SIGMA(p,q,i,j,s,t)*SIGMA(p,q,k,l,m,n)
						       + SIGMA(p,q,i,j,k,l)*SIGMA(p,q,m,n,s,t));
					      }
			 // DU^3DV blocks
			 for (i=0;i<DIM3;i++)
			    for (j=0;j<DIM3;j++)
			       for (k=0;k<DIM3;k++)
				  for (l=0;l<DIM3;l++)
				     for (m=0;m<DIM3;m++)
					for (n=0;n<DIM3;n++)
					   for (s=1;s<INTERNAL_ATOMS;s++)
					      for (t=0;t<DIM3;t++)
					      {
						 Phi[INDUU(i,j,k,l)][INDUV(m,n,s,t)] =
						    Phi[INDUU(i,j,k,l)][INDVU(s,t,m,n)] =
						    Phi[INDUV(i,j,s,t)][INDUU(k,l,m,n)] =
						    Phi[INDVU(s,t,i,j)][INDUU(k,l,m,n)] += (
						       phi*(
							  PI(Dx,DX,i,j)
							  *PI(Dx,DX,k,l)
							  *PI(Dx,DX,m,n)
							  *OMEGA(Dx,p,q,s,t))
						       +phi1*(
							  PI(Dx,DX,k,l)*PI(Dx,DX,m,n)*GAMMA(Dx,DX,p,q,i,j,s,t)
							  + PI(Dx,DX,i,j)*PI(Dx,DX,m,n)*GAMMA(Dx,DX,p,q,k,l,s,t)
							  + PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*GAMMA(Dx,DX,p,q,m,n,s,t)
							  +0.5*(
							     PI(Dx,DX,k,l)*OMEGA(Dx,p,q,s,t)*PSI(DX,i,j,m,n)
							     + PI(Dx,DX,i,j)*OMEGA(Dx,p,q,s,t)*PSI(DX,k,l,m,n)
							     + PI(Dx,DX,m,n)*OMEGA(Dx,p,q,s,t)*PSI(DX,i,j,k,l)))
						       +phi2*(
							  PI(Dx,DX,k,l)*THETA(DX,p,q,i,j,m,n,s,t)
							  + PI(Dx,DX,i,j)*THETA(DX,p,q,k,l,m,n,s,t)
							  + PI(Dx,DX,m,n)*THETA(DX,p,q,i,j,k,l,s,t)
							  +0.5*(
							     GAMMA(Dx,DX,p,q,k,l,s,t)*PSI(DX,i,j,m,n)
							     + GAMMA(Dx,DX,p,q,i,j,s,t)*PSI(DX,k,l,m,n)
							     + GAMMA(Dx,DX,p,q,m,n,s,t)*PSI(DX,i,j,k,l))));
					      }
			 // DV^3DU blocks
			 for (i=1;i<INTERNAL_ATOMS;i++)
			    for (j=0;j<DIM3;j++)
			       for (k=1;k<INTERNAL_ATOMS;k++)
				  for (l=0;l<DIM3;l++)
				     for (m=1;m<INTERNAL_ATOMS;m++)
					for (n=0;n<DIM3;n++)
					   for (s=0;s<DIM3;s++)
					      for (t=0;t<DIM3;t++)
					      {
						 Phi[INDVV(i,j,k,l)][INDVU(m,n,s,t)] =
						    Phi[INDVV(i,j,k,l)][INDUV(s,t,m,n)] =
						    Phi[INDVU(i,j,s,t)][INDVV(k,l,m,n)] =
						    Phi[INDUV(s,t,i,j)][INDVV(k,l,m,n)] += (
						       phi*(
							  OMEGA(Dx,p,q,i,j)
							  *OMEGA(Dx,p,q,k,l)
							  *OMEGA(Dx,p,q,m,n)
							  *PI(Dx,DX,s,t))
						       +phi1*(
							  OMEGA(Dx,p,q,k,l)*OMEGA(Dx,p,q,m,n)*GAMMA(Dx,DX,p,q,s,t,i,j)
							  + OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,m,n)*GAMMA(Dx,DX,p,q,s,t,k,l)
							  + OMEGA(Dx,p,q,i,j)*OMEGA(Dx,p,q,k,l)*GAMMA(Dx,DX,p,q,s,t,m,n)
							  + OMEGA(Dx,p,q,k,l)*PI(Dx,DX,s,t)*SIGMA(p,q,i,j,m,n)
							  + OMEGA(Dx,p,q,i,j)*PI(Dx,DX,s,t)*SIGMA(p,q,k,l,m,n)
							  + OMEGA(Dx,p,q,m,n)*PI(Dx,DX,s,t)*SIGMA(p,q,i,j,k,l))
						       +phi2*(
							  OMEGA(Dx,p,q,k,l)*XI(p,q,i,j,m,n,s,t)
							  + OMEGA(Dx,p,q,i,j)*XI(p,q,k,l,m,n,s,t)
							  + OMEGA(Dx,p,q,m,n)*XI(p,q,i,j,k,l,s,t)
							  + SIGMA(p,q,i,j,m,n)*GAMMA(Dx,DX,p,q,s,t,k,l)
							  + SIGMA(p,q,k,l,m,n)*GAMMA(Dx,DX,p,q,s,t,i,j)
							  + SIGMA(p,q,i,j,k,l)*GAMMA(Dx,DX,p,q,s,t,m,n)));
					      }
			 // DU^2DV^2 blocks
			 for (i=0;i<DIM3;i++)
			    for (j=0;j<DIM3;j++)
			       for (k=0;k<DIM3;k++)
				  for (l=0;l<DIM3;l++)
				     for (m=1;m<INTERNAL_ATOMS;m++)
					for (n=0;n<DIM3;n++)
					   for (s=1;s<INTERNAL_ATOMS;s++)
					      for (t=0;t<DIM3;t++)
					      {
						 Phi[INDUU(i,j,k,l)][INDVV(m,n,s,t)] =
						    Phi[INDUV(i,j,m,n)][INDUV(k,l,s,t)] =
						    Phi[INDUV(i,j,m,n)][INDVU(s,t,k,l)] =
						    Phi[INDVU(m,n,i,j)][INDUV(k,l,s,t)] =
						    Phi[INDVU(m,n,i,j)][INDVU(s,t,k,l)] =
						    Phi[INDVV(m,n,s,t)][INDUU(i,j,k,l)] += (
						       phi*(
							  PI(Dx,DX,i,j)
							  *PI(Dx,DX,k,l)
							  *OMEGA(Dx,p,q,m,n)
							  *OMEGA(Dx,p,q,s,t))
						       +phi1*(
							  PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*SIGMA(p,q,m,n,s,t)
							  + PI(Dx,DX,k,l)*OMEGA(Dx,p,q,m,n)*GAMMA(Dx,DX,p,q,i,j,s,t)
							  + PI(Dx,DX,i,j)*OMEGA(Dx,p,q,m,n)*GAMMA(Dx,DX,p,q,k,l,s,t)
							  + PI(Dx,DX,k,l)*OMEGA(Dx,p,q,s,t)*GAMMA(Dx,DX,p,q,i,j,m,n)
							  + PI(Dx,DX,i,j)*OMEGA(Dx,p,q,s,t)*GAMMA(Dx,DX,p,q,k,l,m,n)
							  + OMEGA(Dx,p,q,m,n)*OMEGA(Dx,p,q,s,t)*0.5*PSI(DX,i,j,k,l))
						       +phi2*(
							  PI(Dx,DX,k,l)*XI(p,q,m,n,s,t,i,j)
							  + PI(Dx,DX,i,j)*XI(p,q,m,n,s,t,k,l)
							  + OMEGA(Dx,p,q,m,n)*THETA(DX,p,q,i,j,k,l,s,t)
							  + OMEGA(Dx,p,q,s,t)*THETA(DX,p,q,i,j,k,l,m,n)
							  + GAMMA(Dx,DX,p,q,i,j,m,n)*GAMMA(Dx,DX,p,q,k,l,s,t)
							  + GAMMA(Dx,DX,p,q,i,j,s,t)*GAMMA(Dx,DX,p,q,k,l,m,n)
							  + SIGMA(p,q,m,n,s,t)*0.5*PSI(DX,i,j,k,l))
						       +phi3*LAMDA(p,q,i,j,k,l,m,n,s,t));
					      }
			 break;
		  }
	       }
	    }
	 }
      }
   }

   // BodyForce[i] = BodyForce[i] / ForceNorm
   for (i=0;i<INTERNAL_ATOMS;i++)
   {
      for (j=0;j<DIM3;j++)
      {
	 BodyForce_[i][j] /= ForceNorm;
      }
   }
   
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLen_*RefLen_*RefLen_*LatticeVec_.Det()*ShearMod_));

   if (moduliflag)
   {
      return Phi;
   }

   // Add in External Work Terms
   // Recall that Pressure_ is applied stress
   // Thus E = W - pJ

   if (dt == PairPotentials::T0)
   {
      switch (dy)
      {
	 case PairPotentials::Y0:
	    Phi[0][0] = Phi[0][0] - Pressure_*J;
	    break;
	 case PairPotentials::DY:
	 {
	    Matrix Uinv = U.Inverse();
	    for (i=0;i<DIM3;i++)
	       for (j=0;j<DIM3;j++)
	       {
		  Phi[0][INDU(i,j)] -= Pressure_*J*Uinv[i][j];
	       }
	 }
	 break;
	 case PairPotentials::D2Y:
	    for (i=0;i<DIM3;i++)
	       for (j=0;j<DIM3;j++)
		  for (k=0;k<DIM3;k++)
		     for (l=0;l<DIM3;l++)
		     {
			for (q=0;q<DIM3;q++)
			   for (s=0;s<DIM3;s++)
			   {
			      Phi[INDU(i,j)][INDU(k,l)] -=
				 (Pressure_/8.0)*(
				    Alt[i][k][q]*Alt[j][l][s] + Alt[j][k][q]*Alt[i][l][s] +
				    Alt[i][l][q]*Alt[j][k][s] + Alt[j][l][q]*Alt[i][k][s] +
				    Alt[i][q][k]*Alt[j][s][l] + Alt[j][q][k]*Alt[i][s][l] +
				    Alt[i][q][l]*Alt[j][s][k] + Alt[j][q][l]*Alt[i][s][k]
				    )*U[q][s];
			   }
		     }
	    break;
	 case PairPotentials::D3Y:
	    for (i=0;i<DIM3;i++)
	       for (j=0;j<DIM3;j++)
		  for (k=0;k<DIM3;k++)
		     for (l=0;l<DIM3;l++)
			for (q=0;q<DIM3;q++)
			   for (s=0;s<DIM3;s++)
			   {
			      Phi[INDUU(i,j,k,l)][INDU(q,s)] -=
				 (Pressure_/16.0)*(
				    Alt[i][k][q]*Alt[j][l][s] + Alt[j][k][q]*Alt[i][l][s] +
				    Alt[i][l][q]*Alt[j][k][s] + Alt[j][l][q]*Alt[i][k][s] +
				    Alt[i][q][k]*Alt[j][s][l] + Alt[j][q][k]*Alt[i][s][l] +
				    Alt[i][q][l]*Alt[j][s][k] + Alt[j][q][l]*Alt[i][s][k] +
				    Alt[i][k][s]*Alt[j][l][q] + Alt[j][k][s]*Alt[i][l][q] +
				    Alt[i][l][s]*Alt[j][k][q] + Alt[j][l][s]*Alt[i][k][q] +
				    Alt[i][s][k]*Alt[j][q][l] + Alt[j][s][k]*Alt[i][q][l] +
				    Alt[i][s][l]*Alt[j][q][k] + Alt[j][s][l]*Alt[i][q][k]);
			   }
	    break;
	 case PairPotentials::D4Y:
	    // Terms are zero
	    break;
      }
   }

   return Phi;
}

Matrix NiTi9TPPLat::CondensedModuli()
{
   Matrix stiffness = Phi(1,PairPotentials::D2Y);
   int intrn = DOFS-6;
   Matrix CM(6,6), IM(intrn,intrn);
   
   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
      {
	 CM[i][j] = stiffness[i][j];
      }
   
   for (int i=0;i<intrn;i++)
      for (int j=0;j<intrn;j++)
      {
	 IM[i][j] = stiffness[6+i][6+j];
      }
   IM = IM.Inverse();

   // Set up Condensed Moduli
   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
      {
	 for (int m=0;m<intrn;m++)
	    for (int n=0;n<intrn;n++)
	    {
	       CM[i][j] -= stiffness[i][6+m]*IM[m][n]*stiffness[6+n][j];
	    }
      }

   // Remove 2's and 4's
   for (int i=3;i<6;i++)
   {
      for (int j=0;j<3;j++)
      {
	 CM[i][j] /= 2.0;
	 CM[j][i] /= 2.0;
      }

      for (int j=3;j<6;j++)
      {
	 CM[i][j] /= 4.0;
      }
   }

   return CM;
}

CMatrix NiTi9TPPLat::DynamicalStiffness(Vector &Y)
{
   static CMatrix Cy;
   static Matrix U(DIM3,DIM3);
   static Matrix V(INTERNAL_ATOMS,3);
   static Matrix Eigvals(1,DIM3);
   static double X[DIM3];
   static Vector DX(DIM3),Dx(DIM3);
   static Vector Direction(DIM3);
   static double pi = 4.0*atan(1.0);
   static double J;
   static int p,q;
   static int i,j,k,l,m,n,s,t;
   static double r2,phi,phi1,Influancedist[DIM3],tmp;
   static int Top[DIM3],Bottom[DIM3],CurrentInfluanceDist;
   static interaction Inter;
   static complex<double> Ic(0,1);

   Cy.Resize(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3,0.0);
   
   U[0][0] = DOF_[0];
   U[1][1] = DOF_[1];
   U[2][2] = DOF_[2];
   U[0][1] = U[1][0] = DOF_[3];
   U[0][2] = U[2][0] = DOF_[4];
   U[1][2] = U[2][1] = DOF_[5];

   V[0][0] = 0.0;
   V[0][1] = 0.0;
   V[0][2] = 0.0;
   V[1][0] = DOF_[6];
   V[1][1] = DOF_[7];
   V[1][2] = DOF_[8];

   // find largest eigenvalue of the inverse transformation
   // (i.e. from current to ref) and use influence cube of
   // that size...
   //
   // Use the fact that eigs of Uinv = 1/ eigs of U.
   J = U.Det();
   Eigvals = SymEigVal(U);
   tmp = Eigvals[0][0];
   for (i=0;i<DIM3;i++)
      if (Eigvals[0][i] < tmp) tmp = Eigvals[0][i];
   
   // Set to inverse eigenvalue
   tmp = 1.0/tmp;
   for (i=0;i<DIM3;i++)
   {
      Influancedist[i]=tmp*InfluanceDist_;
   }
   
   for (p=0;p<DIM3;p++)
   {
      // set influance distance based on cube size
      //
      // also setup to be large enough to encompass Eulerian sphere
      CurrentInfluanceDist = int(ceil(Influancedist[p]));

      Top[p] = CurrentInfluanceDist;
      Bottom[p] = -CurrentInfluanceDist;
   }

   for (p=0;p<INTERNAL_ATOMS;p++)
   {
      for (q=0;q<INTERNAL_ATOMS;q++)
      {
	 Inter = INTER[p][q];
	 
	 for (X[0] = Bottom[0];X[0] <= Top[0];X[0]++)
	 {
	    for (X[1] = Bottom[1];X[1] <= Top[1];X[1]++)
	    {
	       for (X[2] = Bottom[2];X[2] <= Top[2];X[2]++)
	       {
		  for (i=0;i<DIM3;i++)
		  {
		     DX[i] = 0.0;
				     
		     for (j=0;j<DIM3;j++)
		     {
			DX[i] += (X[j] + A[q][j] - A[p][j] + V[q][j] - V[p][j])*
			   LatticeVec_[j][i]*RefLen_;

		     }
		  }

		  r2 = 0.0;
		  for (i=0;i<DIM3;i++)
		  {
		     Dx[i] = 0.0;

		     for (j=0;j<DIM3;j++)
		     {
			Dx[i] += U[i][j] * DX[j];
		     }
		     r2 += Dx[i]*Dx[i];
		  }
		  // Only use Sphere of Influance (current)
		  if (r2==0 || r2 > InfluanceDist_*InfluanceDist_)
		  {
		     continue;
		  }

		  // Calculate Cy
		  if (p != q)
		  {
		     for (i=0;i<DIM3;i++)
			for (j=0;j<DIM3;j++)
			{
			   Cy[DIM3*p+i][DIM3*p+j] +=
			      (2.0*Del(i,j)*Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY)
			       +4.0*Dx[i]*Dx[j]*Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y))
			      *exp(-2.0*pi*Ic * (Y*Dx));

			   Cy[DIM3*p+i][DIM3*q+j] +=
			      (-2.0*Del(i,j)*Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY)
			       -4.0*Dx[i]*Dx[j]*Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y))
			      *exp(-2.0*pi*Ic * (Y*Dx));
			}
		     }
	       }
	    }
	 }

	 // Normalize through the Mass Matrix
	 for (i=0;i<DIM3;i++)
	    for (j=0;j<DIM3;j++)
	    {
	       Cy[DIM3*p+i][DIM3*q+j] /= sqrt(AtomicMass_[p]*AtomicMass_[q]);
	    }
      }
   }
   return Cy;
}

#include <fstream.h>

int NiTi9TPPLat::DynamicallyStable(Vector &Y)
{
   int Depth = 4;
   double offset;
   double onehalf = 0.5;

   fstream out;

   out.open("test.out",ios::out);

   // I have changed the k values so as to only
   // cover 1/2 of the cube (positive Z)

   
   // Place nodes at corners of cube
   for (int k=1;k<2;k++)
      for (int j=0;j<2;j++)
	 for (int i=0;i<2;i++)
	 {
	    out << setw(20) << i-1/2.0
		<< setw(20) << j-1/2.0
		<< setw(20) << k-1/2.0
		<< endl;
	 }

   for (int n=0;n<Depth;n++)
   {
      double twon = pow(2,n),
	 twon1 = pow(2,n+1);
      int twonm1 = int((n!=0) ? pow(2,n-1) : 0);
      
      offset = (twon-1)/twon1;

      // Place nodes at centroid of n level squares.
      for (int k=twonm1;k<twon;k++)
      {
	 for (int j=0;j<twon;j++)
	 {
	    for (int i=0;i<twon;i++)
	    {
	       out << setw(20) << (i/twon)-offset
		   << setw(20) << (j/twon)-offset
		   << setw(20) << (k/twon)-offset
		   << endl;
	    }
	 }
      }

      // Place nodes at corners of n+1 level squares.
      for (int k=twonm1;k<twon;k++)
      {
	 for (int j=0;j<=twon;j++)
	 {
	    for (int i=0;i<=twon1;i++)
	    {
	       out << setw(20) << (i/twon1) - onehalf
		   << setw(20) << (j/twon) - onehalf
		   << setw(20) << (k/twon) - offset
		   << endl;
	    }
	 }
	 
	 for (int j=0;j<twon;j++)
	 {
	    for (int i=0;i<=twon;i++)
	    {
	       out << setw(20) << (i/twon) - onehalf
		   << setw(20) << (j/twon) - offset
		   << setw(20) << (k/twon) - offset
		   << endl;
	    }
	 }
      }

      out << endl << endl;

      if (twonm1 == 0) twonm1++;
      for (int k=twonm1;k<=twon;k++)
      {
	 for (int j=0;j<twon;j++)
	 {
	    for (int i=0;i<=twon1;i++)
	    {
	       out << setw(20) << (i/twon1) - onehalf
		   << setw(20) << (j/twon) - offset
		   << setw(20) << (k/twon) - onehalf
		   << endl;
	    }
	 }

	 for (int j=0;j<=twon;j++)
	 {
	    for (int i=0;i<twon;i++)
	    {
	       out << setw(20) << (i/twon) - offset
		   << setw(20) << (j/twon) - onehalf
		   << setw(20) << (k/twon) - onehalf
		   << endl;
	    }
	 }
      }

      out << endl << endl;
   }
   
   
   out.close();
   return 1;
}

void NiTi9TPPLat::Print(ostream &out,PrintDetail flag)
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
      EigenValues(1,DOFS),
      CondEV(1,6);

   EigenValues=SymEigVal(stiffness);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<DOFS;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   Matrix
      CondModuli = CondensedModuli();

   CondEV=SymEigVal(CondModuli);
   int RankOneConvex = FullScanRank1Convex3D(CondModuli,ConvexityDX_);

   switch (flag)
   {
      case PrintLong:
	 out << "NiTi9TPPLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	     << "Atomic Mass 0        : " << setw(W) << AtomicMass_[0] << endl
	     << "Atomic Mass 1        : " << setw(W) << AtomicMass_[1] << endl
	     << "Potential Parameters : " << endl
	     << "AA -- " << setw(W) << Potential_[aa] << endl
	     << "BB -- " << setw(W) << Potential_[bb] << endl
	     << "AB -- " << setw(W) << Potential_[ab] << endl
	     << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 cout << "NiTi9TPPLat:" << endl << endl
	      << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	      << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	      << "Atomic Mass 0        : " << setw(W) << AtomicMass_[0] << endl
	      << "Atomic Mass 1        : " << setw(W) << AtomicMass_[1] << endl
	      << "Potential Parameters : " << endl
	      << "AA -- " << setw(W) << Potential_[aa] << endl
	      << "BB -- " << setw(W) << Potential_[bb] << endl
	      << "AB -- " << setw(W) << Potential_[ab] << endl
	      << "Shear Modulus : " << setw(W) << ShearMod_ << endl;	 
	 // passthrough to short
      case PrintShort:
	 out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	     << "Pressure (G Normalized): " << setw(W) << Pressure_ << endl
	     << "DOF's :" << endl << setw(W) << DOF_ << endl
	     << "Potential Value (G Normalized):" << setw(W) << energy << endl
	     << "BodyForce Value 0 (Inf Normalized):" << setw(W) << BodyForce_[0] << endl
	     << "BodyForce Value 1 (Inf Normalized):" << setw(W) << BodyForce_[1] << endl
	     << "Stress (G Normalized):" << setw(W) << stress << endl
	     << "Stiffness (G Normalized):" << setw(W) << stiffness
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl
	     << "Condensed Moduli (G Normalized):" << setw(W) << CondModuli
	     << "CondEV Info:" << setw(W) << CondEV
	     << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl;
	 cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	      << "Pressure (G Normalized): " << setw(W) << Pressure_ << endl
	      << "DOF's :" << endl << setw(W) << DOF_ << endl
	      << "Potential Value (G Normalized):" << setw(W) << energy << endl
	      << "BodyForce Value 0 (Inf Normalized):" << setw(W) << BodyForce_[0] << endl
	      << "BodyForce Value 1 (Inf Normalized):" << setw(W) << BodyForce_[1] << endl
	      << "Stress (G Normalized):" << setw(W) << stress << endl
	      << "Stiffness (G Normalized):" << setw(W) << stiffness
	      << "Eigenvalue Info:"  << setw(W) << EigenValues
	      << "Bifurcation Info:" << setw(W) << MinEigVal
	      << setw(W) << NoNegEigVal << endl
	      << "Condensed Moduli (G Normalized):" << setw(W) << CondModuli
	      << "CondEV Info:" << setw(W) << CondEV
	      << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl;
	 Vector Z(DIM3,0.0);
	 out << "Dynamical Stability" << endl << setw(W) << DynamicalStiffness(Z) << endl;
	 cout << "Dynamical Stability" << endl << setw(W) << DynamicalStiffness(Z) << endl;

	 DynamicallyStable(Z);
	 break;
   }
}

ostream &operator<<(ostream &out,NiTi9TPPLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double NiTi9TPPLat::Alt[DIM3][DIM3][DIM3] = {0.0, 0.0, 0.0,
						    0.0, 0.0, 1.0,
						    0.0, -1.0, 0.0,
						    0.0, 0.0, -1.0,
						    0.0, 0.0, 0.0,
						    1.0, 0.0, 0.0,
						    0.0, 1.0, 0.0,
						    -1.0, 0.0, 0.0,
						    0.0, 0.0, 0.0};

const double NiTi9TPPLat::A[INTERNAL_ATOMS][DIM3] = {0.0, 0.0, 0.0,
						     0.5, 0.5, 0.5};

const NiTi9TPPLat::interaction NiTi9TPPLat::INTER[INTERNAL_ATOMS][INTERNAL_ATOMS] = {
   aa, ab,
   ab, bb};
