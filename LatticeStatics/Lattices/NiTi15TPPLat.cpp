#include "NiTi15TPPLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTi15TPPLat::NiTi15TPPLat(char *datafile)
{
   // First Size DOF
   DOF_.Resize(DOFS);
   // Set LatticeVec_
   LatticeVec_.Resize(DIM3,DIM3,0.0);
   LatticeVec_[0][0] = 1.0;
   LatticeVec_[0][1] = -1.0;
   
   LatticeVec_[1][0] = 1.0;
   LatticeVec_[1][1] = 1.0;

   LatticeVec_[2][2] = 1.0;
   
   // Setup Bodyforce_
   BodyForce_[0].Resize(DIM3,0.0);
   BodyForce_[1].Resize(DIM3,0.0);
   BodyForce_[2].Resize(DIM3,0.0);
   BodyForce_[3].Resize(DIM3,0.0);

   // Get Potential Parameters
   double Tref,A0,B0,Alpha,Rref,Rtheta,Tmelt;
   GetParameter("^Tref",datafile,"%lf",&Tref);
   GetParameter("^A0_aa",datafile,"%lf",&A0);
   GetParameter("^B0_aa",datafile,"%lf",&B0);
   GetParameter("^Alpha_aa",datafile,"%lf",&Alpha);
   GetParameter("^Rref_aa",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_aa",datafile,"%lf",&Rtheta);
   GetParameter("Tmelt_aa",datafile,"%lf",&Tmelt);
   Potential_[aa]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);

   GetParameter("^A0_bb",datafile,"%lf",&A0);
   GetParameter("^B0_bb",datafile,"%lf",&B0);
   GetParameter("^Alpha_bb",datafile,"%lf",&Alpha);
   GetParameter("^Rref_bb",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_bb",datafile,"%lf",&Rtheta);
   GetParameter("Tmelt_bb",datafile,"%lf",&Tmelt);
   Potential_[bb]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);

   GetParameter("^A0_ab",datafile,"%lf",&A0);
   GetParameter("^B0_ab",datafile,"%lf",&B0);
   GetParameter("^Alpha_ab",datafile,"%lf",&Alpha);
   GetParameter("^Rref_ab",datafile,"%lf",&Rref);
   GetParameter("^Rtheta_ab",datafile,"%lf",&Rtheta);
   GetParameter("Tmelt_ab",datafile,"%lf",&Tmelt);
   Potential_[ab]=RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);
   
   // Get Lattice parameters
   GetParameter("^RefLen",datafile,"%lf",&RefLen_);
   GetParameter("^InfluanceDist",datafile,"%u",&InfluanceDist_);
   GetParameter("^NTemp",datafile,"%lf",&NTemp_);
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

int NiTi15TPPLat::FindLatticeSpacing(int iter,double dx)
{
   double oldPressure=Pressure_,
      oldTemp=NTemp_;

   Pressure_=0.0;
   NTemp_=1.0;
   ShearMod_=1.0;
   DOF_[0] = DOF_[1] = DOF_[2] = 1.0;
   for (int i=3;i<15;i++)
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

   ShearMod_ = 0.25*fabs(Stiffness()[5][5]);
   NTemp_=oldTemp;
   Pressure_=oldPressure;

   return 0;
}   

// Lattice Routines

inline double NiTi15TPPLat::PI(const Vector &Dx,const Vector &DX,
			       int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

inline double NiTi15TPPLat::PSI(const Vector &DX,
				int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}

inline double NiTi15TPPLat::OMEGA(const Vector &Dx,int p,int q,int i, int j)
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

inline double NiTi15TPPLat::SIGMA(int p,int q,int i,int j,int k,int l)
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

inline double NiTi15TPPLat::GAMMA(const Vector &Dx,const Vector &DX,
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

inline double NiTi15TPPLat::THETA(const Vector &DX,int p,int q,int i,int j,int k,int l,
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

inline double NiTi15TPPLat::XI(int p,int q,int i,int j,int k,int l,int m,int n)
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

inline double NiTi15TPPLat::LAMDA(int p,int q,int i,int j,int k,int l,int m,int n,int a,
				  int b)
{
   return (0.5*DELTA(m,p,q)*DELTA(a,p,q)*(Del(i,k)*LatticeVec_[n][j]*LatticeVec_[b][l]
					  + Del(i,k)*LatticeVec_[b][j]*LatticeVec_[n][l]
					  + Del(i,l)*LatticeVec_[n][j]*LatticeVec_[b][k]
					  + Del(i,l)*LatticeVec_[b][j]*LatticeVec_[n][k]
					  + Del(j,k)*LatticeVec_[n][i]*LatticeVec_[b][l]
					  + Del(j,k)*LatticeVec_[b][i]*LatticeVec_[n][l]
					  + Del(j,l)*LatticeVec_[n][j]*LatticeVec_[b][j]
					  + Del(j,l)*LatticeVec_[b][i]*LatticeVec_[n][k]));
}

double NiTi15TPPLat::pwr(const double &x,const unsigned y)
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

inline int NiTi15TPPLat::INDU(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

inline int NiTi15TPPLat::INDUU(int k,int l,int m,int n)
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

inline int NiTi15TPPLat::INDV(int i,int j)
{
   if (!i)
   {
      cerr << "Error: INDV(i,j) i==0!!!!!" << endl;
      exit(-1);
   }

   return 6 + (i-1)*3 + j;
}

inline int NiTi15TPPLat::INDVV(int k,int l,int m,int n)
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

inline int NiTi15TPPLat::INDUV(int i,int j,int m,int n)
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

inline int NiTi15TPPLat::INDVU(int m,int n,int i,int j)
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

Matrix NiTi15TPPLat::Phi(unsigned moduliflag,PairPotentials::YDeriv dy,
			 PairPotentials::TDeriv dt)
{
   static Matrix U(3,3);
   static Matrix V(4,3);
   static Matrix Eigvals(1,3);
   static double X[3];
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
	 Phi.Resize(1,15,0.0);
	 break;
      case PairPotentials::D2Y:
	 Phi.Resize(15,15,0.0);
	 break;
      case PairPotentials::D3Y:
	 Phi.Resize(225,15,0.0);
	 break;
      case PairPotentials::D4Y:
	 Phi.Resize(225,225,0.0);
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
   V[2][0] = DOF_[9];
   V[2][1] = DOF_[10];
   V[2][2] = DOF_[11];
   V[3][0] = DOF_[12];
   V[3][1] = DOF_[13];
   V[3][2] = DOF_[14];

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
			//Lower Diag Block (9,9)
			for (i=1;i<4;i++)
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


double NiTi15TPPLat::Energy()
{
   return Phi()[0][0];
}

Matrix NiTi15TPPLat::Stress()
{
   return Phi(0,PairPotentials::DY);
}

Matrix NiTi15TPPLat::StressDT()
{
   return Phi(0,PairPotentials::DY,PairPotentials::DT);
}

Matrix NiTi15TPPLat::Stiffness()
{
   return Phi(0,PairPotentials::D2Y);
}

Matrix NiTi15TPPLat::Moduli()
{
   return Phi(1,PairPotentials::D2Y);
}

int NiTi15TPPLat::StiffnessNulity(double *Min)
{
   int NoNegEigVal = 0;
   int index = 0;

   Matrix EigenValues(1,DOFS);

   EigenValues=SymEigVal(Stiffness());
   if (Min != NULL) *Min = fabs(EigenValues[0][0]);
   for (int i=0;i<DOFS;i++)
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

void NiTi15TPPLat::CriticalPointInfo(double Tolerance,int Width,ostream &out)
{
   Matrix E3=Phi(0,PairPotentials::D3Y),
      E2=Stiffness(),
      E2T=Phi(0,PairPotentials::D2Y,PairPotentials::DT),
      EigVec,
      EigVal=SymEigVal(E2,&EigVec);
   
   // Matrix E4=Phi(0,PairPotential::D4Y);

   Matrix Mode;
   double Eijk[DOFS][DOFS][DOFS],
      EijT[DOFS][DOFS];

   // Find the modes
   int count = 0,
      Ind[DOFS];

   for (int i=0;i<DOFS;i++)
      if (fabs(EigVal[0][i]) < Tolerance)
      {
	 Ind[count++]=i;
      }
   
   Mode.Resize(count,DOFS);

   for (int i=0;i<count;i++)
   {
      for (int j=0;j<DOFS;j++)
      {
	 Mode[i][j] = EigVec[j][Ind[i]];
      }
   }

   // Eijk
   for (int i=0;i<count;i++)
      for (int j=0;j<count;j++)
	 for (int k=0;k<count;k++)
	 {
	    Eijk[i][j][k] = 0.0;
	    for (int a=0;a<DOFS;a++)
	       for (int b=0;b<DOFS;b++)
		  for (int c=0;c<DOFS;c++)
		  {
		     Eijk[i][j][k] += E3[a*DOFS + b][c]*Mode[i][a]*Mode[j][b]*Mode[k][c];
		  }
	 }

   //EijT
   for (int i=0;i<count;i++)
      for (int j=0;j<count;j++)
      {
	 EijT[i][j] = 0.0;
	 for (int a=0;a<DOFS;a++)
	    for (int b=0;b<DOFS;b++)
	    {
	       EijT[i][j] += E2T[a][b]*Mode[i][a]*Mode[j][b];
	    }
      }
   
   // Print out results
   for (int i=0;i<70;i++)
   {
      cout << "-";
      out << "-";
   }
   cout << endl << endl << "Bifurcation Equations:" << endl;
   out << endl << endl << "Bifurcation Equations:" << endl;

   for (int i=0;i<count;i++)
   {
      for (int j=0;j<count;j++)
	 for (int k=0;k<count;k++)
	 {
	    cout << "("
		 << Eijk[i][j][k]
		 << ")a_" << j
		 << "a_"  << k
		 << " + ";
	    out << "("
		<< Eijk[i][j][k]
		<< ")a_" << j
		<< "a_"  << k
		<< " + ";
	 }

      cout << "2T_1( ";
      out << "2T_1( ";
      
      for (int j=0;j<count-1;j++)
      {
	 cout << "("
	      << EijT[i][j]
	      << ")a_" << j
	      << " + ";
	 out << "("
	     << EijT[i][j]
	     << ")a_" << j
	     << " + ";
      }
      cout << "(" << EijT[i][count-1] << ")a_" << count-1
	   << ") = 0" << endl;
      out << "(" << EijT[i][count-1] << ")a_" << count-1
	  << ") = 0" << endl;
   }

   for (int i=0;i<70;i++)
   {
      cout << "-";
      out << "-";
   }
   cout << endl;
   out << endl;
   
   return;
}

void NiTi15TPPLat::Print(ostream &out,PrintDetail flag)
{
   int W=out.width();

   out.width(0);

   double MinEigVal;
   int NoNegEigVal=0;
   
   Matrix
      stiffness = Stiffness(),
      moduli = Moduli(),
      EigenValues(1,DOFS);

   EigenValues=SymEigVal(stiffness);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<DOFS;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   switch (flag)
   {
      case PrintLong:
	 out << "NiTi15TPPLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
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
	     << "Potential Value (G Normalized):" << setw(W) << Energy() << endl
	     << "BodyForce Value 0 (Inf Normalized):" << setw(W) << BodyForce_[0] << endl
	     << "BodyForce Value 1 (Inf Normalized):" << setw(W) << BodyForce_[1] << endl
	     << "BodyForce Value 2 (Inf Normalized):" << setw(W) << BodyForce_[2] << endl
	     << "BodyForce Value 3 (Inf Normalized):" << setw(W) << BodyForce_[3] << endl
	     << "Stress (G Normalized):" << setw(W) << Stress() << endl
	     << "Stiffness (G Normalized):" << setw(W) << stiffness
	     << "Rank 1 Convex:" << setw(W)
	     << Rank1Convex3D(moduli,ConvexityDX_) << endl
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl;
	 break;
   }
}

ostream &operator<<(ostream &out,NiTi15TPPLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double NiTi15TPPLat::Alt[DIM3][DIM3][DIM3] = {0.0, 0.0, 0.0,
						    0.0, 0.0, 1.0,
						    0.0, -1.0,0.0,
						    0.0, 0.0, -1.0,
						    0.0, 0.0, 0.0,
						    1.0, 0.0, 0.0,
						    0.0, 1.0, 0.0,
						    -1.0,0.0, 0.0,
						    0.0, 0.0, 0.0};

const double NiTi15TPPLat::A[INTERNAL_ATOMS][DIM3] = {0.0, 0.0, 0.0,
						      0.5, 0.5, 0.0,
						      0.5, 0.0, 0.5,
						      0.0, 0.5, 0.5};

const NiTi15TPPLat::interaction NiTi15TPPLat::INTER[INTERNAL_ATOMS][INTERNAL_ATOMS] = {
   aa, aa, ab, ab,
   aa, aa, ab, ab,
   ab, ab, bb, bb,
   ab, ab, bb, bb};
