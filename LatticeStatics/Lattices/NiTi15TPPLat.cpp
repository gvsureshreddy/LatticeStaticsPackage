#include "NiTi15TPPLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTi15TPPLat::NiTi15TPPLat(char *datafile)
{
   // First Size DOF
   DOF_.Resize(DOFS);
   // Set LatticeVec_
   LatticeVec_[0].Resize(DIM3,0.0);
   LatticeVec_[1].Resize(DIM3,0.0);
   LatticeVec_[2].Resize(DIM3,0.0);
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

   ShearMod_ = fabs(Stiffness()[5][5]);
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
   double ret;
   
   if (DELTA(i,p,q))
   {
      ret=0;
      for (int s=0;s<DIM3;s++)
      {
	 ret += DOF_[INDU(j,s)]*Dx[s] + Dx[s]*DOF_[INDU(s,j)];
      }
   }
   else
   {
      ret=0;
   }
   
   return ret;
}

inline double NiTi15TPPLat::SIGMA(int p,int q,int i,int j,int k,int l)
{
   double U2=0;
   for (int s=0;s<DIM3;s++)
      U2 += DOF_[INDU(j,s)]*DOF_[INDU(s,l)];
   
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*U2;
}

inline double NiTi15TPPLat::GAMMA(const Vector &Dx,const Vector &DX,
				  int p,int q,int i,int j,int k,int l)
{
   return DELTA(k,p,q)*(Del(j,l)*Dx[i] + Del(i,l)*Dx[j] +
			DOF_[INDU(l,i)]*DX[j] + DOF_[INDU(l,j)]*DX[i]);
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

inline int NiTi15TPPLat::INDU(int k,int l,int m,int n)
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

inline int NiTi15TPPLat::INDV(int i,int j)
{
   if (!i)
   {
      cerr << "Error: INDV(i,j) i==0!!!!!" << endl;
      exit(-1);
   }

   return 6 + (i-1)*3 + j;
}

inline int NiTi15TPPLat::INDV(int k,int l,int m,int n)
{
   cerr << "ERROR : INDV(k,l,m,n) Not known!!!!!!" << endl;
   exit(-1);
}

inline NiTi15TPPLat::interaction NiTi15TPPLat::INTER(int p,int q)
{
   interaction ans;
   switch (p)
   {
      case 0:
      case 1:
	 switch (q)
	 {
	    case 0:
	    case 1:
	       ans = aa;
	       break;
	    case 2:
	    case 3:
	       ans = ab;
	       break;
	 }
	 break;
      case 2:
      case 3:
	 switch (q)
	 {
	    case 0:
	    case 1:
	       ans = ab;
	       break;
	    case 2:
	    case 3:
	       ans = bb;
	       break;
	 }
	 break;
   }

   return ans;
}

Matrix NiTi15TPPLat::Phi(unsigned moduliflag,RadiiMorse::YDeriv dy,
			 RadiiMorse::TDeriv dt)
{
   static Matrix U(3,3);
   static Matrix V(4,3);
   static Matrix Eigvals(1,3);
   static double X[3];
   static Vector DX(DIM3),Dx(DIM3);
   static Vector Direction(DIM3);
   static double ForceNorm;
   static double J;
   static int i,j,k,l,p,q,m,n,s,z;
   static double r2,phi,phi1,phi2,Influancedist[DIM3],tmp;
   static int Top[DIM3],Bottom[DIM3],CurrentInfluanceDist;
   static interaction Inter;
   static Matrix Phi;

   switch (dy)
   {
      case RadiiMorse::Y0:
	 Phi.Resize(1,1,0.0);
	 break;
      case RadiiMorse::DY:
	 Phi.Resize(1,15,0.0);
	 break;
      case RadiiMorse::D2Y:
	 Phi.Resize(15,15,0.0);
	 break;
      case RadiiMorse::D3Y:
	 Phi.Resize(175,15,0.0);
	 break;
      case RadiiMorse::D4Y:
	 Phi.Resize(175,175,0.0);
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
      for (z=0;z<INTERNAL_ATOMS;z++)
      {
	 BodyForce_[z][p] = 0.0;
      }
   }
   // misc initialization
   ForceNorm = 0.0;

   for (p=0;p<INTERNAL_ATOMS;p++)
   {
      for (q=0;q<INTERNAL_ATOMS;q++)
      {
	 Inter = INTER(p,q);
	 
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
			DX[i] += (X[j] + A[q][j] - A[p][j] + V[q][j] - V[p][j])
			   *LatticeVec_[j][i]*RefLen_;
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
		  phi1 = Potential_[Inter].PairPotential(NTemp_,r2,RadiiMorse::DY,
							 RadiiMorse::T0);
		  if (ForceNorm < fabs(-phi1/(2.0*sqrt(r2))))
		  {
		     ForceNorm = fabs(-phi1/(2.0*sqrt(r2)));
		  }
		  for (i=0;i<DIM3;i++)
		  {
		     BodyForce_[p][i] += -phi1*Dx[i]/(2.0*sqrt(r2));
		  }
		  
		  // Calculate Phi
		  phi = Potential_[Inter].PairPotential(NTemp_,r2,dy,dt);
		  switch (dy)
		  {
		     case RadiiMorse::Y0:
			Phi[0][0]+=phi;
			break;
		     case RadiiMorse::DY:
			for (i=0;i<DIM3;i++)
			{
			   for (j=i;j<DIM3;j++)
			   {
			      Phi[0][INDU(i,j)] += phi*PI(Dx,DX,i,j);
			   }
			}
			for (i=1;i<4;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      Phi[0][INDV(i,j)] += phi*OMEGA(Dx,p,q,i,j);
			   }
			}
			break;
		     case RadiiMorse::D2Y:
			phi1=Potential_[Inter].PairPotential(NTemp_,r2,RadiiMorse::DY,
							     dt);

			//Upper Diag Block (6,6)
			for (i=0;i<DIM3;i++)
			{
			   for (j=i;j<DIM3;j++)
			   {
			      for (k=0;k<DIM3;k++)
			      {
				 for (l=k;l<DIM3;l++)
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
			      for (k=1;k<4;k++)
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
			      for (k=1;k<4;k++)
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
		     case RadiiMorse::D3Y:
			// phi1=PairPotential(Inter,r2,D2Y,dt);
			// 
			// for (i=0;i<DIM3;i++)
			//    for (j=i;j<DIM3;j++)
			// 	 for (k=0;k<DIM3;k++)
			// 	    for (l=k;l<DIM3;l++)
			// 	       for (m=0;m<DIM3;m++)
			// 		  for (n=m;n<DIM3;n++)
			// 		  {
			// 		     Phi[INDU(i,j,k,l)][INDU(m,n)] +=
			// 			phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PI(Dx,DX,m,n)
			// 			+0.5*phi1*(PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
			// 				   PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
			// 				   PI(Dx,DX,m,n)*PSI(DX,i,j,k,l));
			// 		  }
			cerr << "D3Y not programmed !" << endl;
			exit(1);
			break;
		     case RadiiMorse::D4Y:
			// phi1=PairPotential(Inter,r2,D3Y,dt);
			// phi2=PairPotential(Inter,r2,D2Y,dt);
			// 
			// for (i=0;i<DIM3;i++)
			//    for (j=i;j<DIM3;j++)
			// 	 for (k=0;k<DIM3;k++)
			// 	    for (l=k;l<DIM3;l++)
			// 	       for (m=0;m<DIM3;m++)
			// 		  for (n=m;n<DIM3;n++)
			// 		     for (p=0;p<DIM3;p++)
			// 			for (q=p;q<DIM3;q++)
			// 			{
			// 			   Phi[INDU(i,j,k,l)][INDU(m,n,p,q)]+=
			// 			      phi*(PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*
			// 				   PI(Dx,DX,m,n)*PI(Dx,DX,p,q)) +
			// 			      0.5*phi1*(0.5*(
			// 				 PI(Dx,DX,p,q)*(
			// 				    PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
			// 				    PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
			// 				    PI(Dx,DX,m,n)*PSI(DX,i,j,k,l)) +
			// 				 PI(Dx,DX,m,n)*(
			// 				    PI(Dx,DX,i,j)*PSI(DX,k,l,p,q) +
			// 				    PI(Dx,DX,k,l)*PSI(DX,i,j,p,q) +
			// 				    PI(Dx,DX,p,q)*PSI(DX,i,j,k,l)) +
			// 				 PI(Dx,DX,k,l)*(
			// 				    PI(Dx,DX,i,j)*PSI(DX,p,q,m,n) +
			// 				    PI(Dx,DX,p,q)*PSI(DX,i,j,m,n) +
			// 				    PI(Dx,DX,m,n)*PSI(DX,i,j,p,q)) +
			// 				 PI(Dx,DX,i,j)*(
			// 				    PI(Dx,DX,p,q)*PSI(DX,k,l,m,n) +
			// 				    PI(Dx,DX,k,l)*PSI(DX,p,q,m,n) +
			// 				    PI(Dx,DX,m,n)*PSI(DX,p,q,k,l))
			// 				 )) +
			// 			      0.25*phi2*(
			// 				 PSI(DX,i,j,m,n)*PSI(DX,k,l,p,q) +
			// 				 PSI(DX,k,l,m,n)*PSI(DX,i,j,p,q) +
			// 				 PSI(DX,i,j,k,l)*PSI(DX,m,n,p,q));
			// 			}
			cerr << "D4Y not programmed !" << endl;
			exit(1);
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
   Phi *= 1.0/(2.0*
	       (2.0*RefLen_*RefLen_*RefLen_*
		(LatticeVec_[0]%LatticeVec_[1])*LatticeVec_[2])
	       *ShearMod_);

   if (moduliflag)
   {
      return Phi;
   }

   // Add in External Work Terms
   // Recall that Pressure_ is applied stress
   // Thus E = W - pJ

   if (dt == RadiiMorse::T0)
   {
      switch (dy)
      {
	 case RadiiMorse::Y0:
	    Phi[0][0] = Phi[0][0] - Pressure_*J;
	    break;
	 case RadiiMorse::DY:
	 {
	    Matrix Uinv = U.Inverse();
	    for (i=0;i<DIM3;i++)
	       for (j=i;j<DIM3;j++)
	       {
		  Phi[0][INDU(i,j)] -= Pressure_*J*Uinv[i][j];
	       }
	 }
	 break;
	 case RadiiMorse::D2Y:
	    for (i=0;i<DIM3;i++)
	       for (j=i;j<DIM3;j++)
		  for (k=0;k<DIM3;k++)
		     for (l=k;l<DIM3;l++)
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
	 case RadiiMorse::D3Y:
	    for (i=0;i<DIM3;i++)
	       for (j=i;j<DIM3;j++)
		  for (k=0;k<DIM3;k++)
		     for (l=k;l<DIM3;l++)
			for (q=0;q<DIM3;q++)
			   for (s=q;s<DIM3;s++)
			   {
			      Phi[INDU(i,j,k,l)][INDU(q,s)] -=
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
	 case RadiiMorse::D4Y:
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
   return Phi(0,RadiiMorse::DY);
}

Matrix NiTi15TPPLat::StressDT()
{
   return Phi(0,RadiiMorse::DY,RadiiMorse::DT);
}

Matrix NiTi15TPPLat::Stiffness()
{
   return Phi(0,RadiiMorse::D2Y);
}

Matrix NiTi15TPPLat::Moduli()
{
   return Phi(1,RadiiMorse::D2Y);
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

void NiTi15TPPLat::CriticalPointInfo(int Width,ostream &out)
{
   // Matrix L4 = Phi(0,D4Y,T0),
   // 	 L3 = Phi(0,D3Y,T0),
   // 	 L2T = Phi(0,D2Y,DT),
   // 	 L = Stiffness();
   // 
   // double E,Ehat,Et1,Et2,eta,etahat,a,b,c;
   // 
   // E = 8.0*L3[INDU(0,1,0,2)][INDU(1,2)];
   // 
   // // Determine eta and etahat
   // a = L[0][0]*L3[INDU(2,2,0,1)][INDU(0,1)] +
   // 	 L[0][1]*(L3[INDU(2,2,0,1)][INDU(0,1)] - 2.0*L3[INDU(0,0,0,1)][INDU(0,1)]);
   // b = (L[0][0] - L[0][1])*(L[0][0] + 2.0*L[0][1]);
   // c = L[0][0]*L3[INDU(0,0,0,1)][INDU(0,1)] - L[0][1]*L3[INDU(2,2,0,1)][INDU(0,1)];
   // 
   // eta = -4.0*a/b;
   // etahat = -4.0*c/b;
   // //
   // 
   // Ehat = 16.0*L4[INDU(0,1,0,1)][INDU(0,1,0,1)] +
   // 	 12.0*(eta*L3[INDU(2,2,0,1)][INDU(0,1)] + 2.0*etahat*L3[INDU(0,0,0,1)][INDU(0,1)]);
   // 
   // Et1 = 4.0*(L3[INDU(0,1,0,1)][INDU(2,2)] + 2.0*L3[INDU(0,1,0,1)][INDU(0,0)]);
   // Et2 = 4.0*L2T[INDU(0,1)][INDU(0,1)];
   // 
   // // Print out results
   // for (int i=0;i<70;i++)
   // {
   // 	 cout << "=";
   // 	 out << "=";
   // }
   // cout << endl; out << endl;
   // cout << "Asymptotic Results for Principal Branch" << endl;
   // out << "Asymptotic Results for Principal Branch" << endl;
   // cout << "E =" << setw(Width) << E << endl;
   // out << "E =" << setw(Width) << E << endl;
   // cout << "Ehat =" << setw(Width) << Ehat << endl;
   // out << "Ehat =" << setw(Width) << Ehat << endl;
   // cout << "Etheta_1 =" << setw(Width) << Et1 << endl;
   // out << "Etheta_1 =" << setw(Width) << Et1 << endl;
   // cout << "Etheta_2 =" << setw(Width) << Et2 << endl;
   // out << "Etheta_2 =" << setw(Width) << Et2 << endl;
   // for (int i=0;i<70;i++)
   // {
   // 	 cout << "-";
   // 	 out << "-";
   // }
   // cout << endl; out << endl;
   // cout << "Etheta = Etheta_1 * dlamda/dtheta + Etheta_2" << endl;
   // out << "Etheta = Etheta_1 * dlamda/dtheta + Etheta_2" << endl;
   // cout << "Rhombohedral Tangent ==> -E/Etheta" << endl;
   // out << "Rhombohedral Tangent ==> -E/Etheta" << endl;
   // cout << "Orthorhombic Curveature ==> -Ehat/(3*Etheta)" << endl;
   // out << "Orthorhombic Curveature ==> -Ehat/(3*Etheta)" << endl;
   // for (int i=0;i<70;i++)
   // {
   // 	 cout << "=";
   // 	 out << "=";
   // }
   // cout << endl; out << endl;
   // 
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
