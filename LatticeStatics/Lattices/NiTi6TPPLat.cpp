#include "NiTi6TPPLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTi6TPPLat::NiTi6TPPLat(char *datafile)
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

int NiTi6TPPLat::FindLatticeSpacing(int iter,double dx)
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

double NiTi6TPPLat::PI(const Vector &Dx,const Vector &DX,
			int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double NiTi6TPPLat::PSI(const Vector &DX,
			 int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}

double NiTi6TPPLat::pwr(const double &x,const unsigned y)
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

inline int NiTi6TPPLat::INDU(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

inline int NiTi6TPPLat::INDUU(int k,int l,int m,int n)
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

Matrix NiTi6TPPLat::Phi(unsigned moduliflag,PairPotentials::YDeriv dy,
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
			DX[i] += (X[j] + A[q][j] - A[p][j])*
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

void NiTi6TPPLat::Print(ostream &out,PrintDetail flag)
{
   int W=out.width();

   out.width(0);
   cout.width(0);

   double MinEigVal;
   int NoNegEigVal=0;

   double energy=Energy();
   
   Matrix
      stress = Stress(),
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

   int RankOneConvex = FullScanRank1Convex3D(moduli,ConvexityDX_);

   switch (flag)
   {
      case PrintLong:
	 out << "NiTi6TPPLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	     << "Potential Parameters : " << endl
	     << "AA -- " << setw(W) << Potential_[aa] << endl
	     << "BB -- " << setw(W) << Potential_[bb] << endl
	     << "AB -- " << setw(W) << Potential_[ab] << endl
	     << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 cout << "NiTi6TPPLat:" << endl << endl
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
	     << "Potential Value (G Normalized):" << setw(W) << energy << endl
	     << "BodyForce Value 0 (Inf Normalized):" << setw(W) << BodyForce_[0] << endl
	     << "BodyForce Value 1 (Inf Normalized):" << setw(W) << BodyForce_[1] << endl
	     << "Stress (G Normalized):" << setw(W) << stress << endl
	     << "Stiffness (G Normalized):" << setw(W) << stiffness
	     << "Rank 1 Convex:" << setw(W) << RankOneConvex << endl
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl;
	 cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	      << "Pressure (G Normalized): " << setw(W) << Pressure_ << endl
	      << "DOF's :" << endl << setw(W) << DOF_ << endl
	      << "Potential Value (G Normalized):" << setw(W) << energy << endl
	      << "BodyForce Value 0 (Inf Normalized):" << setw(W) << BodyForce_[0] << endl
	      << "BodyForce Value 1 (Inf Normalized):" << setw(W) << BodyForce_[1] << endl
	      << "Stress (G Normalized):" << setw(W) << stress << endl
	      << "Stiffness (G Normalized):" << setw(W) << stiffness
	      << "Rank 1 Convex:" << setw(W) << RankOneConvex << endl
	      << "Eigenvalue Info:"  << setw(W) << EigenValues
	      << "Bifurcation Info:" << setw(W) << MinEigVal
	      << setw(W) << NoNegEigVal << endl;	 
	 break;
   }
}

ostream &operator<<(ostream &out,NiTi6TPPLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double NiTi6TPPLat::Alt[DIM3][DIM3][DIM3] = {0.0, 0.0, 0.0,
						    0.0, 0.0, 1.0,
						    0.0, -1.0,0.0,
						    0.0, 0.0, -1.0,
						    0.0, 0.0, 0.0,
						    1.0, 0.0, 0.0,
						    0.0, 1.0, 0.0,
						    -1.0,0.0, 0.0,
						    0.0, 0.0, 0.0};

const double NiTi6TPPLat::A[INTERNAL_ATOMS][DIM3] = {0.0, 0.0, 0.0,
						     0.5, 0.5, 0.5};

const NiTi6TPPLat::interaction NiTi6TPPLat::INTER[INTERNAL_ATOMS][INTERNAL_ATOMS] = {
   aa, ab,
   ab, bb};
