#include "NiTiShuffle1TPPLat.h"
#include <math.h>

#include "UtilityFunctions.h"

NiTiShuffle1TPPLat::NiTiShuffle1TPPLat(char *datafile)
{
   // First Size DOF
   DOF_.Resize(7);
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

int NiTiShuffle1TPPLat::FindLatticeSpacing(int iter,double dx)
{
   double oldPressure=Pressure_,
      oldTemp=NTemp_;

   Pressure_=0.0;
   NTemp_=1.0;
   ShearMod_=1.0;
   DOF_[0] = DOF_[1] = DOF_[2] = 1.0;
   DOF_[3] = DOF_[4] = DOF_[5] = DOF_[6] = 0.0;

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

void NiTiShuffle1TPPLat::GetLatticeVectorInfo(double *SX,double *DXPrime,
					     interaction &Inter,int p,int q)
{
   static double Basis[INTERNAL_ATOMS][DIM3] = {0.0,0.0,0.0,
						0.5,0.5,0.0,
						0.5,0.0,0.5,
						0.0,0.5,0.5};
   double dxprime=0;

   Basis[1][0] = 0.5 + DOF_[6];
   Basis[3][0] = DOF_[6];

   switch (p)
   {
      case 0:
	 switch (q)
	 {
	    case 0:
	       dxprime = 0.0;
	       Inter = aa;
	       break;
	    case 1:
	       dxprime = 1.0;
	       Inter = aa;
	       break;
	    case 2:
	       dxprime = 0.0;
	       Inter = ab;
	       break;
	    case 3:
	       dxprime = 1.0;
	       Inter = ab;
	       break;
	 }
	 break;
      case 1:
	 switch (q)
	 {
	    case 0:
	       dxprime = -1.0;
	       Inter = aa;
	       break;
	    case 1:
	       dxprime = 0.0;
	       Inter = aa;
	       break;
	    case 2:
	       dxprime = -1.0;
	       Inter = ab;
	       break;
	    case 3:
	       dxprime = 0.0;
	       Inter = ab;
	       break;
	 }
	 break;
      case 2:
	 switch (q)
	 {
	    case 0:
	       dxprime = 0.0;
	       Inter = ab;
	       break;
	    case 1:
	       dxprime = 1.0;
	       Inter = ab;
	       break;
	    case 2:
	       dxprime = 0.0;
	       Inter = bb;
	       break;
	    case 3:
	       dxprime = 1.0;
	       Inter = bb;
	       break;
	 }
	 break;
      case 3:
	 switch (q)
	 {
	    case 0:
	       dxprime = -1.0;
	       Inter = ab;
	       break;
	    case 1:
	       dxprime = 0.0;
	       Inter = ab;
	       break;
	    case 2:
	       dxprime = -1.0;
	       Inter = bb;
	       break;
	    case 3:
	       dxprime = 0.0;
	       Inter = bb;
	       break;
	 }
   }

   for (int i=0;i<DIM3;i++)
   {
      SX[i] = Basis[q][i] - Basis[p][i];
      DXPrime[i] = 0.0;
   }
   DXPrime[0] = dxprime;

}

inline double NiTiShuffle1TPPLat::PI(const Vector &Dx,const Vector &DX,
				    int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

inline double NiTiShuffle1TPPLat::PSI(const Vector &DX,
				     int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}


double NiTiShuffle1TPPLat::pwr(const double &x,const unsigned y)
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

inline int NiTiShuffle1TPPLat::IND(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

inline int NiTiShuffle1TPPLat::IND(int k,int l,int m,int n)
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

Matrix NiTiShuffle1TPPLat::Phi(unsigned moduliflag,PairPotentials::YDeriv dy,PairPotentials::TDeriv dt)
{
   static Matrix U(3,3);
   static Matrix Eigvals(1,3);
   static double X[3],SX[3],DXP[3];
   static Vector DX(DIM3),Dx(DIM3);
   static Vector DXPrime(DIM3),DxPrime(DIM3),Direction(DIM3);
   static double ForceNorm;
   static double J;
   static int i,j,k,l,p,q,m,n,s,z;
   static double r2,phi,phi1,phi2,Influancedist[DIM3],tmp;
   static int Top[DIM3],Bottom[DIM3],CurrentInfluanceDist;
   static interaction Inter;
   static Matrix Phi;

   switch (dy)
   {
      case PairPotentials::Y0:
	 Phi.Resize(1,1,0.0);
	 break;
      case PairPotentials::DY:
	 Phi.Resize(1,7,0.0);
	 break;
      case PairPotentials::D2Y:
	 Phi.Resize(7,7,0.0);
	 break;
      case PairPotentials::D3Y:
	 Phi.Resize(49,7,0.0);
	 break;
      case PairPotentials::D4Y:
	 Phi.Resize(49,49,0.0);
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
      for (z=0;z<INTERNAL_ATOMS;z++)
	 BodyForce_[z][p] = 0.0;
   }
   // misc initialization
   ForceNorm = 0.0;

   for (p=0;p<INTERNAL_ATOMS;p++)
   {
      for (q=0;q<INTERNAL_ATOMS;q++)
      {
	 GetLatticeVectorInfo(SX,DXP,Inter,p,q);
	 for (X[0] = Bottom[0];X[0] <= Top[0];X[0]++)
	 {
	    for (X[1] = Bottom[1];X[1] <= Top[1];X[1]++)
	    {
	       for (X[2] = Bottom[2];X[2] <= Top[2];X[2]++)
	       {
		  for (i=0;i<DIM3;i++)
		  {
		     DX[i] = 0.0;
		     DXPrime[i] = 0.0;
		     
		     for (j=0;j<DIM3;j++)
		     {
			DX[i] += (X[j] + SX[j])*LatticeVec_[j][i]*RefLen_;
			DXPrime[i] += DXP[j]*RefLen_*LatticeVec_[j][i];
		     }
		  }

		  r2 = 0.0;
		  for (i=0;i<DIM3;i++)
		  {
		     Dx[i] = 0.0;
		     DxPrime[i] = 0.0;

		     for (j=0;j<DIM3;j++)
		     {
			Dx[i] += U[i][j] * DX[j];
			DxPrime[i] += U[i][j] * DXPrime[j];
		     }
		     r2 += Dx[i]*Dx[i];
		  }
		  // Only use Sphere of Influance (current)
		  if (r2==0 || r2 > InfluanceDist_*InfluanceDist_)
		  {
		     continue;
		  }

		  // Calculate bodyforce
		  phi1 = 2.0*sqrt(r2)*
		     Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,PairPotentials::T0);
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
			      Phi[0][IND(i,j)] += phi*PI(Dx,DX,i,j);
			   }
			   Phi[0][6] += phi*(2.0*DxPrime[i]*Dx[i]);
			}
			break;
		     case PairPotentials::D2Y:
			phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::DY,dt);
		     
			for (i=0;i<DIM3;i++)
			{
			   for (j=0;j<DIM3;j++)
			   {
			      for (k=0;k<DIM3;k++)
			      {
				 for (l=0;l<DIM3;l++)
				 {
				    Phi[IND(i,j)][IND(k,l)]+=
				       phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)
				       +phi1*(0.5)*PSI(DX,i,j,k,l);
				 }
			      }
			      Phi[6][IND(i,j)] = Phi[IND(i,j)][6]
				 += phi*(PI(Dx,DX,i,j)*(2.0*DxPrime*Dx))
				 + phi1*(DxPrime[i]*DX[j] + Dx[i]*DXPrime[j] +
					       DXPrime[i]*Dx[j] + DX[i]*DxPrime[j]);
			   }
			}
			Phi[6][6] += phi*(2.0*(DxPrime*Dx))*(2.0*(DxPrime*Dx))
			   + phi1*(2.0*DxPrime*DxPrime);
			break;
		     case PairPotentials::D3Y:
			// phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y,dt);
			// 
			// for (i=0;i<DIM3;i++)
			//    for (j=0;j<DIM3;j++)
			// 	 for (k=0;k<DIM3;k++)
			// 	    for (l=0;l<DIM3;l++)
			// 	       for (m=0;m<DIM3;m++)
			// 		  for (n=0;n<DIM3;n++)
			// 		  {
			// 		     Phi[IND(i,j,k,l)][IND(m,n)] +=
			// 			phi*PI(Dx,DX,i,j)*PI(Dx,DX,k,l)*PI(Dx,DX,m,n)
			// 			+0.5*phi1*(PI(Dx,DX,i,j)*PSI(DX,k,l,m,n) +
			// 				   PI(Dx,DX,k,l)*PSI(DX,i,j,m,n) +
			// 				   PI(Dx,DX,m,n)*PSI(DX,i,j,k,l));
			// 		  }
			cerr << "D3Y not programmed !" << endl;
			exit(1);
			break;
		     case PairPotentials::D4Y:
			// phi1=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D3Y,dt);
			// phi2=Potential_[Inter].PairPotential(NTemp_,r2,PairPotentials::D2Y,dt);
			// 
			// for (i=0;i<DIM3;i++)
			//    for (j=0;j<DIM3;j++)
			// 	 for (k=0;k<DIM3;k++)
			// 	    for (l=0;l<DIM3;l++)
			// 	       for (m=0;m<DIM3;m++)
			// 		  for (n=0;n<DIM3;n++)
			// 		     for (p=0;p<DIM3;p++)
			// 			for (q=0;q<DIM3;q++)
			// 			{
			// 			   Phi[IND(i,j,k,l)][IND(m,n,p,q)]+=
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
   Phi *= 1.0/(2.0*(2.0*RefLen_*RefLen_*RefLen_)*ShearMod_);

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
		  Phi[0][IND(i,j)] -= Pressure_*J*Uinv[i][j];
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
			      Phi[IND(i,j)][IND(k,l)] -=
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
			      Phi[IND(i,j,k,l)][IND(q,s)] -=
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

void NiTiShuffle1TPPLat::Print(ostream &out,PrintDetail flag)
{
   int W=out.width();

   out.width(0);

   double MinEigVal;
   int NoNegEigVal=0;
   
   Matrix
      stiffness = Stiffness(),
      moduli = Moduli(),
      EigenValues(1,7);

   EigenValues=SymEigVal(stiffness);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<7;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   switch (flag)
   {
      case PrintLong:
	 out << "NiTiShuffle1TPPLat:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_ << endl
	     << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl
	     << "Potential Parameters : "
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

ostream &operator<<(ostream &out,NiTiShuffle1TPPLat &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}

const double NiTiShuffle1TPPLat::Alt[DIM3][DIM3][DIM3]= {0.0, 0.0, 0.0,
							 0.0, 0.0, 1.0,
							 0.0, -1.0,0.0,
							 0.0, 0.0, -1.0,
							 0.0, 0.0, 0.0,
							 1.0, 0.0, 0.0,
							 0.0, 1.0, 0.0,
							 -1.0,0.0, 0.0,
							 0.0, 0.0, 0.0};
