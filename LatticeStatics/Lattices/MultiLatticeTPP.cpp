#include "MultiLatticeTPP.h"
#include <math.h>

#include "UtilityFunctions.h"

MultiLatticeTPP::~MultiLatticeTPP()
{
   delete [] BodyForce_;
   delete [] AtomicMass_;
   for (int i=0;i<INTERNAL_ATOMS;++i)
      for (int j=i;j<INTERNAL_ATOMS;++j)
	 delete Potential_[i][j];
   delete [] Potential_[0];
   delete [] Potential_;
   delete [] A_;
}

MultiLatticeTPP::MultiLatticeTPP(char *datafile)
{
   // Get Lattice definition
   char tmp[LINELENGTH];
   GetParameter("^InternalAtoms",datafile,"%u",&INTERNAL_ATOMS);
   DOFS = 6 + 3*(INTERNAL_ATOMS-1);
   for (int i=0;i<DIM3;++i)
   {
      LatticeBasis[i].Resize(DIM3);
      sprintf(tmp,"^LatticeBasis_%u",i);
      GetVectorParameter(tmp,datafile,&(LatticeBasis[i]));
   }
   
   // First Size DOF
   DOF_.Resize(DOFS,0.0);
   DOF_[0]=DOF_[1]=DOF_[2] = 1.0;
   // Set RefLattice_
   RefLattice_.Resize(DIM3,DIM3);
   // Set A
   A_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      A_[i].Resize(DIM3);
      sprintf(tmp,"^AtomPosition_%u",i);
      GetVectorParameter(tmp,datafile,&(A_[i]));
   }
   
   // Setup Bodyforce_
   BodyForce_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
      BodyForce_[i].Resize(DIM3,0.0);

   // Get Potential Parameters
   Potential_ = new PairPotentials**[INTERNAL_ATOMS];
   Potential_[0] = new PairPotentials*[INTERNAL_ATOMS*INTERNAL_ATOMS];
   for (int i=1;i<INTERNAL_ATOMS;++i)
   {
      Potential_[i] = Potential_[i-1] + INTERNAL_ATOMS;
   }

   AtomicMass_ = new double[INTERNAL_ATOMS];

   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      for (int j=i;j<INTERNAL_ATOMS;++j)
      {
	 Potential_[i][j] = Potential_[j][i] = InitializePairPotential(datafile,i,j);
      }

      sprintf(tmp,"^AtomicMass_%u",i);
      GetParameter(tmp,datafile,"%lf",&(AtomicMass_[i]));
   }

	 
   // Get Lattice parameters
   
   // all 3 components of RefLen could be set by the input file
   // but for now just take 1
   GetParameter("^RefLen",datafile,"%lf",&(RefLen_[0]));
   RefLen_[2] = RefLen_[1] = RefLen_[0];
   
   GetParameter("^InfluanceDist",datafile,"%u",&InfluanceDist_);
   GetParameter("^NTemp",datafile,"%lf",&NTemp_);
   GetParameter("^Pressure",datafile,"%lf",&Pressure_);
   GetParameter("^ConvexityDX",datafile,"%lf",&ConvexityDX_);
   
   // needed to initialize reference length
   int iter;
   double DX;
   GetParameter("^MaxIterations",datafile,"%u",&iter);
   GetParameter("^InitializeStepSize",datafile,"%lf",&DX);
   GetParameter("^BlochWaveGridSize",datafile,"%u",&GridSize_);

   // Initialize RefLattice_
   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
	 RefLattice_[i][j] = LatticeBasis[i][j]*RefLen_[i];

   // Initiate the Lattice Sum object
   LatSum_(&DOF_,&RefLattice_,INTERNAL_ATOMS,A_,&InfluanceDist_);

   int err=0;
   err=FindLatticeSpacing(iter,DX);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << endl;
      exit(-1);
   }

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   UCIter_(GridSize_,&RefLattice_);
   
}

// This routing assumes Cubic symm, i.e., RefLen_[0--2] are equal
int MultiLatticeTPP::FindLatticeSpacing(int iter,double dx)
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

   LatSum_.Recalc();
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
	 cout << setw(20) << RefLen_[0] << setw(20) << s << endl;
	 RefLen_[0] -= dx*sign;
	 // Update RefLattice_
	 for (int p=0;p<DIM3;++p)
	    for (int q=0;q<DIM3;++q)
	       RefLattice_[p][q] = LatticeBasis[p][q]*RefLen_[0];
	 
	 LatSum_.Recalc();
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
	    cout << setw(20) << RefLen_[0] << setw(20) << s << endl;
	    i++;
	    RefLen_[0] -= dx*newsign/(pow(2,i));
	    // Update RefLattice_
	    for (int p=0;p<DIM3;++p)
	       for (int q=0;q<DIM3;++q)
		  RefLattice_[p][q] = LatticeBasis[p][q]*RefLen_[0];

	    LatSum_.Recalc();
	    s = Stress()[0][0];
	    newsign=s/fabs(s);
	 }
	 if (i > iter)
	    return 1;
      }
   }

   RefLen_[2]=RefLen_[1] = RefLen_[0];
   ShearMod_ = 0.25*fabs(Stiffness()[5][5]);
   NTemp_=oldTemp;
   Pressure_=oldPressure;

   return 0;
}

// Lattice Routines

double MultiLatticeTPP::PI(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double MultiLatticeTPP::PSI(double *DX,int r, int s, int t, int u)
{
   return (Del(r,t)*DX[s]*DX[u] +
	   Del(r,u)*DX[s]*DX[t] +
	   Del(s,t)*DX[r]*DX[u] +
	   Del(s,u)*DX[r]*DX[t]);
}

double MultiLatticeTPP::OMEGA(double *Dx,int p,int q,int i, int j)
{
   double ret=0;
   
   ret=0;
   for (int s=0;s<DIM3;s++)
   {
      for (int t=0;t<DIM3;t++)
      {
	 ret += (LatticeBasis[j][s]*DOF_[INDU(s,t)]*Dx[t] +
		 Dx[t]*DOF_[INDU(t,s)]*LatticeBasis[j][s]);
      }
   }
   ret *= DELTA(i,p,q);

   return ret;
}

double MultiLatticeTPP::SIGMA(int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   for (int s=0;s<DIM3;s++)
   {
      for (int t=0;t<DIM3;t++)
      {
	 for (int r=0;r<DIM3;r++)
	 {
	    tmp += (LatticeBasis[j][s]*DOF_[INDU(s,t)]*DOF_[INDU(t,r)]*LatticeBasis[l][r] +
		    LatticeBasis[l][s]*DOF_[INDU(s,t)]*DOF_[INDU(t,r)]*LatticeBasis[j][r]);
	 }
      }
   }
   
   return DELTA(i,p,q)*DELTA(k,p,q)*tmp;
}

double MultiLatticeTPP::GAMMA(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   
   for (int s=0;s<DIM3;s++)
   {
      tmp += (LatticeBasis[l][s]*DOF_[INDU(s,i)]*DX[j] +
	      LatticeBasis[l][s]*DOF_[INDU(s,j)]*DX[i] +
	      DX[i]*DOF_[INDU(j,s)]*LatticeBasis[l][s] +
	      DX[j]*DOF_[INDU(i,s)]*LatticeBasis[l][s]);
   }
   
   
   return (0.5*DELTA(k,p,q)*(LatticeBasis[l][i]*Dx[j] + LatticeBasis[l][j]*Dx[i] + tmp));
}

double MultiLatticeTPP::THETA(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return (0.5*DELTA(m,p,q)*(Del(i,k)*LatticeBasis[n][j]*DX[l]
			     + Del(i,k)*DX[j]*LatticeBasis[n][l]
			     + Del(i,l)*LatticeBasis[n][j]*DX[k]
			     + Del(i,l)*DX[j]*LatticeBasis[n][k]
			     + Del(j,k)*LatticeBasis[n][i]*DX[l]
			     + Del(j,k)*DX[i]*LatticeBasis[n][l]
			     + Del(j,l)*LatticeBasis[n][i]*DX[k]
			     + Del(j,l)*DX[i]*LatticeBasis[n][k]));
}

double MultiLatticeTPP::XI(int p,int q,int i,int j,int k,int l,int m,int n)
{
   double tmp=0;

   for (int s=0;s<DIM3;s++)
   {
      tmp += (LatticeBasis[j][m]*DOF_[INDU(n,s)]*LatticeBasis[l][s]
	      + LatticeBasis[j][n]*DOF_[INDU(m,s)]*LatticeBasis[l][s]
	      + LatticeBasis[j][s]*DOF_[INDU(s,m)]*LatticeBasis[l][n]
	      + LatticeBasis[j][s]*DOF_[INDU(s,n)]*LatticeBasis[l][m]
	      + LatticeBasis[l][m]*DOF_[INDU(n,s)]*LatticeBasis[j][s]
	      + LatticeBasis[l][n]*DOF_[INDU(m,s)]*LatticeBasis[j][s]
	      + LatticeBasis[l][s]*DOF_[INDU(s,m)]*LatticeBasis[j][n]
	      + LatticeBasis[l][s]*DOF_[INDU(s,n)]*LatticeBasis[j][m]);
   }

   return (0.5*DELTA(i,p,q)*DELTA(k,p,q)*tmp);
}

double MultiLatticeTPP::LAMDA(int p,int q,int i,int j,int k,int l,int m,int n,int a,
			  int b)
{
   return (0.5*DELTA(m,p,q)*DELTA(a,p,q)*
	   (Del(i,k)*LatticeBasis[n][j]*LatticeBasis[b][l]
	    + Del(i,k)*LatticeBasis[b][j]*LatticeBasis[n][l]
	    + Del(i,l)*LatticeBasis[n][j]*LatticeBasis[b][k]
	    + Del(i,l)*LatticeBasis[b][j]*LatticeBasis[n][k]
	    + Del(j,k)*LatticeBasis[n][i]*LatticeBasis[b][l]
	    + Del(j,k)*LatticeBasis[b][i]*LatticeBasis[n][l]
	    + Del(j,l)*LatticeBasis[n][i]*LatticeBasis[b][k]
	    + Del(j,l)*LatticeBasis[b][i]*LatticeBasis[n][k]));
}

inline int MultiLatticeTPP::INDU(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

inline int MultiLatticeTPP::INDUU(int k,int l,int m,int n)
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

inline int MultiLatticeTPP::INDV(int i,int j)
{
   if (!i)
   {
      cerr << "Error: INDV(i,j) i==0!!!!!" << endl;
      exit(-1);
   }

   return 6 + (i-1)*3 + j;
}

inline int MultiLatticeTPP::INDVV(int k,int l,int m,int n)
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

inline int MultiLatticeTPP::INDUV(int i,int j,int m,int n)
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

inline int MultiLatticeTPP::INDVU(int m,int n,int i,int j)
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

double MultiLatticeTPP::Energy()
{
   double Phi = 0.0;

   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Phi
      Phi += Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::Y0,PairPotentials::T0);
   }

   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));

   Phi = Phi - Pressure_*LatSum_.J();;

   return Phi;
}

Matrix MultiLatticeTPP::stress(PairPotentials::TDeriv dt)
{
   static Matrix S;
   Matrix Uinv(DIM3,DIM3);
   double ForceNorm = 0.0;
   double phi,J;
   int i,j;

   S.Resize(1,DOFS,0.0);

   for (i=0;i<INTERNAL_ATOMS;++i)
   {
      for (j=0;j<DIM3;++j)
      {
	 BodyForce_[i][j] = 0.0;
      }
   }

   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate bodyforce
      // NOTE: phi1 = d(phi)/d(r2)
      // We need d(phi)/dr = 2*r*d(phi)/d(r2)
      phi = 2.0*sqrt(LatSum_.r2())*
	 Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	    NTemp_,LatSum_.r2(),PairPotentials::DY,PairPotentials::T0);
      if (ForceNorm < fabs(-phi/2.0))
      {
	 ForceNorm = fabs(-phi/2.0);
      }
      for (i=0;i<DIM3;i++)
      {
	 BodyForce_[LatSum_.Atom(0)][i] += -phi*LatSum_.Dx(i)/(2.0*sqrt(LatSum_.r2()));
      }

      // Claculate the Stress
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::DY,dt);
      for (i=0;i<DIM3;i++)
      {
	 for (j=0;j<DIM3;j++)
	 {
	    S[0][INDU(i,j)] += phi*PI(LatSum_.pDx(),LatSum_.pDX(),i,j);
	 }
      }
      for (i=1;i<INTERNAL_ATOMS;i++)
      {
	 for (j=0;j<DIM3;j++)
	 {
	    S[0][INDV(i,j)] += phi*OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					 LatSum_.Atom(1),i,j);
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
   
   // S = S/(2*Vr*ShearMod)
   S *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));

   // Pressure terms
   Uinv = LatSum_.UInv();
   J = LatSum_.J();
   for (i=0;i<DIM3;i++)
      for (j=0;j<DIM3;j++)
      {
	 S[0][INDU(i,j)] -= Pressure_*J*Uinv[i][j];
      }
   
   return S;
}
      
Matrix MultiLatticeTPP::stiffness(int moduliflag,PairPotentials::TDeriv dt)
{
   static Matrix Phi;
   Matrix U(DIM3,DIM3);
   double phi,phi1;
   int i,j,k,l,q,s;

   Phi.Resize(DOFS,DOFS,0.0);

   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D2Y,dt);
      phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::DY,dt);
      
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
		     phi*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
			  *PI(LatSum_.pDx(),LatSum_.pDX(),k,l))
		     +phi1*(0.5)*PSI(LatSum_.pDX(),i,j,k,l);
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
		     phi*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
			  *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l))
		     +phi1*SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
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
		     phi*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
			  *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l))
		     +phi1*GAMMA(LatSum_.pDx(),LatSum_.pDX(),
				 LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
	       }
	    }
	 }
      }
   }
      
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));
   
   if ((!moduliflag) && (dt == PairPotentials::T0))
   {
      U = LatSum_.U();
      // Pressure terms
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
			      Alt(i,k,q)*Alt(j,l,s) + Alt(j,k,q)*Alt(i,l,s) +
			      Alt(i,l,q)*Alt(j,k,s) + Alt(j,l,q)*Alt(i,k,s) +
			      Alt(i,q,k)*Alt(j,s,l) + Alt(j,q,k)*Alt(i,s,l) +
			      Alt(i,q,l)*Alt(j,s,k) + Alt(j,q,l)*Alt(i,s,k)
			      )*U[q][s];
		     }
	       }
   }
   
   return Phi;
}

Matrix MultiLatticeTPP::E3()
{
   static Matrix Phi;
   double phi,phi1,phi2;
   int i,j,k,l,m,n,q,s;

   Phi.Resize(DOFS*DOFS,DOFS,0.0);

   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D2Y,PairPotentials::T0);
      phi2=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::DY,PairPotentials::T0);
	 
      // DU^3 block
      for (i=0;i<DIM3;i++)
	 for (j=0;j<DIM3;j++)
	    for (k=0;k<DIM3;k++)
	       for (l=0;l<DIM3;l++)
		  for (m=0;m<DIM3;m++)
		     for (n=0;n<DIM3;n++)
		     {
			Phi[INDUU(i,j,k,l)][INDU(m,n)] +=
			   phi*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				*PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				*PI(LatSum_.pDx(),LatSum_.pDX(),m,n))
			   +0.5*phi1*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				      *PSI(LatSum_.pDX(),k,l,m,n) +
				      PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				      *PSI(LatSum_.pDX(),i,j,m,n) +
				      PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				      *PSI(LatSum_.pDX(),i,j,k,l));
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
			   phi*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
				*OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l)
				*OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
				       LatSum_.Atom(1),m,n))
			   +phi1*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					LatSum_.Atom(1),k,l)
				  *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				  + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					  LatSum_.Atom(1),i,j)
				  *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				  + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					  LatSum_.Atom(1),m,n)
				  *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l));
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
			      phi*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				   *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				   *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					  LatSum_.Atom(1),m,n))
			      +phi1*(PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				     *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					    LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				     + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				     *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					    LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				     + 0.5*OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
						 LatSum_.Atom(1),m,n)
				     *PSI(LatSum_.pDX(),i,j,k,l))
			      +phi2*THETA(LatSum_.pDX(),LatSum_.Atom(0),
					  LatSum_.Atom(1),i,j,k,l,m,n));
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
			      phi*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					 LatSum_.Atom(1),i,j)
				   *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					  LatSum_.Atom(1),k,l)
				   *PI(LatSum_.pDx(),LatSum_.pDX(),m,n))
			      +phi1*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),k,l)
				     *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					    LatSum_.Atom(0),LatSum_.Atom(1),m,n,i,j)
				     + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),i,j)
				     *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					    LatSum_.Atom(0),LatSum_.Atom(1),m,n,k,l)
				     + PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				     *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
			      +phi2*XI(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l,m,n));
		     }
   }

   
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));

   // Pressure terms
   for (i=0;i<DIM3;i++)
      for (j=0;j<DIM3;j++)
	 for (k=0;k<DIM3;k++)
	    for (l=0;l<DIM3;l++)
	       for (q=0;q<DIM3;q++)
		  for (s=0;s<DIM3;s++)
		  {
		     Phi[INDUU(i,j,k,l)][INDU(q,s)] -=
			(Pressure_/16.0)*(
			   Alt(i,k,q)*Alt(j,l,s) + Alt(j,k,q)*Alt(i,l,s) +
			   Alt(i,l,q)*Alt(j,k,s) + Alt(j,l,q)*Alt(i,k,s) +
			   Alt(i,q,k)*Alt(j,s,l) + Alt(j,q,k)*Alt(i,s,l) +
			   Alt(i,q,l)*Alt(j,s,k) + Alt(j,q,l)*Alt(i,s,k) +
			   Alt(i,k,s)*Alt(j,l,q) + Alt(j,k,s)*Alt(i,l,q) +
			   Alt(i,l,s)*Alt(j,k,q) + Alt(j,l,s)*Alt(i,k,q) +
			   Alt(i,s,k)*Alt(j,q,l) + Alt(j,s,k)*Alt(i,q,l) +
			   Alt(i,s,l)*Alt(j,q,k) + Alt(j,s,l)*Alt(i,q,k));
		  }

   return Phi;
}

Matrix MultiLatticeTPP::E4()
{
   static Matrix Phi;
   double phi,phi1,phi2,phi3;
   int i,j,k,l,m,n,s,t;

   Phi.Resize(DOFS*DOFS,DOFS*DOFS,0.0);

   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {  
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D4Y,PairPotentials::T0);
      phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi2=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::D2Y,PairPotentials::T0);
      phi3=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	 NTemp_,LatSum_.r2(),PairPotentials::DY,PairPotentials::T0);
      
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
				 phi*(PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				      *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				      *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				      *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)) +
				 0.5*phi1*(
				    PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				    *PSI(LatSum_.pDX(),i,j,s,t)
				    + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				    *PSI(LatSum_.pDX(),k,l,s,t)
				    + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				    *PSI(LatSum_.pDX(),m,n,s,t)
				    + PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				    *PSI(LatSum_.pDX(),i,j,m,n)
				    + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				    *PSI(LatSum_.pDX(),k,l,m,n)
				    + PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				    *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				    *PSI(LatSum_.pDX(),i,j,k,l))
				 +0.25*phi2*(
				    PSI(LatSum_.pDX(),i,j,m,n)
				    *PSI(LatSum_.pDX(),k,l,s,t)
				    + PSI(LatSum_.pDX(),i,j,s,t)
				    *PSI(LatSum_.pDX(),k,l,m,n)
				    + PSI(LatSum_.pDX(),i,j,k,l)
				    *PSI(LatSum_.pDX(),m,n,s,t));
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
				 phi*(OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),i,j)
				      *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),k,l)
				      *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),m,n)
				      *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),s,t))
				 +phi1*(
				    OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					  LatSum_.Atom(1),k,l)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),m,n)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
				    + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),i,j)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),m,n)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
				    + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),i,j)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),k,l)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
				    + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),k,l)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),s,t)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				    + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),i,j)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),s,t)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				    + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					    LatSum_.Atom(1),m,n)
				    *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					   LatSum_.Atom(1),s,t)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
				 +phi2*(
				    SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
				    + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				    + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l)
				    *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t));
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
				       PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),s,t))
				    +phi1*(
				       PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
				       +0.5*(
					  PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
					  *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
						 LatSum_.Atom(1),s,t)
					  *PSI(LatSum_.pDX(),i,j,m,n)
					  + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
					  *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
						 LatSum_.Atom(1),s,t)
					  *PSI(LatSum_.pDX(),k,l,m,n)
					  + PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
					  *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
						 LatSum_.Atom(1),s,t)
					  *PSI(LatSum_.pDX(),i,j,k,l)))
				    +phi2*(
				       PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *THETA(LatSum_.pDX(),LatSum_.Atom(0),
					      LatSum_.Atom(1),i,j,m,n,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *THETA(LatSum_.pDX(),LatSum_.Atom(0),
					      LatSum_.Atom(1),k,l,m,n,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),m,n)
				       *THETA(LatSum_.pDX(),LatSum_.Atom(0),
					      LatSum_.Atom(1),i,j,k,l,s,t)
				       +0.5*(
					  GAMMA(LatSum_.pDx(),LatSum_.pDX(),
						LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
					  *PSI(LatSum_.pDX(),i,j,m,n)
					  + GAMMA(LatSum_.pDx(),LatSum_.pDX(),
						  LatSum_.Atom(0),LatSum_.Atom(1),
						  i,j,s,t)
					  *PSI(LatSum_.pDX(),k,l,m,n)
					  + GAMMA(LatSum_.pDx(),LatSum_.pDX(),
						  LatSum_.Atom(0),LatSum_.Atom(1),
						  m,n,s,t)
					  *PSI(LatSum_.pDX(),i,j,k,l))));
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
				       OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),i,j)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),k,l)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),s,t))
				    +phi1*(
				       OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),k,l)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,i,j)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),i,j)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,k,l)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),i,j)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),k,l)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,m,n)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),k,l)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				       *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				       *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),m,n)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),s,t)
				       *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
				    +phi2*(
				       OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					     LatSum_.Atom(1),k,l)
				       *XI(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n,s,t)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),i,j)
				       *XI(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n,s,t)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),m,n)
				       *XI(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l,s,t)
				       + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,k,l)
				       + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,i,j)
				       + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),s,t,m,n)));
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
				       PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),s,t))
				    +phi1*(
				       PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),s,t)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),s,t)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),m,n)
				       *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					      LatSum_.Atom(1),s,t)
				       *0.5*PSI(LatSum_.pDX(),i,j,k,l))
				    +phi2*(
				       PI(LatSum_.pDx(),LatSum_.pDX(),k,l)
				       *XI(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t,i,j)
				       + PI(LatSum_.pDx(),LatSum_.pDX(),i,j)
				       *XI(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t,k,l)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),m,n)
				       *THETA(LatSum_.pDX(),LatSum_.Atom(0),
					      LatSum_.Atom(1),i,j,k,l,s,t)
				       + OMEGA(LatSum_.pDx(),LatSum_.Atom(0),
					       LatSum_.Atom(1),s,t)
				       *THETA(LatSum_.pDX(),LatSum_.Atom(0),
					      LatSum_.Atom(1),i,j,k,l,m,n)
				       + GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					       LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
				       + GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					       LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
				       *GAMMA(LatSum_.pDx(),LatSum_.pDX(),
					      LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
				       + SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
				       *0.5*PSI(LatSum_.pDX(),i,j,k,l))
				    +phi3*LAMDA(LatSum_.Atom(0),LatSum_.Atom(1),
						i,j,k,l,m,n,s,t));
			   }
   }

   
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));

   // Pressure terms are zero
   
   return Phi;
}

Matrix MultiLatticeTPP::CondensedModuli()
{
   Matrix stiff = Moduli();
   int intrn = DOFS-6;
   Matrix CM(6,6), IM(intrn,intrn);
   
   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
      {
	 CM[i][j] = stiff[i][j];
      }

   // Make sure there are internal DOF's
   if (intrn)
   {
      for (int i=0;i<intrn;i++)
	 for (int j=0;j<intrn;j++)
	 {
	    IM[i][j] = stiff[6+i][6+j];
	 }
      IM = IM.Inverse();
      
      // Set up Condensed Moduli
      for (int i=0;i<6;i++)
	 for (int j=0;j<6;j++)
	 {
	    for (int m=0;m<intrn;m++)
	       for (int n=0;n<intrn;n++)
	       {
		  CM[i][j] -= stiff[i][6+m]*IM[m][n]*stiff[6+n][j];
	       }
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

CMatrix MultiLatticeTPP::DynamicalStiffness(Vector &Y)
{
   static CMatrix Cy;
   static double pi = 4.0*atan(1.0);
   static complex<double> Ic(0,1);
   int i,j;

   Cy.Resize(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3,0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Cy
      if (LatSum_.Atom(0) != LatSum_.Atom(1))
      {
	 for (i=0;i<DIM3;i++)
	    for (j=0;j<DIM3;j++)
	    {
	       // K != K' terms (i.e., off block (3x3) diagonal terms)
	       Cy[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
		  (-2.0*Del(i,j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::DY)
		   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::D2Y))
		  *exp(-2.0*pi*Ic *
		       (Y[0]*LatSum_.Dx(0) + Y[1]*LatSum_.Dx(1) + Y[2]*LatSum_.Dx(2)));

	       // K==K' components (i.e., Phi(0,k,k) term)
	       Cy[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(0)+j] +=
		  (2.0*Del(i,j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::DY)
		   +4.0*LatSum_.Dx(i)*LatSum_.Dx(j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::D2Y));
	    }
      }
      else
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       Cy[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
		  (-2.0*Del(i,j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::DY)
		   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)
		   *Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
		      NTemp_,LatSum_.r2(),PairPotentials::D2Y))
		  *(exp(-2.0*pi*Ic *
			(Y[0]*LatSum_.Dx(0) + Y[1]*LatSum_.Dx(1)
			 + Y[2]*LatSum_.Dx(2)))
		    - 1.0);
	    }
      }
      // Normalize through the Mass Matrix
      for (i=0;i<DIM3;i++)
	 for (j=0;j<DIM3;j++)
	 {
	    Cy[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j]
	       /= sqrt(AtomicMass_[LatSum_.Atom(0)]*AtomicMass_[LatSum_.Atom(1)]);
	 }
   }
   
   return Cy;
}

int MultiLatticeTPP::BlochWave(Vector &Y)
{
   static CMatrix A(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3),
      U(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3),
      D(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3);

   static Matrix DefGrad(DIM3,DIM3);
   static Vector Z(DIM3);

   // Setup DefGrad
   DefGrad[0][0] = DOF_[0];
   DefGrad[1][1] = DOF_[1];
   DefGrad[2][2] = DOF_[2];
   DefGrad[0][1] = DefGrad[1][0] = DOF_[3];
   DefGrad[0][2] = DefGrad[2][0] = DOF_[4];
   DefGrad[2][1] = DefGrad[1][2] = DOF_[5];
   
   // Iterate over points in cubic unit cell
   for (UCIter_.Reset();!UCIter_.Done();++UCIter_)
   {
      for (int i=0;i<DIM3;++i)
      {
	 Y[i] = UCIter_[i];
      }

      Z = DefGrad*Y;
      A = DynamicalStiffness(Z);

      Cholesky(A,U,D);
      
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 // if w^2 <= 0.0 --> Re(i*w*x) > 0 --> growing solutions --> unstable
	 if ( real(D[i][i]) <= 0.0 )
	 {
	    return 0;
	 }
      }
   }
   return 1;
}

void MultiLatticeTPP::Print(ostream &out,PrintDetail flag)
{
   static int W;
   static int NoNegEigVal;
   static double MinEigVal;
   static double energy;
   static Matrix
      stress(1,DOFS),
      stiff(DOFS,DOFS),
      EigenValues(1,DOFS),
      CondEV(1,6);
   static Matrix
      CondModuli(6,6);
   static int RankOneConvex;
   static Vector Y(DIM3);
   static int BlochWaveStable;
   static Matrix
      RefLatticeInv = RefLattice_.Inverse();
   
   W=out.width();

   out.width(0);
   cout.width(0);

   NoNegEigVal=0;

   energy = Energy();
   stress = Stress();
   stiff = Stiffness();
   
   EigenValues=SymEigVal(stiff);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<DOFS;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   CondModuli = CondensedModuli();

   CondEV=SymEigVal(CondModuli);
   RankOneConvex = FullScanRank1Convex3D(CondModuli,ConvexityDX_);

   Y.Resize(DIM3,0.0);
   if (NoNegEigVal == 0)
   {
      BlochWaveStable = BlochWave(Y);
      Y = RefLatticeInv.Transpose()*Y;
   }
   else
   {
      BlochWaveStable = -1;
   }


   switch (flag)
   {
      case PrintLong:
	 out << "MultiLatticeTPP:" << endl << endl
	     << "Cell Reference Length: " << setw(W) << RefLen_[0]
	     << setw(W) << RefLen_[1] << setw(W) << RefLen_[2] << endl;
	 out << "Lattice Basis_0 : " << setw(W) << LatticeBasis[0] <<endl
	     << "Lattice Basis_1 : " << setw(W) << LatticeBasis[1] <<endl
	     << "Lattice Basis_2 : " << setw(W) << LatticeBasis[2] <<endl << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "Atom_" << i << " Position : " << setw(W) << A_[i] << endl;
	 }
	 out << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "Atomic Mass " << i << "  : "
		<< setw(W) << AtomicMass_[i] << endl;
	 }
	 out << "Potential Parameters : " << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    for (int j=i;j<INTERNAL_ATOMS;j++)
	    {
	       out << "[" << i << "][" << j << "] -- "
		   << setw(W) << Potential_[i][j] << endl;
	    }
	 }
	 out  << "Shear Modulus : " << setw(W) << ShearMod_ << endl;
	 // also send to cout
	 cout << "MultiLatticeTPP:" << endl << endl
	      << "Cell Reference Length: " << setw(W) << RefLen_[0]
	      << setw(W) << RefLen_[1] << setw(W) << RefLen_[2] << endl;
	 cout << "Lattice Basis_0 : " << setw(W) << LatticeBasis[0] <<endl
	      << "Lattice Basis_1 : " << setw(W) << LatticeBasis[1] <<endl
	      << "Lattice Basis_2 : " << setw(W) << LatticeBasis[2] <<endl << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "Atom_" << i << " Position : " << setw(W) << A_[i] << endl;
	 }
	 cout << "Influance Distance   : " << setw(W) << InfluanceDist_ << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "Atomic Mass " << i << "  : "
		 << setw(W) << AtomicMass_[i] << endl;
	 }
	 cout << "Potential Parameters : " << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    for (int j=i;j<INTERNAL_ATOMS;j++)
	    {
	       cout << "[" << i << "][" << j << "] -- "
		    << setw(W) << Potential_[i][j] << endl;
	    }
	 }
	 cout  << "Shear Modulus : " << setw(W) << ShearMod_ << endl;	 
	 // passthrough to short
      case PrintShort:
	 out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	     << "Pressure (G Normalized): " << setw(W) << Pressure_ << endl
	     << "DOF's :" << endl << setw(W) << DOF_ << endl
	     << "Potential Value (G Normalized):" << setw(W) << energy << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "BodyForce Value " << i << " (Inf Normalized):"
		<< setw(W) << BodyForce_[i] << endl;
	 }
	 out << "Stress (G Normalized):" << setw(W) << stress << endl
	     << "Stiffness (G Normalized):" << setw(W) << stiff
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl
	     << "Condensed Moduli (G Normalized):" << setw(W) << CondModuli
	     << "CondEV Info:" << setw(W) << CondEV
	     << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl
	     << "BlochWave Stability:" << setw(W) << BlochWaveStable << ", "
	     << setw(W) << Y << endl;
	 // send to cout also
	 cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	      << "Pressure (G Normalized): " << setw(W) << Pressure_ << endl
	      << "DOF's :" << endl << setw(W) << DOF_ << endl
	      << "Potential Value (G Normalized):" << setw(W) << energy << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "BodyForce Value " << i << " (Inf Normalized):"
		 << setw(W) << BodyForce_[i] << endl;
	 }
	 cout << "Stress (G Normalized):" << setw(W) << stress << endl
	      << "Stiffness (G Normalized):" << setw(W) << stiff
	      << "Eigenvalue Info:"  << setw(W) << EigenValues
	      << "Bifurcation Info:" << setw(W) << MinEigVal
	      << setw(W) << NoNegEigVal << endl
	      << "Condensed Moduli (G Normalized):" << setw(W) << CondModuli
	      << "CondEV Info:" << setw(W) << CondEV
	      << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl
	      << "BlochWave Stability:" << setw(W) << BlochWaveStable << ", "
	      << setw(W) << Y << endl;
	 break;
   }
}

ostream &operator<<(ostream &out,MultiLatticeTPP &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
