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
   delete [] AtomPositions_;
   delete [] MovableAtoms_;
}

MultiLatticeTPP::MultiLatticeTPP(char *datafile,const char *prefix,int Echo)
{
   Echo_ = Echo;
   // Get Lattice definition
   char tmp[LINELENGTH];
   if(!GetParameter(prefix,"InternalAtoms",datafile,"%u",&INTERNAL_ATOMS)) exit(-1);
   DOFS = 6 + 3*(INTERNAL_ATOMS-1);
   for (int i=0;i<DIM3;++i)
   {
      LatticeBasis[i].Resize(DIM3);
      sprintf(tmp,"LatticeBasis_%u",i);
      if(!GetVectorParameter(prefix,tmp,datafile,&(LatticeBasis[i]))) exit(-1);
   }
   
   // First Size DOF
   DOF_.Resize(DOFS,0.0);
   DOF_[0]=DOF_[1]=DOF_[2] = 1.0;
   // Set RefLattice_
   RefLattice_.Resize(DIM3,DIM3);
   // Set AtomPositions_
   AtomPositions_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      AtomPositions_[i].Resize(DIM3);
      sprintf(tmp,"AtomPosition_%u",i);
      if(!GetVectorParameter(prefix,tmp,datafile,&(AtomPositions_[i]))) exit(-1);
   }
   // Set MovableAtoms_ : atoms which can move in initialization.
   if(!GetParameter(prefix,"NoMovableAtoms",datafile,"%u",&NoMovable_)) exit(-1);
   MovableAtoms_ = new int[NoMovable_];
   if((NoMovable_ > 0) &&
      (!GetIntVectorParameter(prefix,"MovableAtoms",datafile,NoMovable_,MovableAtoms_)))
      exit(-1);
   
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
	 Potential_[i][j] = Potential_[j][i]
	    = InitializePairPotential(datafile,prefix,i,j);
      }

      sprintf(tmp,"AtomicMass_%u",i);
      if(!GetParameter(prefix,tmp,datafile,"%lf",&(AtomicMass_[i]))) exit(-1);
   }

	 
   // Get Lattice parameters
   
   if(!GetParameter(prefix,"RefLen_0",datafile,"%lf",&(RefLen_[0]))) exit(-1);
   if(!GetParameter(prefix,"RefLen_1",datafile,"%lf",&(RefLen_[1]))) exit(-1);
   if(!GetParameter(prefix,"RefLen_2",datafile,"%lf",&(RefLen_[2]))) exit(-1);
   
   if(!GetParameter(prefix,"InfluanceDist",datafile,"%u",&InfluanceDist_)) exit(-1);
   if(!GetParameter(prefix,"NTemp",datafile,"%lf",&NTemp_)) exit(-1);
   if(!GetParameter(prefix,"Pressure",datafile,"%lf",&Pressure_)) exit(-1);
   if(!GetParameter(prefix,"ConvexityDX",datafile,"%lf",&ConvexityDX_)) exit(-1);
   
   // needed to initialize reference length
   int iter;
   double DX;
   if(!GetParameter(prefix,"MaxIterations",datafile,"%u",&iter)) exit(-1);
   if(!GetParameter(prefix,"InitializeStepSize",datafile,"%lf",&DX)) exit(-1);
   if(!GetParameter(prefix,"BlochWaveGridSize",datafile,"%u",&GridSize_)) exit(-1);
   
   // Initialize RefLattice_
   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
	 RefLattice_[i][j] = LatticeBasis[i][j]*RefLen_[i];

   // Initiate the Lattice Sum object
   LatSum_(&DOF_,&RefLattice_,INTERNAL_ATOMS,AtomPositions_,Potential_,&InfluanceDist_,
	   &NTemp_);

   int err=0;
   err=FindLatticeSpacing(iter,DX);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << endl;
      exit(-1);
   }

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   UCIter_(GridSize_);
   if (EnterDebugMode())
   {
      cout << setw(15);
      DebugMode();
   }
   
}

// This routing assumes the lattice has fixed angles
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

   Vector RHS(DIM3+NoMovable_,0.0);
   Matrix Stiff(DIM3+NoMovable_,DIM3+NoMovable_,0.0);
   Vector De(DIM3+NoMovable_,0.0);
   Matrix stress=Stress();
   Matrix stiff=stiffness();
   Vector LatBasisNorm[3];

   for (int i=0;i<DIM3;++i)
   {
      LatBasisNorm[i].Resize(DIM3);
      LatBasisNorm[i] = LatticeBasis[i]/LatticeBasis[i].Norm();
   }
   
   for (int i=0;i<DIM3;++i)
   {
      for (int j=0;j<DIM3;++j)
	 // Don't forget the 2's in the stress[][]
	 for (int k=j;k<DIM3;k++)
	 {
	    RHS[i] += LatBasisNorm[i][j]*stress[0][INDU(j,k)]*LatBasisNorm[i][k];
	 }
   }
   for (int i=0;i<NoMovable_;++i)
   {
      RHS[DIM3+i] = stress[0][MovableAtoms_[i]];
   }

   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
      {
	 Stiff[i][j] = 0.0;
	 for (int k=0;k<DIM3;++k)
	    for (int l=k;l<DIM3;++l)
	       for (int m=0;m<DIM3;++m)
		  for (int n=m;n<DIM3;++n)
		  {
		     Stiff[i][j] +=
			LatBasisNorm[i][k]*LatBasisNorm[i][l]
			*stiff[INDU(k,l)][INDU(m,n)]
			*LatBasisNorm[j][m]*LatBasisNorm[j][n];
		  }
      }
   for (int i=0;i<NoMovable_;++i)
      for (int j=0;j<DIM3;++j)
      {
	 Stiff[DIM3+i][j] = 0.0;
	 for (int m=0;m<DIM3;++m)
	    for (int n=m;n<DIM3;++n)
	    {
	       Stiff[DIM3+i][j] += stiff[MovableAtoms_[i]][INDU(m,n)]
		  *LatBasisNorm[j][m]*LatBasisNorm[j][n];
	    }
	 Stiff[j][DIM3+i] = Stiff[DIM3+i][j];
      }
   for (int i=0;i<NoMovable_;++i)
      for (int j=0;j<NoMovable_;++j)
      {
	 Stiff[DIM3+i][DIM3+j] = stiff[MovableAtoms_[i]][MovableAtoms_[j]];
      }

   int itr=0;
   do
   {
      itr++;

#ifdef SOLVE_SVD
      De = SolveSVD(Stiff,RHS,MAXCONDITION,Echo_);
#else
      De = SolvePLU(Stiff,RHS);
#endif

      for (int i=0;i<DIM3;++i)
	 RefLen_[i] -= De[i];
      for (int i=0;i<NoMovable_;++i)
      {
	 AtomPositions_[1+(MovableAtoms_[i]-6)/3][(MovableAtoms_[i]-6)%3] -= De[DIM3+i];
      }
   
      // Update RefLattice_
      for (int p=0;p<DIM3;++p)
	 for (int q=0;q<DIM3;++q)
	    RefLattice_[p][q] = LatticeBasis[p][q]*RefLen_[p];

      LatSum_.Recalc();
      stress = Stress();
      stiff = Stiffness();
      
      for (int i=0;i<DIM3;++i)
      {
	 RHS[i] = 0.0;
	 for (int j=0;j<DIM3;++j)
	    // Don't forget the 2's in stress[][]
	    for (int k=j;k<DIM3;k++)
	    {
	       RHS[i] += LatBasisNorm[i][j]*stress[0][INDU(j,k)]*LatBasisNorm[i][k];
	    }
      }
      for (int i=0;i<NoMovable_;++i)
      {
	 RHS[DIM3+i] = stress[0][MovableAtoms_[i]];
      }
      
      for (int i=0;i<DIM3;++i)
	 for (int j=0;j<DIM3;++j)
	 {
	    Stiff[i][j] = 0.0;
	    for (int k=0;k<DIM3;++k)
	       for (int l=k;l<DIM3;++l)
		  for (int m=0;m<DIM3;++m)
		     for (int n=m;n<DIM3;++n)
		     {
			Stiff[i][j] +=
			   LatBasisNorm[i][k]*LatBasisNorm[i][l]
			   *stiff[INDU(k,l)][INDU(m,n)]
			   *LatBasisNorm[j][m]*LatBasisNorm[j][n];
		     }
	 }
      for (int i=0;i<NoMovable_;++i)
	 for (int j=0;j<DIM3;++j)
	 {
	    Stiff[DIM3+i][j] = 0.0;
	    for (int m=0;m<DIM3;++m)
	       for (int n=m;n<DIM3;++n)
	       {
		  Stiff[DIM3+i][j] += stiff[MovableAtoms_[i]][INDU(m,n)]
		     *LatBasisNorm[j][m]*LatBasisNorm[j][n];
	       }
	    Stiff[j][DIM3+i] = Stiff[DIM3+i][j];
	 }
      for (int i=0;i<NoMovable_;++i)
	 for (int j=0;j<NoMovable_;++j)
	 {
	    Stiff[DIM3+i][DIM3+j] = stiff[MovableAtoms_[i]][MovableAtoms_[j]];
	 }

      if (Echo_)
      {
	 cout << setw(10) << itr << endl
	      << "    RHS=" << setw(20) << RHS << setw(20) << RHS.Norm() << endl
	      << "     De=" << setw(20) << De  << setw(20) << De.Norm() << endl
	      << "RefLen_=" << setw(20) << RefLen_[0] << setw(20) << RefLen_[1]
	      << setw(20) << RefLen_[2] << endl;
	 if (NoMovable_) cout << "AtomPositions_=";
	 for (int i=0;i<NoMovable_;++i)
	 {
	    cout << setw(20)
		 << AtomPositions_[1+(MovableAtoms_[i]-6)/3][(MovableAtoms_[i]-6)%3];
	 }
	 if (NoMovable_) cout << endl;
      }
   }
   while ((itr < iter)
	  && ((RHS.Norm() > 1.0e-13) || (De.Norm() > 1.0e-14)));

   ShearMod_ = 0.25*fabs(Stiffness()[5][5]);
   NTemp_=oldTemp;
   Pressure_=oldPressure;
   //Normalize the Pressure
   // Should be input in normalized form

   LatSum_.Recalc();
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
      phi = 2.0*sqrt(LatSum_.r2())*LatSum_.phi1();
      if (ForceNorm < fabs(-phi/2.0))
      {
	 ForceNorm = fabs(-phi/2.0);
      }
      for (i=0;i<DIM3;i++)
      {
	 BodyForce_[LatSum_.Atom(0)][i] += -phi*LatSum_.Dx(i)/(2.0*sqrt(LatSum_.r2()));
      }

      // Claculate the Stress
      if (dt == PairPotentials::T0)
	 phi=LatSum_.phi1();
      else
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
      if (dt==PairPotentials::T0)
      {
	 phi = LatSum_.phi2();
	 phi1 = LatSum_.phi1();
      }
      else
      {
	 phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	    NTemp_,LatSum_.r2(),PairPotentials::D2Y,dt);
	 phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
	    NTemp_,LatSum_.r2(),PairPotentials::DY,dt);
      }
      
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
      //Lower Diag Block (3*INTERNAL_ATOMS,3*INTERNAL_ATOMS)
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
      phi1=LatSum_.phi2();
      phi2=LatSum_.phi1();
	 
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
      phi2=LatSum_.phi2();
      phi3=LatSum_.phi1();
      
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
   Matrix stiff = Stiffness();
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

int MultiLatticeTPP::comp(const void *a,const void *b)
{
   double t;
   if( *((double*) a) == *((double*) b)) return 0;
   else
   {
      t= *((double*) a) - *((double*) b);
      t/=fabs(t);
      return int(t);
   }
}

int MultiLatticeTPP::abscomp(const void *a,const void *b)
{
   double t;
   if( fabs(*((double*) a)) == fabs(*((double*) b))) return 0;
   else
   {
      t= fabs(*((double*) a)) - fabs(*((double*) b));
      t/=fabs(t);
      return int(t);
   }
}

void MultiLatticeTPP::interpolate(Matrix *EigVals,int zero,int one,int two)
{
   // Calculate expected value for eigvals and store in zero position
   EigVals[zero] = 2.0*EigVals[one] - EigVals[zero];

   double delta,dtmp;
   int i,j,pos;

   for (i=0;i<EigVals[0].Cols();++i)
   {
      pos = i;
      delta = fabs(EigVals[zero][0][i] - EigVals[two][0][i]);
      for (j=i+1;j<EigVals[0].Cols();++j)
      {
	 dtmp = fabs(EigVals[zero][0][i] - EigVals[two][0][j]);
	 if (dtmp < delta)
	 {
	    delta = dtmp;
	    pos = j;
	 }
      }
      // move correct eigval to current pos
      dtmp = EigVals[two][0][i];
      EigVals[two][0][i] = EigVals[two][0][pos];
      EigVals[two][0][pos] = dtmp;
   }
}

CMatrix MultiLatticeTPP::ReferenceDynamicalStiffness(Vector &K)
{
   static CMatrix Dk;
   static double pi = 4.0*atan(1.0);
   static complex<double> Ic(0,1);
   static complex<double> A = 2.0*pi*Ic;
   int i,j;

   Dk.Resize(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3,0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Dk
      if (LatSum_.Atom(0) != LatSum_.Atom(1))
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       // y != y' terms (i.e., off block (3x3) diagonal terms)
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
		  (-2.0*Del(i,j)*LatSum_.phi1()
		   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2())
		  *exp(A*
		       (K[0]*LatSum_.DX(0) + K[1]*LatSum_.DX(1) + K[2]*LatSum_.DX(2)));

	       // y==y' components (i.e., Phi(0,y,y) term)
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(0)+j] +=
		  (2.0*Del(i,j)*LatSum_.phi1()
		   +4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2());
	    }
      }
      else
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
		  (-2.0*Del(i,j)*LatSum_.phi1()
		   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2())
		  *(exp(A*
			(K[0]*LatSum_.DX(0) + K[1]*LatSum_.DX(1)
			 + K[2]*LatSum_.DX(2)))
		    - 1.0);
	    }
      }
   }
   // Normalize through the Mass Matrix
   for (int p=0;p<INTERNAL_ATOMS;++p)
      for (int q=0;q<INTERNAL_ATOMS;++q)
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       Dk[DIM3*p+i][DIM3*q+j] /= sqrt(AtomicMass_[p]*AtomicMass_[q]);
	    }
   
   return Dk;
}

void MultiLatticeTPP::ReferenceDispersionCurves(Vector K,int NoPTS,const char *prefix,
						ostream &out)
{
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);

   Matrix Tmp(DIM3,DIM3),InverseLat(DIM3,DIM3);
   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
      {
	 Tmp[i][j] = RefLattice_[i][j];
      }
   InverseLat = Tmp.Inverse();

   Matrix EigVal[3];
   for (int i=0;i<3;++i) EigVal[i].Resize(1,INTERNAL_ATOMS*DIM3);

   Vector Z1(DIM3),Z2(DIM3);
   for (int k=0;k<DIM3;++k)
   {
      Z1[k] = K[k];
      Z2[k] = K[DIM3 + k];
   }
   Z1 = InverseLat*Z1;
   Z2 = InverseLat*Z2;
   
   Vector Z(DIM3),
      DZ=Z2-Z1;
   double dz = 1.0/(NoPTS-1);
   for (int k=0;k<2;++k)
   {
      Z = Z1 + (k*dz)*DZ;
      EigVal[k] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[k][0],INTERNAL_ATOMS*DIM3,sizeof(double),&comp);
      
      out << prefix << setw(w) << NTemp_ << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << NTemp_ << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 out << setw(w) << EigVal[k][0][i];
	 if (Echo_) cout << setw(w) << EigVal[k][0][i];
      }
      out << endl;
      if (Echo_) cout << endl;
   }
   int zero=0,one=1,two=2;
   for (int k=2;k<NoPTS;++k)
   {
      Z = Z1 + (k*dz)*DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[two][0],INTERNAL_ATOMS*DIM3,sizeof(double),&comp);
      interpolate(EigVal,zero,one,two);
      
      out << prefix << setw(w) << NTemp_ << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << NTemp_ << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 out << setw(w) << EigVal[two][0][i];;
	 if (Echo_) cout << setw(w) << EigVal[two][0][i];;
      }
      out << endl;
      if (Echo_) cout << endl;

      zero = (++zero)%3; one = (zero+1)%3; two = (one+1)%3;
   }
}

int MultiLatticeTPP::ReferenceBlochWave(Vector &K)
{
   static CMatrix A(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3);
   static Matrix EigVals(1,INTERNAL_ATOMS*DIM3);
   static Matrix InverseLat(DIM3,DIM3),Tmp(DIM3,DIM3);
   static Vector Z(DIM3);

   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
      {
	 Tmp[i][j] = RefLattice_[i][j];
      }
   InverseLat = Tmp.Inverse();

   // Iterate over points in cubic unit cell
   for (UCIter_.Reset();!UCIter_.Done();++UCIter_)
   {
      for (int i=0;i<DIM3;++i)
      {
	 K[i] = UCIter_[i];
      }

      Z = InverseLat*K;
      A = ReferenceDynamicalStiffness(Z);

      EigVals = HermiteEigVal(A);
      
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 // if w^2 <= 0.0 --> Re(i*w*x) > 0 --> growing solutions --> unstable
	 if ( EigVals[0][i] <= 0.0 )
	 {
	    return 0;
	 }
      }
   }
   return 1;
}

void MultiLatticeTPP::LongWavelengthModuli(double dk, int gridsize,const char *prefix,
					   ostream &out)
{
   static double pi = 4*atan(1.0);
   static double twopi = 2*pi;
   double GS = double(gridsize);
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);

   Matrix
      Lp=CondensedModuli(),
      Ap(DIM3,DIM3),
      A(DIM3,DIM3);

   //----------------  setup L condensed moduli wrt F NOT U ----------------------
   Matrix Phi(9,9,0.0),
      Dpp((INTERNAL_ATOMS-1)*3,(INTERNAL_ATOMS-1)*3,0.0),
      Dfp((INTERNAL_ATOMS-1)*3,9,0.0);
   double phi,phi1,tmp[3][3][3];
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi = LatSum_.phi2();
      phi1 = LatSum_.phi1();

      // upper 9x9 block
      for (int i=0;i<DIM3;i++)
      {
	 for (int j=0;j<DIM3;j++)
	 {
	    for (int k=0;k<DIM3;k++)
	    {
	       for (int l=0;l<DIM3;l++)
	       {
		  Phi[3*i+j][3*k+l]+=
		     4.0*phi*(LatSum_.Dx(i)*LatSum_.DX(j))
			  *(LatSum_.Dx(k)*LatSum_.DX(l))
		     +2*phi1*(Del(i,k)*LatSum_.DX(j)*LatSum_.DX(l));
	       }
	    }
	 }
      }

      // lower block
      for (int i=1;i<INTERNAL_ATOMS;++i)
      {
	 for (int j=0;j<DIM3;++j)
	 {
	    for (int k=1;k<INTERNAL_ATOMS;++k)
	    {
	       for (int l=0;l<DIM3;++l)
	       {
		  Dpp[3*(i-1)+j][3*(k-1)+l]+=
		     phi*OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
		     *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l)
		     +phi1*SIGMA(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
	       }
	    }
	 }
      }
      
      // off-diagonal block
      for (int i=0;i<DIM3;++i)
	 for (int j=0;j<DIM3;++j)
	    for (int l=0;l<DIM3;++l)
	    {
	       tmp[i][j][l] = 0.0;
	       for (int k=0;k<DIM3;++k)
	       {
		  tmp[i][j][l] += LatticeBasis[l][k]*DOF_[INDU(k,i)]*LatSum_.DX(j)
		     +LatticeBasis[l][k]*DOF_[INDU(i,k)]*LatSum_.DX(j);
	       }
	    }
      for (int k=1;k<INTERNAL_ATOMS;++k)
	 for (int l=0;l<DIM3;++l)
	    for (int i=0;i<DIM3;++i)
	       for (int j=0;j<DIM3;++j)
	       {
		  Dfp[3*(k-1)+l][3*i+j] +=
		     phi*(2.0*LatSum_.Dx(i)*LatSum_.DX(j))
		     *OMEGA(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l) +
		     phi1*0.5*DELTA(k,LatSum_.Atom(0),LatSum_.Atom(1))
		     *( LatticeBasis[l][i]*LatSum_.Dx(j) +
			LatticeBasis[l][j]*LatSum_.Dx(i) +
			tmp[i][j][l]);
	       }
   }

   // Condense moduli.
   Phi -= (Dfp.Transpose())*(Dpp.Inverse())*Dfp;
   // Phi = Phi/(2*Vr*ShearMod)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*ShearMod_));
   
   //-----------------------------------------------------------------------------

   Matrix SymPhi(6,6,0.0);

   for (int i=0;i<DIM3;++i)
      for (int j=i;j<DIM3;++j)
	 for (int k=0;k<DIM3;++k)
	    for (int l=k;l<DIM3;++l)
	    {
	       SymPhi[INDU(i,j)][INDU(k,l)] = 0.25*(
		  Phi[3*i+j][3*k+l] +
		  Phi[3*j+i][3*k+l] +
		  Phi[3*i+j][3*l+k] +
		  Phi[3*j+i][3*l+k]);
	    }

   cout << "Condensed Moduli Check!:" << endl;
   cout << setw(w) << Lp-SymPhi << endl;

   

   Vector K(DIM3),Z(DIM3,0.0);
   Matrix BlkEigVal(1,INTERNAL_ATOMS*DIM3);
   Matrix ModEigVal(1,DIM3);

   double Mc=0.0;
   double Vc=RefLattice_.Det();
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      Mc += AtomicMass_[i];
   }

   for (int phi=0;phi<gridsize;++phi)
   {
      for (int theta=0;theta<gridsize;++theta)
      {
	 K[0] = sin(pi*(phi/GS))*cos(twopi*(theta/GS));
	 K[1] = sin(pi*(phi/GS))*sin(twopi*(theta/GS));
	 K[2] = cos(pi*(phi/GS));

	 Z=dk*K;
	 BlkEigVal = HermiteEigVal(ReferenceDynamicalStiffness(Z));

	 // sort by absolute value
	 qsort(BlkEigVal[0],INTERNAL_ATOMS*DIM3,sizeof(double),&abscomp);

	 for (int i=0;i<DIM3;++i)
	 {
	    // wave speed squared
	    BlkEigVal[0][i] /= (twopi*dk*twopi*dk);
	 }

	 for (int i=0;i<DIM3;++i)
	    for (int j=0;j<DIM3;++j)
	    {
	       A[i][j] = 0.0;
	       for (int k=0;k<DIM3;++k)
		  for (int l=0;l<DIM3;++l)
		  {
		     A[i][j] += Phi[3*i+k][3*j+l]*K[k]*K[l];
		  }
	    }
	 
	 ModEigVal = SymEigVal(A);
	 qsort(ModEigVal[0],DIM3,sizeof(double),&abscomp);
	 for (int i=0;i<3;++i)
	 {
	    // normalize by G/(Mc/Vc)
	    ModEigVal[0][i] *= ShearMod_/(Mc/Vc);
	 }

	 out << prefix << setw(w/2) << phi << setw(w/2) << theta;
	 if (Echo_) cout << prefix << setw(w/2) << phi << setw(w/2) << theta;
	 for (int i=0;i<DIM3;++i)
	 {
	    out << setw(w) << ModEigVal[0][i];
	    if (Echo_) cout << setw(w) << ModEigVal[0][i];
	 }
	 for (int i=0;i<DIM3;++i)
	 {
	    out << setw(w) << BlkEigVal[0][i];
	    if (Echo_) cout << setw(w) << BlkEigVal[0][i];
	 }
	 for (int i=0;i<DIM3;++i)
	 {
	    out << setw(w) << (ModEigVal[0][i]-BlkEigVal[0][i])/ModEigVal[0][i];
	    if (Echo_) cout << setw(w)
			    << (ModEigVal[0][i]-BlkEigVal[0][i])/ModEigVal[0][i];
	 }
	 out << endl;
	 if (Echo_) cout << endl;
      }
      out << endl;
      if (Echo_) cout << endl;
   }
}

void MultiLatticeTPP::NeighborDistances(int cutoff,ostream &out)
{
   Matrix NeighborDist = LatSum_.NeighborDistances(cutoff,pow(10,-(out.precision()-1)));
   
   int W=out.width();
   int types = (INTERNAL_ATOMS*(INTERNAL_ATOMS+1))/2;
   for (int i=0;i<cutoff;++i)
   {
      out << setw(W) << NTemp_ << setw(W) << NeighborDist[i][0];
      for (int j=0;j<types;++j)
      {
	 out << setw(W/4) << int(NeighborDist[i][1+j]);
      }
      out << endl;
   }
   out << endl;
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
   static Vector K(DIM3);
   static int BlochWaveStable;
   
   W=out.width();

   out.width(0);
   if (Echo_) cout.width(0);

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

   K.Resize(DIM3,0.0);
   if (RankOneConvex)
   {
      BlochWaveStable = BlochWave(K);
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
	    out << "Atom_" << i << " Position : "
		<< setw(W) << AtomPositions_[i] << endl;
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
	 if (Echo_)
	 {
	    cout << "MultiLatticeTPP:" << endl << endl
		 << "Cell Reference Length: " << setw(W) << RefLen_[0]
		 << setw(W) << RefLen_[1] << setw(W) << RefLen_[2] << endl;
	    cout << "Lattice Basis_0 : " << setw(W) << LatticeBasis[0] <<endl
		 << "Lattice Basis_1 : " << setw(W) << LatticeBasis[1] <<endl
		 << "Lattice Basis_2 : " << setw(W) << LatticeBasis[2] <<endl << endl;
	    for (int i=0;i<INTERNAL_ATOMS;++i)
	    {
	       cout << "Atom_" << i << " Position : "
		    << setw(W) << AtomPositions_[i] << endl;
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
	 }
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
	     << setw(W) << K << endl;
	 // send to cout also
	 if (Echo_)
	 {
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
		 << "BlochWave Stability (GridSize=" << GridSize_ << "):"
		 << setw(W) << BlochWaveStable << ", "
		 << setw(W) << K << endl;
	 }
	 break;
   }
   // check for debug mode request
   if (EnterDebugMode())
      {
	 cout << setw(15);
	 DebugMode();
      }
}

ostream &operator<<(ostream &out,MultiLatticeTPP &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}


//---------------------- Debug Mode Handler --------------------------


void MultiLatticeTPP::DebugMode()
{
   char *Commands[] = {
      "INTERNAL_ATOMS",                // 0
      "DOFS",                          // 1
      "RefLen_",                       // 2
      "InfluanceDist_",                // 3
      "NTemp_",                        // 4
      "DOF_",                          // 5
      "LatticeBasis",                  // 6
      "RefLattice_",                   // 7
      "ShearMod_",                     // 8
      "Pressure_",                     // 9
      "BodyForce_",                    // 10
      "AtomicMass_",                   // 11
      "GridSize_",                     // 12
      "Potential_",                    // 13
      "ConvexityDX_",                  // 14
      "NoMovable_",                    // 15
      "MovableAtoms_",                 // 16
      "stress",                        // 17
      "stiffness",                     // 18
      "CondensedModuli",               // 19
      "ReferenceDispersionCurves",     // 20
      "ReferenceBlochWave",            // 21
      "ReferenceDynamicalStiffness",   // 22
      "SetDOF",                        // 23
      "StressDT",                      // 24
      "StiffnessDT",                   // 25
      "SetTemp",                       // 26
      "Energy",                        // 27
      "Moduli",                        // 28
      "E3",                            // 29
      "E4",                            // 30
      "SetGridSize",                   // 31
      "NeighborDistances",             // 32
      "Print-short",                   // 33
      "Print-long",                    // 34
      "SetPressure",                   // 35
      "FindLatticeSpacing"             // 36
   };
   int NOcommands=37;
   
   char response[LINELENGTH];
   char prompt[] = "Debug > ";
   int W=cout.width();

   cout << setw(0) << prompt;

   cin.getline(response,LINELENGTH);

   while (strcasecmp(response,"q") &&
	  strcasecmp(response,"quit") &&
	  strcasecmp(response,"exit"))
   {
      if (!strcmp(response,Commands[0]))
	 cout << "INTERNAL_ATOMS = " << INTERNAL_ATOMS << endl;
      else if (!strcmp(response,Commands[1]))
	 cout << "DOFS = " << DOFS << endl;
      else if (!strcmp(response,Commands[2]))
	 cout << "RefLen_[0] = " << RefLen_[0] << endl
	      << "RefLen_[1] = " << RefLen_[1] << endl
	      << "RefLen_[2] = " << RefLen_[2] << endl;
      else if (!strcmp(response,Commands[3]))
	 cout << "InfluanceDist_ = " << InfluanceDist_ << endl;
      else if (!strcmp(response,Commands[4]))
	 cout << "NTemp_ = " << NTemp_ << endl;
      else if (!strcmp(response,Commands[5]))
      {
	 for (int i=0;i<DOFS;++i)
	    cout << "DOF_[" << i << "] = " << DOF_[i] << endl;
      }
      else if (!strcmp(response,Commands[6]))
      {
	 cout << "LatticeBasis[0]= " << setw(W) << LatticeBasis[0] << endl;
	 cout << "LatticeBasis[1]= " << setw(W) << LatticeBasis[1] << endl;
	 cout << "LatticeBasis[2]= " << setw(W) << LatticeBasis[2] << endl;
      }
      else if (!strcmp(response,Commands[7]))
	 cout << "RefLattice_= " << setw(W) << RefLattice_;
      else if (!strcmp(response,Commands[8]))
	 cout << "ShearMod_= " << ShearMod_ << endl;
      else if (!strcmp(response,Commands[9]))
	 cout << "Pressure_= " << Pressure_ << endl;
      else if (!strcmp(response,Commands[10]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "BodyForce_[" << i << "]= " << setw(W)
		 << BodyForce_[i] << endl;
	 }
      }
      else if (!strcmp(response,Commands[11]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "AtomicMass_[" << i << "]= " << setw(W)
		 << AtomicMass_[i] << endl;
	 }
      }
      else if (!strcmp(response,Commands[12]))
	 cout << "GridSize_= " << GridSize_ << endl;
      else if (!strcmp(response,Commands[13]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	    for (int j=i;j<INTERNAL_ATOMS;++j)
	    {
	       cout << "Potential_[" << i << "][" << j << "]= "
		    << setw(W) << Potential_[i][j] << endl;
	    }
      }
      else if (!strcmp(response,Commands[14]))
	 cout << "ConvexityDX_= " << ConvexityDX_ << endl;
      else if (!strcmp(response,Commands[15]))
	 cout << "NoMovable_= " << NoMovable_ << endl;
      else if (!strcmp(response,Commands[16]))
	 for (int i=0;i<NoMovable_;++i)
	    cout << "MovableAtoms_[" << i << "]= " << MovableAtoms_[i] << endl;
      else if (!strcmp(response,Commands[17]))
	 cout << "stress= " << setw(W) << stress();
      else if (!strcmp(response,Commands[18]))
	 cout << "stiffness= " << setw(W) << stiffness();
      else if (!strcmp(response,Commands[19]))
	 cout << "CondensedModuli= " << setw(W) << CondensedModuli();
      else if (!strcmp(response,Commands[20]))
      {
	 Vector K(DIM3,0.0);
	 int NoPTS;
	 char prefix[LINELENGTH];
	 int oldEcho_=Echo_;
	 cout << "\tK > ";
	 cin >> K;
	 cin.sync(); // clear input
	 cout << "\tNoPTS > ";
	 cin >> NoPTS;
	 cin.sync(); // clear input
	 cout << "\tprefix > ";
	 cin >> prefix;
	 cin.sync(); // clear input
	 Echo_=0;
	 cout << "ReferenceDispersionCurves= ";
	 ReferenceDispersionCurves(K,NoPTS,prefix,cout);
	 Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[21]))
      {
	 Vector K(DIM3,0.0);
	 cout << "\tK > ";
	 cin >> K;
	 cin.sync(); // clear input
	 cout << "ReferenceBlochWave= " << ReferenceBlochWave(K);
      }
      else if (!strcmp(response,Commands[22]))
      {
	 cout << "\tK > ";
	 Vector K(3,0.0);
	 cin >> K;
	 cin.sync(); // clear input
	 cout << "ReferenceDynamicalStiffness= "
	      << setw(W) << ReferenceDynamicalStiffness(K) << endl;
      }
      else if (!strcmp(response,Commands[23]))
      {
	 Vector DOF(DOFS,0.0);
	 cout << "\tDOF > ";
	 cin >> DOF;
	 cin.sync(); // clear input
	 SetDOF(DOF);
      }
      else if (!strcmp(response,Commands[24]))
	 cout << "StressDT= " << setw(W) << StressDT();
      else if (!strcmp(response,Commands[25]))
	 cout << "StiffnessDT= " << setw(W) << StiffnessDT();
      else if (!strcmp(response,Commands[26]))
      {
	 double Temp;
	 cout << "\tTemp > ";
	 cin >> Temp;
	 cin.sync(); // clear input
	 SetTemp(Temp);
      }
      else if (!strcmp(response,Commands[27]))
	 cout << "Energy= " << Energy() << endl;
      else if (!strcmp(response,Commands[28]))
	 cout << "Moduli= " << setw(W) << Moduli();
      else if (!strcmp(response,Commands[29]))
	 cout << "E3= " << setw(W) << E3();
      else if (!strcmp(response,Commands[30]))
	 cout << "E4= " << setw(W) << E4();
      else if (!strcmp(response,Commands[31]))
      {
	 int GridSize;
	 cout << "\tGridSize > ";
	 cin >> GridSize;
	 cin.sync(); // clear input
	 SetGridSize(GridSize);
      }
      else if (!strcmp(response,Commands[32]))
      {
	 int oldEcho_=Echo_;
	 int cutoff;
	 cout << "\tcutoff > ";
	 cin >> cutoff;
	 cin.sync(); // clear input
	 Echo_ = 0;
	 NeighborDistances(cutoff,cout);
	 Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[33]))
      {
	 int oldEcho_=Echo_;
	 Echo_=0;
	 cout << setw(W) << *this;
	 Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[34]))
      {
	 int oldEcho_=Echo_;
	 Echo_=0;
	 Print(cout,PrintLong);
	 Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[35]))
      {
	 double pressure;
	 cout << "\tPressure > ";
	 cin >> pressure;
	 cin.sync(); // clear input
	 SetPressure(pressure);
      }
      else if (!strcmp(response,Commands[36]))
      {
	 int iter;
	 double dx;
	 cout << "\titer > ";
	 cin >> iter;
	 cin.sync(); // clear input
	 cout << "\tdx > ";
	 cin >> dx;
	 cin.sync(); // clear input
	 FindLatticeSpacing(iter,dx);
      }
      else if (!strcmp(response,"?") ||
	       !strcasecmp(response,"help"))
      {
	 cout << setiosflags(ios::left);
	 for (int i=0;i<NOcommands/2 + NOcommands%2;++i)
	 {
	    cout << "  " << setw(30) << Commands[i];
	    if (i==NOcommands/2 && !NOcommands%2)
	       cout << endl;
	    else
	       cout << setw(30) << Commands[NOcommands/2+i] << endl;
	    
	    if (!((i+1)%30))
	    {
	       cout << "more...." << endl;
	       char ans;
	       cin.sync(); // clear input
	       ans=kbhitWait();
	       if (ans=='q') break;
	    }
	 }
	 cout << resetiosflags(ios::left) << endl;
      }
      else
      {
	 cout << "!--- Error - Unknown command ---!" << endl << endl;
      }
      
      cout << endl << prompt;
      cin.getline(response,LINELENGTH);
   }  
}
