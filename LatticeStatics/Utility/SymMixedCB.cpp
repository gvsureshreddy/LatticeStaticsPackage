#include "SymMixedCB.h"

SymMixedCB::SymMixedCB(int InternalAtoms,const char* prefix,const char* datafile):
   CBKinematics(InternalAtoms,prefix,datafile)
{
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);

   SetReferenceDOFs();
   Reset();
}

void SymMixedCB::Reset()
{
   int i,j,q,p;
   for (i=0;i<DIM3;++i)
   {
      for (j=0;j<DIM3;++j)
      {
	 F_[i][j] = DOF_[INDF(i,j)];
      }
   }
   
   S_[0][0] = 0.0;
   S_[0][1] = 0.0;
   S_[0][2] = 0.0;
   i=Fsize();
   for (q=1;q<InternalAtoms_;++q)
   {
      for (p=0;p<DIM3;p++)
      {
         S_[q][p] = DOF_[i++];
      }
   }
}

Vector SymMixedCB::FractionalPosVec(int p)
{
   Vector fracpos(DIM3,0.0),tmp(DIM3,0.0);
   Matrix CurrentLattice(DIM3,DIM3,0.0),InverseLattice(DIM3,DIM3);

   for (int i=0;i<DIM3;++i)
   {
      tmp = CurrentLatticeVec(i);
      for (int j=0;j<DIM3;++j)
      {
	 CurrentLattice[i][j] = tmp[j];
      }
   }
   InverseLattice = CurrentLattice.Inverse();
   

   for (int i=0;i<DIM3;++i)
   {
      fracpos[i] += InternalPOS_[p][i];;
      for (int j=0;j<DIM3;++j)
	 fracpos[i] += S_[p][j]*InverseLattice[j][i];
   }

   return fracpos;;
}

double SymMixedCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += (X[j] + InternalPOS_[q][j] - InternalPOS_[p][j])*RefLattice_[j][i];
   }

   return tmp;
}

double SymMixedCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += F_[i][j]*DX(X,p,q,j);
   }
   tmp += S_[q][i] - S_[p][i];

   return tmp;
}

double SymMixedCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + Dx[s]*DX[r]);
}

double SymMixedCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 0.5*(Del(r,t)*DX[s]*DX[u] + Del(r,t)*DX[s]*DX[t]
	       + Del(s,t)*DX[r]*DX[u] + Del(s,u)*DX[r]*DX[t]);
}

double SymMixedCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   return 2.0*DELTA(i,p,q)*Dx[j];
}

double SymMixedCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*Del(j,l);
}

double SymMixedCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   return DELTA(k,p,q)*(Del(i,l)*DX[j] + Del(j,l)*DX[i]);
}

double SymMixedCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 0.0;
}

double SymMixedCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   return 0.0;
}

double SymMixedCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return 0.0;
}
