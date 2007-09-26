#include "MixedCB.h"

MixedCB::MixedCB(int InternalAtoms,const char* prefix,const char* datafile):
   CBKinematics(InternalAtoms,prefix,datafile)
{
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);

   SetReferenceDOFs();
   Reset();
}

void MixedCB::Reset()
{
   int i,j,q,p;
   for (i=0;i<DIM3;++i)
   {
      for (j=0;j<DIM3;++j)
      {
	 F_[i][j] = DOF_[INDF(i,j)];
      }
   }
   
   i=Fsize();
   for (q=0;q<InternalAtoms_;++q)
   {
      for (p=0;p<DIM3;p++)
      {
         S_[q][p] = DOF_[i++];
      }
   }
}

Vector MixedCB::FractionalPosVec(int p)
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

double MixedCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += (X[j] + InternalPOS_[q][j] - InternalPOS_[p][j])*RefLattice_[j][i];
   }

   return tmp;
}

double MixedCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += F_[i][j]*DX(X,p,q,j);
   }
   tmp += S_[q][i] - S_[p][i];

   return tmp;
}

double MixedCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return 2.0*Dx[r]*DX[s];
}

double MixedCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 2.0*Del(r,t)*DX[s]*DX[u];
}

double MixedCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   return 2.0*DELTA(i,p,q)*Dx[j];
}

double MixedCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*Del(j,l);
}

double MixedCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(k,p,q)*Del(i,l)*DX[j];
}

double MixedCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 0.0;
}

double MixedCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   return 0.0;
}

double MixedCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return 0.0;
}
