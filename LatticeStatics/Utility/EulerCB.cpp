#include "EulerCB.h"

EulerCB::EulerCB(int InternalAtoms,const char* prefix,const char* datafile):
   CBKinematics(InternalAtoms,prefix,datafile)
{
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);
   
   SetReferenceDOFs();
   Reset();
   
}

void EulerCB::Reset()
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

Vector EulerCB::FractionalPosVec(int p)
{
   Vector pos(DIM3,0.0),fracpos(DIM3,0.0),tmp(DIM3,0.0);
   Matrix CurrentLattice(DIM3,DIM3,0.0),InverseLattice(DIM3,DIM3);
   
   for (int i=0;i<DIM3;++i)
   {
      tmp = CurrentLatticeVec(i);
      for (int j=0;j<DIM3;++j)
      {
         CurrentLattice[i][j] = tmp[j];
      }
      
      pos[i] = InternalPOS_[p][i] + S_[p][i];
   }
   InverseLattice = CurrentLattice.Inverse();
   
   
   for (int i=0;i<DIM3;++i)
   {
      for (int j=0;j<DIM3;++j)
         fracpos[i] += pos[j]*InverseLattice[j][i];
   }
   
   return fracpos;
}

double EulerCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;
   
   for (int j=0;j<DIM3;++j)
   {
      tmp += X[j]*RefLattice_[j][i];
   }
   
   return tmp;
}

double EulerCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;
   
   for (int j=0;j<DIM3;++j)
   {
      tmp += F_[i][j]*DX(X,p,q,j);
   }
   tmp += InternalPOS_[q][i] - InternalPOS_[p][i] + S_[q][i] - S_[p][i];
   
   return tmp;
}

double EulerCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return 2.0*Dx[r]*DX[s];
}

double EulerCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 2.0*Del(r,t)*DX[s]*DX[u];
}

double EulerCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   return 2.0*DELTA(i,p,q)*Dx[j];
}

double EulerCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*Del(j,l);
}

double EulerCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(k,p,q)*Del(i,l)*DX[j];
}

double EulerCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 0.0;
}

double EulerCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   return 0.0;
}

double EulerCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return 0.0;
}
