#include "MixedCB.h"

MixedCB::MixedCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   U_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);
   Reset();
}

double MixedCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += (X[j] + InternalPOS_[q][j] - InternalPOS_[p][j])*(*RefLattice_)[j][i];
   }

   return tmp;
}

double MixedCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += U_[i][j]*DX(X,p,q,j);
   }
   tmp += S_[q][i] - S_[p][i];

   return tmp;
}

double MixedCB::DyDU(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + Dx[s]*DX[r]);
}

double MixedCB::D2yDUU(double *DX,int r, int s, int t, int u)
{
   return 0.5*(DX[u]*DX[s]*Del(t,r) + DX[t]*DX[s]*Del(u,r)
	       + DX[u]*DX[r]*Del(t,s) + DX[t]*DX[r]*Del(u,s));
}

double MixedCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   return 2.0*DELTA(i,p,q)*Dx[j];
}

double MixedCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*Del(j,l);
}

double MixedCB::D2yDUS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   return DELTA(k,p,q)*(Del(i,l)*DX[j] + Del(j,l)*DX[i]);
}

double MixedCB::D3yDUUS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 0.0;
}

double MixedCB::D3yDUSS(int p,int q,int i,int j,int k,int l,int m,int n)
{
   return 0.0;
}

double MixedCB::D4yDUUSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return 0.0;
}
