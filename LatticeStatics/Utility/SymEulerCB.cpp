#include "SymEulerCB.h"

SymEulerCB::SymEulerCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);
   Reset();
}

void SymEulerCB::Reset()
{
   int i,q,p;
   F_[0][0] = (*DOF_)[0];
   F_[1][1] = (*DOF_)[1];
   F_[2][2] = (*DOF_)[2];
   F_[0][1] = F_[1][0] = (*DOF_)[3];
   F_[0][2] = F_[2][0] = (*DOF_)[4];
   F_[1][2] = F_[2][1] = (*DOF_)[5];
   
   S_[0][0] = 0.0;
   S_[0][1] = 0.0;
   S_[0][2] = 0.0;
   i=6;
   for (q=1;q<InternalAtoms_;++q)
   {
      for (p=0;p<DIM3;p++)
      {
         S_[q][p] = (*DOF_)[i++];
      }
   }
}

double SymEulerCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += X[j]*(*RefLattice_)[j][i];
   }

   return tmp;
}

double SymEulerCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += F_[i][j]*DX(X,p,q,j);
   }
   tmp += InternalPOS_[q][i] - InternalPOS_[p][i] + S_[q][i] - S_[p][i];

   return tmp;
}

double SymEulerCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + Dx[s]*DX[r]);
}

double SymEulerCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 0.5*(Del(r,t)*DX[s]*DX[u] + Del(r,u)*DX[s]*DX[t]
	       + Del(s,t)*DX[r]*DX[u] + Del(s,u)*DX[r]*DX[t]);
}

double SymEulerCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   return 2.0*DELTA(i,p,q)*Dx[j];
}

double SymEulerCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   return 2.0*DELTA(i,p,q)*DELTA(k,p,q)*Del(j,l);
}

double SymEulerCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   return DELTA(k,p,q)*(Del(i,l)*DX[j] + Del(j,l)*DX[i]);
}

double SymEulerCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 0.0;
}

double SymEulerCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   return 0.0;
}

double SymEulerCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return 0.0;
}
