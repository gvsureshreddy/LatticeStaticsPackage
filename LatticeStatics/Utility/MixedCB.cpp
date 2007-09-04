#include "MixedCB.h"

MixedCB::MixedCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);
   Reset();
}

void MixedCB::Reset()
{
   int i,q,p;
   F_[0][0] = (*DOF_)[0];
   F_[0][1] = (*DOF_)[1];
   F_[0][2] = (*DOF_)[2];
   F_[1][0] = (*DOF_)[3];
   F_[1][1] = (*DOF_)[4];
   F_[1][2] = (*DOF_)[5];
   F_[2][0] = (*DOF_)[6];
   F_[2][1] = (*DOF_)[7];
   F_[2][2] = (*DOF_)[8];
   
   i=9;
   for (q=0;q<InternalAtoms_;++q)
   {
      for (p=0;p<DIM3;p++)
      {
         S_[q][p] = (*DOF_)[i++];
      }
   }
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
