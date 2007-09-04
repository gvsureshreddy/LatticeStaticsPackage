#include "SymLagrangeCB.h"

SymLagrangeCB::SymLagrangeCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms,DIM3);
   Reset();
}

void SymLagrangeCB::Reset()
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

double SymLagrangeCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += (X[j] + ((InternalPOS_[q][j] + S_[q][j])
		      - (InternalPOS_[p][j] + S_[p][j])))
	 *(*RefLattice_)[j][i];
   }

   return tmp;
}

double SymLagrangeCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int k=0;k<DIM3;++k)
   {
      tmp += F_[i][k]*DX(X,p,q,k);
   }
   
   return tmp;
}

double SymLagrangeCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + Dx[s]*DX[r]);
}

double SymLagrangeCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 0.5*(Del(r,t)*DX[s]*DX[u] +
	       Del(r,u)*DX[s]*DX[t] +
	       Del(s,t)*DX[r]*DX[u] +
	       Del(s,u)*DX[r]*DX[t]);
}

double SymLagrangeCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   double ret=0;
   
   ret=0;
   if (DELTA(i,p,q))
   {
      for (int r=0;r<DIM3;r++)
      {
	 for (int k=0;k<DIM3;k++)
	 {
	    ret += F_[k][r]*(*RefLattice_)[j][r]*Dx[k];
	 }
      }
      ret *= 2.0*DELTA(i,p,q);
   }

   return ret;
}

double SymLagrangeCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   if (DELTA(i,p,q)*DELTA(k,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 for (int t=0;t<DIM3;t++)
	 {
	    for (int r=0;r<DIM3;r++)
	    {
	       tmp += F_[t][r]*(*RefLattice_)[j][r]*F_[t][s]*(*RefLattice_)[l][s];
	    }
	 }
      }
      tmp *= 2.0*DELTA(i,p,q)*DELTA(k,p,q);
   }
   
   return tmp;
}

double SymLagrangeCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   double tmp=0;

   if (DELTA(k,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 tmp += F_[i][s]*(*RefLattice_)[l][s]*DX[j] + F_[j][s]*(*RefLattice_)[l][s]*DX[i];
      }
      tmp = DELTA(k,p,q)*(tmp + Dx[i]*(*RefLattice_)[l][j] + Dx[j]*(*RefLattice_)[l][i]);
   }

   return tmp;
}

double SymLagrangeCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return (0.5*DELTA(m,p,q)*(Del(i,k)*(*RefLattice_)[n][l]*DX[j]
			     + Del(i,k)*DX[l]*(*RefLattice_)[n][j]
			     + Del(i,l)*(*RefLattice_)[n][k]*DX[j]
			     + Del(i,l)*DX[k]*(*RefLattice_)[n][j]
			     + Del(j,k)*(*RefLattice_)[n][l]*DX[i]
			     + Del(j,k)*DX[l]*(*RefLattice_)[n][i]
			     + Del(j,l)*(*RefLattice_)[n][k]*DX[i]
			     + Del(j,l)*DX[k]*(*RefLattice_)[n][i]));
}

double SymLagrangeCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   double tmp=0;

   if (DELTA(i,p,q)*DELTA(k,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 tmp += ((*RefLattice_)[j][n]*F_[m][s]*(*RefLattice_)[l][s]
		 + (*RefLattice_)[j][m]*F_[n][s]*(*RefLattice_)[l][s]
		 + F_[m][s]*(*RefLattice_)[j][s]*(*RefLattice_)[l][n]
		 + F_[n][s]*(*RefLattice_)[j][s]*(*RefLattice_)[l][m]);
      }
      tmp *= DELTA(i,p,q)*DELTA(k,p,q);
   }

   return tmp;
}

double SymLagrangeCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return (0.5*DELTA(m,p,q)*DELTA(a,p,q)*
	   (Del(i,k)*(*RefLattice_)[n][l]*(*RefLattice_)[b][j]
	    + Del(i,k)*(*RefLattice_)[b][l]*(*RefLattice_)[n][j]
	    + Del(i,l)*(*RefLattice_)[n][k]*(*RefLattice_)[b][j]
	    + Del(i,l)*(*RefLattice_)[b][k]*(*RefLattice_)[n][j]
	    + Del(j,k)*(*RefLattice_)[n][l]*(*RefLattice_)[b][i]
	    + Del(j,k)*(*RefLattice_)[b][l]*(*RefLattice_)[n][i]
	    + Del(j,l)*(*RefLattice_)[n][k]*(*RefLattice_)[b][i]
	    + Del(j,l)*(*RefLattice_)[b][k]*(*RefLattice_)[n][i]));
}
