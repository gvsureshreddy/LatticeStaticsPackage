#include "LagrangeCB.h"

LagrangeCB::LagrangeCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   U_.Resize(DIM3,DIM3);
   V_.Resize(InternalAtoms,DIM3);
   Reset();
}

double LagrangeCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int j=0;j<DIM3;++j)
   {
      tmp += (X[j] + ((InternalPOS_[q][j] + V_[q][j])
		      - (InternalPOS_[p][j] + V_[p][j])))
	 *(*RefLattice_)[j][i];
   }

   return tmp;
}

double LagrangeCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;

   for (int k=0;k<DIM3;++k)
      for (int j=0;j<DIM3;++j)
      {
	 tmp += U_[i][k]*(X[j] + ((InternalPOS_[q][j] + V_[q][j])
				  - (InternalPOS_[p][j] + V_[p][j])))
	    *(*RefLattice_)[j][k];
      }
   
   return tmp;
}

double LagrangeCB::DyDU(double *Dx,double *DX,int r, int s)
{
   return (Dx[r]*DX[s] + DX[r]*Dx[s]);
}

double LagrangeCB::D2yDUU(double *DX,int r, int s, int t, int u)
{
   return 0.5*(Del(r,t)*DX[s]*DX[u] +
	       Del(r,u)*DX[s]*DX[t] +
	       Del(s,t)*DX[r]*DX[u] +
	       Del(s,u)*DX[r]*DX[t]);
}

double LagrangeCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   double ret=0;
   
   ret=0;
   if (DELTA(i,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 for (int t=0;t<DIM3;t++)
	 {
	    ret += ((*RefLattice_)[j][s]*U_[s][t]*Dx[t] +
		    Dx[t]*U_[t][s]*(*RefLattice_)[j][s]);
	 }
      }
      ret *= DELTA(i,p,q);
   }

   return ret;
}

double LagrangeCB::D2yDSS(int p,int q,int i,int j,int k,int l)
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
	       tmp += ((*RefLattice_)[j][s]*U_[s][t]*U_[t][r]*(*RefLattice_)[l][r] +
		       (*RefLattice_)[l][s]*U_[s][t]*U_[t][r]*(*RefLattice_)[j][r]);
	    }
	 }
      }
      tmp *= DELTA(i,p,q)*DELTA(k,p,q);
   }
   
   return tmp;
}

double LagrangeCB::D2yDUS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   double tmp=0;

   if (DELTA(k,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 tmp += ((*RefLattice_)[l][s]*U_[s][i]*DX[j] +
		 (*RefLattice_)[l][s]*U_[s][j]*DX[i] +
		 DX[i]*U_[j][s]*(*RefLattice_)[l][s] +
		 DX[j]*U_[i][s]*(*RefLattice_)[l][s]);
      }
      tmp = (0.5*DELTA(k,p,q)*(
		2.0*(*RefLattice_)[l][i]*Dx[j] + 2.0*(*RefLattice_)[l][j]*Dx[i] + tmp));
   }

   return tmp;
}

double LagrangeCB::D3yDUUS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return (0.5*DELTA(m,p,q)*(Del(i,k)*(*RefLattice_)[n][j]*DX[l]
			     + Del(i,k)*DX[j]*(*RefLattice_)[n][l]
			     + Del(i,l)*(*RefLattice_)[n][j]*DX[k]
			     + Del(i,l)*DX[j]*(*RefLattice_)[n][k]
			     + Del(j,k)*(*RefLattice_)[n][i]*DX[l]
			     + Del(j,k)*DX[i]*(*RefLattice_)[n][l]
			     + Del(j,l)*(*RefLattice_)[n][i]*DX[k]
			     + Del(j,l)*DX[i]*(*RefLattice_)[n][k]));
}

double LagrangeCB::D3yDUSS(int p,int q,int i,int j,int k,int l,int m,int n)
{
   double tmp=0;

   if (DELTA(i,p,q)*DELTA(k,p,q))
   {
      for (int s=0;s<DIM3;s++)
      {
	 tmp += ((*RefLattice_)[j][m]*U_[n][s]*(*RefLattice_)[l][s]
		 + (*RefLattice_)[j][n]*U_[m][s]*(*RefLattice_)[l][s]
		 + (*RefLattice_)[j][s]*U_[s][m]*(*RefLattice_)[l][n]
		 + (*RefLattice_)[j][s]*U_[s][n]*(*RefLattice_)[l][m]
		 + (*RefLattice_)[l][m]*U_[n][s]*(*RefLattice_)[j][s]
		 + (*RefLattice_)[l][n]*U_[m][s]*(*RefLattice_)[j][s]
		 + (*RefLattice_)[l][s]*U_[s][m]*(*RefLattice_)[j][n]
		 + (*RefLattice_)[l][s]*U_[s][n]*(*RefLattice_)[j][m]);
      }
      tmp *= 0.5*DELTA(i,p,q)*DELTA(k,p,q);
   }

   return tmp;
}

double LagrangeCB::D4yDUUSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return (0.5*DELTA(m,p,q)*DELTA(a,p,q)*
	   (Del(i,k)*(*RefLattice_)[n][j]*(*RefLattice_)[b][l]
	    + Del(i,k)*(*RefLattice_)[b][j]*(*RefLattice_)[n][l]
	    + Del(i,l)*(*RefLattice_)[n][j]*(*RefLattice_)[b][k]
	    + Del(i,l)*(*RefLattice_)[b][j]*(*RefLattice_)[n][k]
	    + Del(j,k)*(*RefLattice_)[n][i]*(*RefLattice_)[b][l]
	    + Del(j,k)*(*RefLattice_)[b][i]*(*RefLattice_)[n][l]
	    + Del(j,l)*(*RefLattice_)[n][i]*(*RefLattice_)[b][k]
	    + Del(j,l)*(*RefLattice_)[b][i]*(*RefLattice_)[n][k]));
}
