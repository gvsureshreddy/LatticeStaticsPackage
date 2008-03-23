#include "LagrangeCB.h"

LagrangeCB::LagrangeCB(unsigned InternalAtoms,Matrix &RefLattice,Vector *AtomPositions)
   : CBKinematics(InternalAtoms,RefLattice,AtomPositions)
{
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms_,DIM3);
   
   SetReferenceDOFs();
   Reset();
}

LagrangeCB::LagrangeCB(PerlInput &Input,PerlInput::HashStruct *ParentHash)
   : CBKinematics(Input,ParentHash)
{
   F_.Resize(DIM3,DIM3);
   S_.Resize(InternalAtoms_,DIM3);
   
   SetReferenceDOFs();
   Reset();
}

void LagrangeCB::Reset()
{
   unsigned i,j,q,p;
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

Vector LagrangeCB::FractionalPosVec(int p)
{
   Vector pos(DIM3,0.0);
   
   for (unsigned i=0;i<DIM3;++i)
   {
      pos[i] = InternalPOS_[p][i] + S_[p][i];
   }
   
   return pos;
}

double LagrangeCB::DX(double *X,int p,int q,int i)
{
   double tmp=0.0;
   
   for (unsigned j=0;j<DIM3;++j)
   {
      tmp += (X[j] + ((InternalPOS_[q][j] + S_[q][j])
                      - (InternalPOS_[p][j] + S_[p][j])))
         *RefLattice_[j][i];
   }
   
   return tmp;
}

double LagrangeCB::Dx(double *X,int p,int q,int i)
{
   double tmp=0.0;
   
   for (unsigned k=0;k<DIM3;++k)
   {
      tmp += F_[i][k]*DX(X,p,q,k);
   }
   
   return tmp;
}

double LagrangeCB::DyDF(double *Dx,double *DX,int r, int s)
{
   return 2.0*Dx[r]*DX[s];
}

double LagrangeCB::D2yDFF(double *DX,int r, int s, int t, int u)
{
   return 2.0*Del(r,t)*DX[s]*DX[u];
}

double LagrangeCB::DyDS(double *Dx,int p,int q,int i, int j)
{
   double ret=0;
   
   ret=0;
   if (DELTA(i,p,q))
   {
      for (unsigned r=0;r<DIM3;r++)
      {
         for (unsigned k=0;k<DIM3;k++)
         {
            ret += F_[k][r]*RefLattice_[j][r]*Dx[k];
         }
      }
      ret *= 2.0*DELTA(i,p,q);
   }
   
   return ret;
}

double LagrangeCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   if (DELTA(i,p,q)*DELTA(k,p,q))
   {
      for (unsigned t=0;t<DIM3;t++)
      {
         for (unsigned r=0;r<DIM3;r++)
         {
            for (unsigned s=0;s<DIM3;s++)
            {
               tmp += F_[t][r]*RefLattice_[j][r]*F_[t][s]*RefLattice_[l][s];
            }
         }
      }
      tmp *= 2.0*DELTA(i,p,q)*DELTA(k,p,q);
   }
   
   return tmp;
}

double LagrangeCB::D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
   double tmp=0;
   
   if (DELTA(k,p,q))
   {
      for (unsigned r=0;r<DIM3;r++)
      {
         tmp += F_[i][r]*RefLattice_[l][r];
      }
      tmp = 2.0*DELTA(k,p,q)*(Dx[i]*RefLattice_[l][j] + tmp*DX[j]);
   }
   
   return tmp;
}

double LagrangeCB::D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
   return 2.0*DELTA(m,p,q)*Del(i,k)*(RefLattice_[n][j]*DX[l] + DX[j]*RefLattice_[n][l]);
}

double LagrangeCB::D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n)
{
   double tmp=0;
   
   if (DELTA(i,p,q)*DELTA(k,p,q))
   {
      for (unsigned s=0;s<DIM3;s++)
      {
         tmp += F_[m][s]*(RefLattice_[j][n]*RefLattice_[l][s]
                          + RefLattice_[j][s]*RefLattice_[l][n]);
      }
      tmp *= 2.0*DELTA(i,p,q)*DELTA(k,p,q);
   }
   
   return tmp;
}

double LagrangeCB::D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
   return (2.0*DELTA(m,p,q)*DELTA(a,p,q)*Del(i,k)*
           (RefLattice_[n][j]*RefLattice_[b][l]
            + RefLattice_[n][l]*RefLattice_[b][j]));
}
