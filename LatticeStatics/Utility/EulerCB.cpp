#include "EulerCB.h"

EulerCB::EulerCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS)
{
   DOF_=DOF;
   RefLattice_=RefLat;
   InternalAtoms_=InternalAtoms;
   InternalPOS_=InternalPOS;
   U_.Resize(DIM3,DIM3);
   V_.Resize(InternalAtoms,DIM3);
   Reset();
}

double EulerCB::DX(double *X,int p,int q,int i)
{
}

double EulerCB::Dx(double *X,int p,int q,int i)
{
}

double EulerCB::DyDU(double *Dx,double *DX,int r, int s)
{
}

double EulerCB::D2yDUU(double *DX,int r, int s, int t, int u)
{
}

double EulerCB::DyDS(double *Dx,int p,int q,int i, int j)
{
}

double EulerCB::D2yDSS(int p,int q,int i,int j,int k,int l)
{
}

double EulerCB::D2yDUS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l)
{
}

double EulerCB::D3yDUUS(double *DX,int p,int q,int i,int j,int k,int l,int m, int n)
{
}

double EulerCB::D3yDUSS(int p,int q,int i,int j,int k,int l,int m,int n)
{
}

double EulerCB::D4yDUUSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b)
{
}
