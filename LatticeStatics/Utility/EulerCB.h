#ifndef __EulerCB
#define __EulerCB

#include <Matrix.h>
#include <Vector.h>
#include "CBKinematics.h"

using namespace std;

class EulerCB: public CBKinematics
{
private:
   
public:
   EulerCB(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS);
   virtual ~EulerCB() {};

   virtual double DX(double *X,int p,int q,int i);
   virtual double Dx(double *X,int p,int q,int i);

   virtual double DyDU(double *Dx,double *DX,int r,int s);
   virtual double D2yDUU(double *DX,int r,int s,int t,int u);
   virtual double DyDS(double *Dx,int p,int q,int i,int j);
   virtual double D2yDSS(int p,int q,int i,int j,int k, int l);
   virtual double D2yDUS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l);
   virtual double D3yDUUS(double *DX,int p,int q,int i,int j,int k,int l,int m,int n);
   virtual double D3yDUSS(int p,int q,int i,int j,int k,int l,int m,int n);
   virtual double D4yDUUSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b);
};

#endif
