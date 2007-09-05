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

   virtual void Reset();

   virtual double DX(double *X,int p,int q,int i);
   virtual double Dx(double *X,int p,int q,int i);

   virtual double DyDF(double *Dx,double *DX,int r,int s);
   virtual double D2yDFF(double *DX,int r,int s,int t,int u);
   virtual double DyDS(double *Dx,int p,int q,int i,int j);
   virtual double D2yDSS(int p,int q,int i,int j,int k, int l);
   virtual double D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l);
   virtual double D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m,int n);
   virtual double D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n);
   virtual double D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b);

   virtual char *IDString() {return "EulerCB";}

};

#endif
