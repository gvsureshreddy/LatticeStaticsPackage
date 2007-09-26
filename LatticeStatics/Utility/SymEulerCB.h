#ifndef __SymEulerCB
#define __SymEulerCB

#include <Matrix.h>
#include <Vector.h>
#include "CBKinematics.h"

using namespace std;

class SymEulerCB: public CBKinematics
{
private:
   virtual void Reset();
   
public:
   SymEulerCB(int InternalAtoms,const char* prefix,const char* datafile);
   virtual ~SymEulerCB() {};

   #include "UwithoutTransMapping.def"

   virtual Vector FractionalPosVec(int p);
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

   virtual char *IDString() {return "SymEulerCB";}

};

#endif
