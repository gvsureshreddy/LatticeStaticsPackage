#ifndef __CBKinematics
#define __CBKinematics

#include <Matrix.h>
#include <Vector.h>

using namespace std;

class CBKinematics
{
private:
public:
   const static int DIM3 = 3;
   
   Vector *DOF_;
   Matrix *RefLattice_;
   int InternalAtoms_;
   Vector *InternalPOS_;

   Matrix U_;
   Matrix V_;
   
   virtual ~CBKinematics() {};

   virtual void Reset();

   virtual void InfluenceRegion(double *InfluenceRegion);

   virtual double DX(double *X,int p,int q,int i) = 0;
   virtual double Dx(double *X,int p,int q,int i) = 0;

   virtual double DyDU(double *Dx,double *DX,int r,int s) = 0;
   virtual double D2yDUU(double *DX,int r,int s,int t,int u) = 0;
   virtual double DyDS(double *Dx,int p,int q,int i,int j) = 0;
   virtual double D2yDSS(int p,int q,int i,int j,int k, int l) = 0;
   virtual double D2yDUS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l) = 0;
   virtual double D3yDUUS(double *DX,int p,int q,int i,int j,int k,int l,int m,int n) = 0;
   virtual double D3yDUSS(int p,int q,int i,int j,int k,int l,int m,int n) = 0;
   virtual double D4yDUUSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b) = 0;

   inline double Del(int i,int j) {return i==j;}
   inline double DELTA(int s,int p,int q) {return Del(s,q) - Del(s,p);}
};

#endif
