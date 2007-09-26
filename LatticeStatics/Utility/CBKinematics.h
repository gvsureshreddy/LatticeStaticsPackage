#ifndef __CBKinematics
#define __CBKinematics

#include <Matrix.h>
#include <Vector.h>
#include "UtilityFunctions.h"

using namespace std;

class CBKinematics
{
private:
   virtual void Reset() = 0;

public:
   const static int DIM3 = 3;
   
   Vector DOF_;
   Matrix RefLattice_;
   int InternalAtoms_;
   Vector *InternalPOS_;

   Matrix F_;
   Matrix S_;

   CBKinematics(int InternalAtoms,const char* prefix,const char* datafile);
   virtual ~CBKinematics() {delete [] InternalPOS_;}

   virtual void InfluenceRegion(double *InfluenceRegion);

   virtual void SetReferenceToCurrent();
   virtual void SetReferenceDOFs();
   const Matrix& RefLattice() {return RefLattice_;}
   const double RefVolume() {return RefLattice_.Det();}
   const Vector &DOF() {return DOF_;}
   void SetDOF(const Vector &dof) {DOF_ = dof; Reset();}
   const Vector& AtomPositions(int i) {return InternalPOS_[i];}
   
   virtual inline int DOFS() {return Fsize() + Ssize();}
   virtual int Fsize() = 0;
   virtual int Ssize() = 0;
   virtual int NoTrans() = 0;
   virtual int INDF(int i,int j) = 0;
   virtual int INDS(int i,int j) = 0;
   virtual int INDFF(int k,int l,int m,int n) = 0;
   virtual int INDSS(int k,int l,int m,int n) = 0;
   virtual int INDFS(int i,int j,int m,int n) = 0;
   virtual int INDSF(int m,int n,int i,int j) = 0;

   virtual Vector CurrentLatticeVec(int p);
   virtual Vector FractionalPosVec(int p) {return Vector(DIM3,0.0);}
   virtual double DX(double *X,int p,int q,int i) = 0;
   virtual double Dx(double *X,int p,int q,int i) = 0;

   virtual double DyDF(double *Dx,double *DX,int r,int s) = 0;
   virtual double D2yDFF(double *DX,int r,int s,int t,int u) = 0;
   virtual double DyDS(double *Dx,int p,int q,int i,int j) = 0;
   virtual double D2yDSS(int p,int q,int i,int j,int k, int l) = 0;
   virtual double D2yDFS(double *Dx,double *DX,int p,int q,int i,int j,int k,int l) = 0;
   virtual double D3yDFFS(double *DX,int p,int q,int i,int j,int k,int l,int m,int n) = 0;
   virtual double D3yDSSF(int p,int q,int i,int j,int k,int l,int m,int n) = 0;
   virtual double D4yDFFSS(int p,int q,int i,int j,int k,int l,int m,int n,int a,int b) = 0;

   inline double Del(int i,int j) {return i==j;}
   inline double DELTA(int s,int p,int q) {return Del(s,q) - Del(s,p);}

   virtual char *IDString() = 0;
   friend ostream &operator<<(ostream &out,CBKinematics &CBK)
   {out << CBK.IDString(); return out;}
};

#endif
