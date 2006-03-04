#ifndef __Lattice
#define __Lattice

#include <Matrix.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>

using namespace std;

#define DOFMAX 256
#define BIFMAX 25

class Lattice
{
public:
   int Echo_;
   int dbg_;

   enum LoadType {Temperature,Load};
   LoadType LoadParameter_;
   LoadType LoadParameter() {return LoadParameter_;}
   
   virtual ~Lattice() {}
   
   virtual Vector DOF() = 0;
   virtual void SetDOF(const Vector &dof) = 0;
   virtual double Entropy() = 0;
   virtual double HeatCapacity() = 0;
   virtual Matrix StressDT() = 0;
   virtual Matrix StiffnessDT() = 0;
   virtual double Temp() = 0;
   virtual void SetTemp(const double &temp) = 0;
   virtual Matrix StressDL() = 0;
   virtual Matrix StiffnessDL() = 0;
   virtual double Lambda() = 0;
   virtual void SetLambda(const double &lambda) = 0;
   
   
   virtual double Energy() = 0;
   double E0() {return Energy();}
   virtual Matrix Stress() = 0;
   Matrix E1() {return Stress();}
   virtual Matrix Stiffness() = 0;
   Matrix E2() {return Stiffness();}
   virtual Matrix E3() = 0;
   virtual Matrix E4() = 0;
   virtual int StiffnessNulity(double *Min=NULL);
   virtual void DispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out) {};
   virtual int BlochWave(Vector &K) {return -1;}
   virtual void LongWavelengthModuli(double dk,int gridsize,const char *prefix,
				     ostream &out) {};
   virtual void SetGridSize(int Grid) = 0;
   virtual void NeighborDistances(int cutoff,ostream &out) {};
   virtual void CriticalPointInfo(const Vector &DrDt,int NumZeroEigenVals,
				  double Tolerance,char *datafile,const char *prefix,
				  int Width,ostream &out);
   void ConsistencyCheck(double ConsistencyEpsilon,int Width,ostream &out);
   virtual void DebugMode() {};

   enum PrintDetail {PrintLong,PrintShort};
   virtual void Print(ostream &out,PrintDetail flag) = 0;
   friend ostream &operator<<(ostream &out,Lattice *L)
   {L->Print(out,PrintShort); return out;}
};

#endif
