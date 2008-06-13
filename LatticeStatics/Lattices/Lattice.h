#ifndef __Lattice
#define __Lattice

#include <Matrix.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>

#include "PerlInput.h"

using namespace std;

#define LINELENGTH 600
#define DOFMAX 256
#define BIFMAX 25

class Lattice
{
public:
   int Echo_;
   int dbg_;
   
   enum LoadType {Temperature,Load};
   enum StateType {LHS,RHS,CRITPT};
   LoadType LoadParameter_;
   LoadType LoadParameter() {return LoadParameter_;}
   int OrderedTFs_;
   int ThirdOrder_;
   
   Lattice(PerlInput &Input);
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
   void SetLoadParameter(const double &load);
   
   virtual double E0() = 0;
   virtual Matrix E1() = 0;
   virtual Matrix E1DLoad() = 0;
   virtual Matrix E2() = 0;
   virtual Matrix E3() = 0;
   virtual Matrix E4() = 0;
   virtual int TestFunctions(Vector &TF1, StateType State = LHS, Vector *EV2= NULL);
   virtual void DispersionCurves(Vector K,int NoPTS,const char *prefix,ostream &out) {};
   virtual int BlochWave(Vector &K) {return -1;}
   virtual void LongWavelengthModuli(double dk,int gridsize,const char *prefix,
                                     ostream &out) {};
   virtual void SetParameters(double *Vals,int ResetRef = 1) = 0;
   virtual void SetGridSize(int Grid) = 0;
   virtual void NeighborDistances(int cutoff,ostream &out) {};
   virtual void CriticalPointInfo(const Vector &DrDt,int NumZeroEigenVals,
                                  double Tolerance,int Width,ostream &out);
   void ConsistencyCheck(double ConsistencyEpsilon,int Width,ostream &out);
   virtual void DebugMode() {};
   
   enum PrintDetail {PrintLong,PrintShort};
   virtual void Print(ostream &out,PrintDetail flag) = 0;
   friend ostream &operator<<(ostream &out,Lattice *L)
   {L->Print(out,PrintShort); return out;}

private:
   // "static" member variables
   // TestFunctions
   int test_flag_static;
   Matrix Stiffness_1_static;
   Matrix Stiffness_2_static;
   Matrix Stiffness_3_static;
   Matrix Stiffness_temp_static;
   Matrix Stiffness_diagonalized_static;
   Matrix EigVect_static;
   Matrix EigVectRHS_static;
   Matrix EigVectLHS_static;
   Matrix EV1_static;
   Matrix EV2_static;
};

#endif
