#ifndef __ScanningSolution
#define __ScanningSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class ScanningSolution : public SolutionMethod
{
private:
   int Echo_;
   LatticeMode *Mode_;
   unsigned MaxIter_;
   double Tolerance_;
   double NewtonTolerance_;
   enum YN {No,Yes};
   YN ScanFullField_;
   YN OnSolution_;

   Vector InitialDef_;

   enum ScanDir {Loading,Deformation};
   ScanDir Direction_;
   double ScanStart_;
   double ScanEnd_;
   double ScanStep_;

   double LineStart_;
   double LineEnd_;
   double LineStep_;

   double CurrentScanLine_;

   double ScanningNewton(int &good);
   void InitializeLine();

public:
   ScanningSolution(LatticeMode *Mode,char *datafile,const char *prefix,int Echo=1);
   ~ScanningSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(int LHN,int RHN,Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out) {return 1;}
   
};

#endif
