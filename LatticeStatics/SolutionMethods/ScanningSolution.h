#ifndef __ScanningSolution
#define __ScanningSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

class Lattice;

class ScanningSolution : public SolutionMethod
{
private:
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

   int ScanningNewton();
   void InitializeLine();

public:
   ScanningSolution(LatticeMode *Mode,char *datafile);
   ~ScanningSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual int FindNextSolution();
   virtual int BisectAlert(Lattice *Lat,int Width,fstream &out) {return 1;}
   
};

#endif
