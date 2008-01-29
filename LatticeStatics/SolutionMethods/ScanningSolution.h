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
   int ModeDOFS_;
   unsigned MaxIter_;
   double Tolerance_;
   double NewtonTolerance_;
   enum YN {No,Yes};
   YN ScanFullField_;
   YN OnSolution_;
   
   Vector InitialDef_;
   
   enum ScanDir {Loading,Deformation};
   int ScnDefParam_;
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
   
   //----------------------------------------------------------------
   virtual double ScanningDefParameter();
   virtual void ScanningDefParamSet(const double val);
   virtual void ScanningDefParamUpdate(const double newval);
   virtual double ScanningLoadParameter();
   virtual void ScanningLoadParamSet(const double val);
   virtual void ScanningLoadParamUpdate(const double newval);
   virtual double ScanningStressParameter();
   
   virtual Vector ScanningForce();
   virtual Vector ScanningDef();
   virtual void ScanningSet(const Vector &val);
   virtual void ScanningUpdate(const Vector &newval);
   virtual Matrix ScanningStiffness();
   
public:
   ScanningSolution(LatticeMode *Mode,char *datafile,const char *prefix,int Echo=1);
   ~ScanningSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int FindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				 char *datafile,const char *prefix,int Width,fstream &out)
   {return 1;}
};

#endif
