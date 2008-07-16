#ifndef RSE__ScanningSolution
#define RSE__ScanningSolution

#include "PerlInput.h"
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
   Vector ModeDOF_;
   int MaxIter_;
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

   void ScanningNewton(int& good);
   void InitializeLine();
   
   //----------------------------------------------------------------
   virtual double const& ScanningDefParameter() const;
   virtual void ScanningDefParamSet(double const& val);
   virtual void ScanningDefParamUpdate(double const& newval);
   virtual double const& ScanningLoadParameter() const;
   virtual void ScanningLoadParamSet(double const& val);
   virtual void ScanningLoadParamUpdate(double const& newval);
   virtual double const& ScanningStressParameter() const;
   
   virtual Vector const& ScanningForce() const;
   virtual Vector const& ScanningDef() const;
   virtual void ScanningSet(Vector const& val);
   virtual void ScanningUpdate(Vector const& newval);
   virtual Matrix const& ScanningStiffness() const;
   
public:
   ScanningSolution(LatticeMode* const Mode,PerlInput const& Input,int const& Echo=1);
   ScanningSolution(LatticeMode* const Mode,
                    int const& MaxIter,double const& Tolerance,double const& NewtonTolerance,
                    YN const& ScanFullField,Vector const& InitialDef,int const& ScnDefParam,
                    ScanDir const& Direction,double const& ScanStart,double const& ScanEnd,
                    double const& ScanStep,double const& LineStart,double const& LineEnd,
                    double const& LineStep,YN const& OnSolution=No,int const& Echo=1);
   ~ScanningSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const;
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice* const Lat,PerlInput const& Input,int const& Width,
                                 fstream& out)
   {return 1;}

private:
   // "static" member variables
   // ScanningDef
   mutable Vector DEF_static;
   // ScanningForce
   mutable Vector stress_static;
   mutable Vector force_static;
   // ScanningStiffness
   mutable Matrix ModeK_static;
   mutable Matrix K_static;
};

#endif
