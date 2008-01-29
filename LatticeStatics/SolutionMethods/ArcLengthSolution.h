#ifndef __ArcLengthSolution
#define __ArcLengthSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class ArcLengthSolution : public SolutionMethod
{
private:
   int Echo_;
   LatticeMode *Mode_;
   int ModeDOFS_;
   unsigned MaxIter_;
   double Tolerance_;
   double BisectTolerance_;
   
   double DSMax_;
   double DSMin_;
   double CurrentDS_;
   double AngleCutoff_;
   double AngleIncrease_;
   double Aspect_;
   
   unsigned NumSolutions_;
   unsigned CurrentSolution_;
   int ClosedLoopStart_;
   Vector FirstSolution_;
   
   Vector Difference_;
   
   double ArcLengthNewton(int &good);
   void ConsistencyCheck(Vector &Solution1,Vector &Solution2,
			 double ConsistencyEpsilon,int Width,fstream &out);
   virtual int OldFindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				    char *datafile,const char *prefix,int Width,fstream &out);
   
   Vector ArcLenForce(double DS,const Vector &Diff,double Aspect);
   Vector ArcLenDef() {return Mode_->ModeDOF();}
   void ArcLenSet(const Vector &val) {Mode_->SetModeDOF(val);}
   void ArcLenUpdate(const Vector &newval) {Mode_->UpdateModeDOF(newval);}
   double ArcLenAngle(Vector Old,Vector New,double Aspect);
   Matrix ArcLenStiffness(const Vector &Diff,double Aspect);
   
public:
   ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
		     const Vector &one,const Vector &two,int Echo=1);
   ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
		     char *startfile,fstream &out,int Echo=1);
   ~ArcLengthSolution() {}
   
   double GetAspect() {return Aspect_;}
   void SetCurrentDS(double ds) {CurrentDS_ = ds;}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int FindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				 char *datafile,const char *prefix,int Width,fstream &out);
};

#endif
