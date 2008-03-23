#ifndef __ArcLengthSolution
#define __ArcLengthSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

#define CLOSEDDEFAULT 30

class Lattice;

class ArcLengthSolution : public SolutionMethod
{
private:
   int Echo_;
   LatticeMode *Mode_;
   unsigned ModeDOFS_;
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
   unsigned ClosedLoopStart_;
   Vector FirstSolution_;
   
   Vector Difference_;
   
   void ArcLengthNewton(int &good);
   void ConsistencyCheck(Vector &Solution1,Vector &Solution2,
                         double ConsistencyEpsilon,int Width,fstream &out);
   virtual int OldFindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
                                    PerlInput &Input,int Width,fstream &out);
   void ZBrent(Lattice *Lat,int track,double fa,double fb,const Vector &OriginalDiff,
               const double OriginalDS,Vector &CurrentTF);
   
   Vector ArcLenForce(double DS,const Vector &Diff,double Aspect);
   Vector ArcLenDef() {return Mode_->ModeDOF();}
   void ArcLenSet(const Vector &val) {Mode_->SetModeDOF(val);}
   void ArcLenUpdate(const Vector &newval) {Mode_->UpdateModeDOF(newval);}
   double ArcLenAngle(Vector Old,Vector New,double Aspect);
   Matrix ArcLenStiffness(const Vector &Diff,double Aspect);
   
public:
   ArcLengthSolution(LatticeMode *Mode,const Vector &dofs,
                     unsigned MaxIter,double Tolerance,double BisectTolerance,double DSMax,
                     double DSMin,double CurrentDS,double AngleCutoff,double AngleIncrease,
                     double Aspect,unsigned NumSolutions,unsigned CurrentSolution,
                     const Vector &FirstSolution,const Vector &Difference,
                     unsigned ClosedLoopStart=CLOSEDDEFAULT,int Echo=1);
   ArcLengthSolution(LatticeMode *Mode,PerlInput &Input,
                     const Vector &one,const Vector &two,int Echo=1);
   ArcLengthSolution(LatticeMode *Mode,PerlInput &Input,int Echo=1);
   ~ArcLengthSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out);
};

#endif
