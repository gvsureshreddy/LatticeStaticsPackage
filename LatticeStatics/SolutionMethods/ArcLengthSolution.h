#ifndef RSE__ArcLengthSolution
#define RSE__ArcLengthSolution

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
   int ModeDOFS_;
   int MaxIter_;
   double Tolerance_;
   double BisectTolerance_;
   
   double DSMax_;
   double DSMin_;
   double CurrentDS_;
   double AngleCutoff_;
   double AngleIncrease_;
   double Aspect_;
   
   int NumSolutions_;
   int CurrentSolution_;
   int ClosedLoopStart_;
   Vector FirstSolution_;
   int StopAtCPNum_;
   int TotalNumCPs_;
   
   Vector Difference_;
   
   void ArcLengthNewton(int& good);
   void ConsistencyCheck(Vector const& Solution1,Vector const& Solution2,
                         double const& ConsistencyEpsilon,int const& Width,fstream& out);
   virtual int OldFindCriticalPoint(int const& LHN,double const& LHEV,int const& RHN,
                                    double const& RHEV,Lattice* const Lat,
                                    PerlInput const& Input,int const& Width,fstream& out);
   void ZBrent(Lattice* const Lat,int const& track,Vector const& OriginalDiff,
               double const& OriginalDS,double& fa,double& fb,Vector& CurrentTF);
   
   Vector const& ArcLenForce(double const& DS,Vector const& Diff,double const& Aspect) const;
   Vector ArcLenDef() const {return Mode_->ModeDOF();}
   void ArcLenSet(Vector const& val) {Mode_->SetModeDOF(val);}
   void ArcLenUpdate(Vector const& newval) {Mode_->UpdateModeDOF(newval);}
   double ArcLenAngle(Vector const& Old,Vector const& New,double const& Aspect) const;
   Matrix const& ArcLenStiffness(Vector const& Diff,double const& Aspect) const;
   
public:
   ArcLengthSolution(LatticeMode* const Mode,Vector const& dofs,
                     int const& MaxIter,double const& Tolerance,double const& BisectTolerance,
                     double const& DSMax,double const& DSMin,double const& CurrentDS,
                     double const& AngleCutoff,double const& AngleIncrease,
                     double const& Aspect,int const& NumSolutions,int const& CurrentSolution,
                     Vector const& FirstSolution,Vector const& Difference,
                     int const& ClosedLoopStart=CLOSEDDEFAULT,int const& StopAtCPNum=-1,
                     int const& Echo=1);
   ArcLengthSolution(LatticeMode* const Mode,PerlInput const& Input,
                     Vector const& one,Vector const& two,int const& Echo=1);
   ArcLengthSolution(LatticeMode* const Mode,PerlInput const& Input,int const Echo=1);
   ~ArcLengthSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const;
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice* const Lat,PerlInput const& Input,int const& Width,
                                 fstream& out);

private:
   // "static" member variables
   // ArcLenForce
   mutable Vector force_static;
   mutable Vector mdfc_static;
   // ArcLenStiffness
   mutable Matrix K_static;
   mutable Matrix ModeK_static;
   // FindCriticalPoint
   mutable Vector TF_LHS_static;
   mutable Vector TF_RHS_static;
   mutable Vector CurrentTF_static;
};

#endif
