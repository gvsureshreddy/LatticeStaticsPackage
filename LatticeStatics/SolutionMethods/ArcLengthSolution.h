#ifndef RSE__ArcLengthSolution
#define RSE__ArcLengthSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "Restriction.h"

using namespace std;

#define CLOSEDDEFAULT 30

class Lattice;

class ArcLengthSolution : public SolutionMethod
{
public:
   enum ConvergeType {Both,Force,Displacement};
private:
   int Echo_;
   Restriction *Restrict_;
   int DOFS_;
   int MaxIter_;
   ConvergeType ConvergeType_;  //Quantities to check for convergence
   double Tolerance_;
   
   double DSMax_;
   double DSMin_;
   double CurrentDS_;
   double AngleCutoff_;
   double AngleIncrease_;
   double Aspect_;
   
   int NumSolutions_;
   int CurrentSolution_;
   int BifStartFlag_;
   Vector BifTangent_;
   int ClosedLoopStart_;
   int ClosedLoopUseAsFirst_;
   Vector FirstSolution_;
   int StopAtCPCrossingNum_;
   
   Vector Difference_;
   
   void ArcLengthNewton(int& good);
   int ZBrent(Lattice* const Lat,int const& track,Vector const& OriginalDiff,
              double const& OriginalDS,double& fa,double& fb,Vector& CurrentTF);
   
   Vector const& ArcLenForce(double const& DS,Vector const& Diff,double const& Aspect) const;
   Vector ArcLenDef() const {++counter_[3]; return Restrict_->DOF();}
   void ArcLenSet(Vector const& val) {++counter_[4]; Restrict_->SetDOF(val);}
   void ArcLenUpdate(Vector const& newval) {++counter_[5]; Restrict_->UpdateDOF(newval);}
   double ArcLenAngle(Vector const& Old,Vector const& New,double const& Aspect) const;
   Matrix const& ArcLenStiffness(Vector const& Diff,double const& Aspect) const;
   
public:
   ArcLengthSolution(Restriction* const Restrict,Vector const& dofs,
                     int const& MaxIter,double const& Tolerance,ConvergeType CnvrgTyp,
                     double const& DSMax,double const& DSMin,double const& CurrentDS,
                     double const& AngleCutoff,double const& AngleIncrease,double const& Aspect,
                     int const& NumSolutions,int const& CurrentSolution,
                     Vector const& FirstSolution,Vector const& Difference,
                     int const& BifStartFlag_,Vector const& BifTangent,
                     int const& ClosedLoopStart,int const& ClosedLoopUseAsFirst,
                     int const& StopAtCPCrossingNum,int const& Echo);
   ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,
                     Vector const& one,Vector const& two,int const& Echo=1);
   ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,int const Echo=1);
   ~ArcLengthSolution();
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const;
   virtual int FindNextSolution();
   virtual void FindCriticalPoint(Lattice* const Lat,int* const TotalNumCPCrossings,
                                  PerlInput const& Input,int const& Width,ostream& out);
   virtual char const* const Type() const {return "ArcLengthSolution";}
   
private:
   // "static" member variables
   // ArcLenForce
   mutable Vector force_static;
   mutable Vector mdfc_static;
   // ArcLenStiffness
   mutable Matrix K_static;
   mutable Matrix RestrictK_static;
   // FindCriticalPoint
   mutable Vector TF_LHS_static;
   mutable Vector TF_RHS_static;
   mutable Vector CurrentTF_static;

   // counter
   static int const nocounters_ = 11;
   mutable int counter_[nocounters_];
};

#endif
