#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include "ArcLengthSolution.h"

using namespace std;

#define ARCLENEPS 1.0e-15

extern "C" void qcbfb_output_(int& nfree,double* u,double& prop,int& nint,int* intdata,int& ndouble,double* doubledata);

ArcLengthSolution::~ArcLengthSolution()
{
   cout.width(0);
   cout << "ArcLengthSolution Function Calls:\n"
        << "\tArcLengthNewton - " << counter_[0] << "\n"
        << "\tZBrent - " << counter_[1] << "\n"
        << "\tArcLenForce - " << counter_[2] << "\n"
        << "\tArcLenDef - " << counter_[3] << "\n"
        << "\tArcLenSet - " << counter_[4] << "\n"
        << "\tArcLenUpdate - " << counter_[5] << "\n"
        << "\tArcLenAngle - " << counter_[6] << "\n"
        << "\tArcLenStiffness - " << counter_[7] << "\n"
        << "\tAllSolutionsFound - " << counter_[8] << "\n"
        << "\tFindNextSolution - " << counter_[9] << "\n"
        << "\tFindCriticalPoint - " << counter_[10] << "\n";
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,Vector const& dofs,
                                     int const& MaxIter,double const& Tolerance,
                                     ConvergeType CnvrgTyp,double const& DSMax,
                                     double const& DSMin,double const& CurrentDS,
                                     double const& AngleCutoff,double const& AngleIncrease,
                                     double const& Aspect,int const& NumSolutions,
                                     int const& CurrentSolution,Vector const& FirstSolution,
                                     Vector const& Difference,int const& BifStartFlag,
                                     Vector const& BifTangent,int const& ClosedLoopStart,
                                     int const& ClosedLoopUseAsFirst,
                                     int const& StopAtCPCrossingNum,int const& Echo)
   : Echo_(Echo),
     Restrict_(Restrict),
     DOFS_(Restrict_->DOF().Dim()),
     MaxIter_(MaxIter),
     ConvergeType_(CnvrgTyp),
     Tolerance_(Tolerance),
     DSMax_(DSMax),
     DSMin_(DSMin),
     CurrentDS_(CurrentDS),
     AngleCutoff_(AngleCutoff),
     AngleIncrease_(AngleIncrease),
     Aspect_(Aspect),
     NumSolutions_(NumSolutions_),
     CurrentSolution_(CurrentSolution),
     BifStartFlag_(BifStartFlag),
     BifTangent_(BifTangent),
     ClosedLoopStart_(ClosedLoopStart),
     ClosedLoopUseAsFirst_(ClosedLoopUseAsFirst),
     FirstSolution_(FirstSolution),
     StopAtCPCrossingNum_(StopAtCPCrossingNum),
     Difference_(Difference),
     force_static(DOFS_),
     mdfc_static(DOFS_-1),
     K_static(DOFS_,DOFS_),
     RestrictK_static(DOFS_-1,DOFS_)
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
   ArcLenSet(dofs);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,
                                     Vector const& one,Vector const& two,int const& Echo)
   : Echo_(Echo),
     Restrict_(Restrict),
     CurrentSolution_(0),
     BifStartFlag_(0),
     BifTangent_(),
     Difference_(two-one)
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
   
   DOFS_=Restrict_->DOF().Dim();
   // initialize "static" members variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_-1);
   K_static.Resize(DOFS_,DOFS_);
   RestrictK_static.Resize(DOFS_-1,DOFS_);

   
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash,"MaxIterations");
   if (Input.ParameterOK(Hash,"ConvergeType"))
   {
      char const* const cnvrgtyp=Input.getString(Hash,"ConvergeType");
      if (!strcmp("Both",cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force",cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement",cnvrgtyp))
      {
         ConvergeType_ = Displacement;
      }
      else
      {
         cerr << "Unknown ConvergeType: " << cnvrgtyp << "\nExiting!\n";
         exit(-22);
      }
   }
   else
   {
      Input.useString("Both",Hash,"ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
   Tolerance_ = Input.getDouble(Hash,"Tolerance");
   DSMax_ = Input.getDouble(Hash,"DSMax");
   CurrentDS_ = Input.getDouble(Hash,"DSStart");
   DSMin_ = Input.getDouble(Hash,"DSMin");
   AngleCutoff_ = Input.getDouble(Hash,"AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash,"AngleIncrease");
   Aspect_ = Input.getDouble(Hash,"Aspect");
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT,Hash,"ClosedLoopStart"); // Default Value
   }
   if (Input.ParameterOK(Hash,"StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash,"StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1,Hash,"StopAtCPCrossingNum"); // Default Value
   }

   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
   {
      Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
      Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
      FirstSolution_ = Restrict_->RestrictDOF(clfs);
   }
   else
   {
      if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
      {
         ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
         if (ClosedLoopUseAsFirst_ == 0)
         {
            FirstSolution_ = ArcLenDef();
         }
         else if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
         {
            cerr << "Error: ArcLengthSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                 << endl;
            exit(-33);
         }
      }
      else
      {
         // Default Value
         ClosedLoopUseAsFirst_ = Input.usePosInt(0,"StartType","ClosedLoopUseAsFirst");
         FirstSolution_ = ArcLenDef();
      }
   }
   Input.EndofInputSection();

   // Set Lattice to solution "two"
   ArcLenSet(two);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,
                                     int const Echo)
   :  Echo_(Echo),
      Restrict_(Restrict),
      CurrentSolution_(0),
      BifStartFlag_(0),
      BifTangent_()
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
   
   DOFS_=Restrict_->DOF().Dim();
   // initialize "static" memver variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_-1);
   K_static.Resize(DOFS_,DOFS_);
   RestrictK_static.Resize(DOFS_-1,DOFS_);

   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash,"MaxIterations");
   if (Input.ParameterOK(Hash,"ConvergeType"))
   {
      char const* const cnvrgtyp=Input.getString(Hash,"ConvergeType");
      if (!strcmp("Both",cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force",cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement",cnvrgtyp))
      {
         ConvergeType_ = Displacement;
      }
      else
      {
         cerr << "Unknown ConvergeType: " << cnvrgtyp << "\nExiting!\n";
         exit(-22);
      }
   }
   else
   {
      Input.useString("Both",Hash,"ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
   Tolerance_ = Input.getDouble(Hash,"Tolerance");
   DSMax_ = Input.getDouble(Hash,"DSMax");
   CurrentDS_ = Input.getDouble(Hash,"DSStart");
   DSMin_ = Input.getDouble(Hash,"DSMin");
   AngleCutoff_ = Input.getDouble(Hash,"AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash,"AngleIncrease");
   Aspect_ = Input.getDouble(Hash,"Aspect");
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT,Hash,"ClosedLoopStart"); // Default Value
   }
   if (Input.ParameterOK(Hash,"StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash,"StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1,Hash,"StopAtCPCrossingNum"); // Default Value
   }
   Input.EndofInputSection();
   
   const char *starttype = Input.getString("StartType","Type");

   if (!strcmp("Bifurcation",starttype))
   {
      //Bifurcation
      BifStartFlag_ = 1;
      BifTangent_.Resize(DOFS_);
      
      // Set Difference and Lattice state
      double eps = Input.getDouble("StartType","Epsilon");

      Difference_.Resize(DOFS_);
      Vector diff(Input.getArrayLength("StartType","Tangent"));
      Input.getVector(diff,"StartType","Tangent");
      Difference_ = Restrict_->TransformVector(diff);
      Difference_ *= eps;
      // set BifTangent_ for projection output
      BifTangent_ = Difference_/Difference_.Norm();
      BifTangent_[DOFS_-1] = 0.0;
      
      Vector stat(Input.getArrayLength("StartType","BifurcationPoint"));
      Input.getVector(stat,"StartType","BifurcationPoint");
      // Set Lattice state to the bifurcation point
      Vector RestrictedStat = Restrict_->RestrictDOF(stat);
      ArcLenSet(RestrictedStat);
      
      // Set FirstSolution
      FirstSolution_.Resize(DOFS_);
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
         Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
            if (ClosedLoopUseAsFirst_ == 0)
            {
               FirstSolution_ = ArcLenDef();
            }
            else if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
            {
               cerr << "Error: ArcLengthSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                    << endl;
               exit(-33);
            }
         }
         else
         {
            // Default Value
            ClosedLoopUseAsFirst_ = Input.useInt(0,"StartType","ClosedLoopUseAsFirst");
            FirstSolution_ = ArcLenDef();
         }
      }

      cout << "Projection on BifTangent of BifurcationPoint = " << RestrictedStat*BifTangent_ << "\n";
   }
   else if (!strcmp("Continuation",starttype))
   {
      // Get solution1
      Vector one(DOFS_);
      Vector onetmp(Input.getArrayLength("StartType","Solution1"));
      Input.getVector(onetmp,"StartType","Solution1");
      one = Restrict_->RestrictDOF(onetmp);

      // Set Lattice state to Solution2
      Vector two(DOFS_);
      Vector twotmp(Input.getArrayLength("StartType","Solution2"));
      Input.getVector(twotmp,"StartType","Solution2");
      two = Restrict_->RestrictDOF(twotmp);
      ArcLenSet(two);
      
      // Set Difference_ to   two - one
      Difference_.Resize(DOFS_);
      Difference_ = two - one;
      
      // Set FirstSolution
      FirstSolution_.Resize(DOFS_);
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
         Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
            if (ClosedLoopUseAsFirst_ == 0)
            {
               FirstSolution_ = ArcLenDef();
            }
            else if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
            {
               cerr << "Error: ArcLengthSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                    << endl;
               exit(-33);
            }
         }
         else
         {
            // Default Value
            ClosedLoopUseAsFirst_ = Input.useInt(0,"StartType","ClosedLoopUseAsFirst");
            FirstSolution_ = ArcLenDef();
         }
      }
   }
   else if (!strcmp("ConsistencyCheck",starttype))
   {
      double ConsistencyEpsilon;
      int Width;
      Vector Solution(DOFS_);

      Vector onetmp(Input.getArrayLength("StartType","Solution"));
      Input.getVector(onetmp,"StartType","Solution");
      Solution = Restrict_->RestrictDOF(onetmp);
      // Get Epsilon and Width
      ConsistencyEpsilon = Input.getDouble("StartType","Epsilon");
      Width = Input.getPosInt("Main","FieldWidth");

      ostream::fmtflags oldflags=cout.flags();
      cout << scientific;
      Restrict_->ConsistencyCheck(Solution,ConsistencyEpsilon,Width,cout);
      cout.flags(oldflags);
      // We're done
      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      cerr << "Unknown StartType!" << "\n";
      exit(-1);
   }
   Input.EndofInputSection();
}

Vector const& ArcLengthSolution::ArcLenForce(double const& DS,Vector const& Diff,
                                             double const& Aspect) const
{
   ++counter_[2];
   
   mdfc_static = Restrict_->Force();
   
   force_static[DOFS_-1] = DS*DS - Diff[DOFS_-1]*Diff[DOFS_-1]/(Aspect*Aspect);
   for (int i=0;i<DOFS_-1;++i)
   {
      force_static[i] = mdfc_static[i];
      force_static[DOFS_-1] -= Diff[i]*Diff[i];
   }
   
   return force_static;
}

Matrix const& ArcLengthSolution::ArcLenStiffness(Vector const& Diff,double const& Aspect)
   const
{
   ++counter_[7];
   
   RestrictK_static = Restrict_->Stiffness();
   
   for (int i=0;i<DOFS_-1;++i)
   {
      for (int j=0;j<=DOFS_-1;++j)
      {
         K_static[i][j] = RestrictK_static[i][j];
      }
      K_static[DOFS_-1][i] = -2.0*Diff[i];
   }
   K_static[DOFS_-1][DOFS_-1] = -2.0*Diff[DOFS_-1]/(Aspect*Aspect);
   
   return K_static;
}

double ArcLengthSolution::ArcLenAngle(Vector const& Old,Vector const& New,double const& Aspect)
   const
{
   ++counter_[6];
   
   double angle = 0.0;
   double NewNorm = 0.0;
   double OldNorm = 0.0;

   for (int i=0;i<DOFS_-1;++i)
   {
      angle += Old[i]*New[i];
      NewNorm += New[i]*New[i];
      OldNorm += Old[i]*Old[i];
   }

   angle += Old[DOFS_-1]*New[DOFS_-1]/(Aspect*Aspect);
   NewNorm += New[DOFS_-1]*New[DOFS_-1]/(Aspect*Aspect);
   OldNorm += Old[DOFS_-1]*Old[DOFS_-1]/(Aspect*Aspect);
   NewNorm = sqrt(NewNorm);
   OldNorm = sqrt(OldNorm);
   
   return fabs(acos( angle/(NewNorm*OldNorm) ));
}

int ArcLengthSolution::AllSolutionsFound() const
{
   ++counter_[8];
   
   return (CurrentSolution_ >= NumSolutions_);
}

int ArcLengthSolution::FindNextSolution()
{
   ++counter_[9];
   
   int good=0;
   double AngleTest;
   // Assume that first solution should not be strictly restricted to the
   // adaptive steping angle constraint
   double AngleFactor = CurrentSolution_ ? 1.0 : 5.0;
   
   Vector OldDiff = Difference_;
   
   do
   {
      ArcLengthNewton(good);
      
      AngleTest = ArcLenAngle(OldDiff,Difference_,Aspect_);
      
      cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << "\n";
   }
   while (((AngleTest >= AngleFactor*AngleCutoff_) || !good)
          && (CurrentDS_ >= DSMin_)
          && (ArcLenUpdate(-Difference_),// back to previous solution
              Difference_ = OldDiff,
              CurrentDS_=CurrentDS_/2.0));
   
   if ((AngleTest <= AngleIncrease_) && (CurrentDS_ < DSMax_))
   {
      CurrentDS_ *= 2.0;
      if (CurrentDS_ > DSMax_) CurrentDS_ = DSMax_;
   }
   
   if (!good)
   {
      cerr << "ArcLenghtSolution did not converge properly" << "\n";
   }
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((ArcLenDef() - FirstSolution_).Norm() < CurrentDS_))
   {
      // We are done -- set currentsolution to numsolutions
      cout << "NOTE: Closed Loop detected at Solution # " << CurrentSolution_
           << " --- Terminating!" << "\n";
      
      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      CurrentSolution_++;
   }
   if (CurrentSolution_ == ClosedLoopUseAsFirst_)
   {
      // set First Solution for use with Closed Loop check.
      FirstSolution_ = ArcLenDef();
   }
   
   // Always have the current "solution" state printed as a solution point
   good = 1;
   
   return good;
}

void ArcLengthSolution::ArcLengthNewton(int& good)
{
   ++counter_[0];
   
   int itr = 0;
   int Dim=ArcLenDef().Dim();
   double ForceNorm = 0.0;
   double Magnitude1 = 0.0;
   double Magnitude2 = 0.0;
   double Kappa = 0.0;
   double const eta = 0.1;
   
   Vector Dx(Dim),
      RHS(Dim);
   Matrix stif(Dim,Dim);
   
   // Predictor step
   cout << "Taking Predictor Step. CurrentDS = " << CurrentDS_ << "\n";
   ArcLenUpdate(Difference_);
   
   // Iterate until convergence

   itr++;
   // get stiffness first for efficiency
   stif=ArcLenStiffness(Difference_,Aspect_);
   RHS = -ArcLenForce(CurrentDS_,Difference_,Aspect_);
   ForceNorm = RHS.Norm();
#ifdef SOLVE_SVD
   Dx = SolveSVD(
      stif,
      RHS,MAXCONDITION,Echo_);
#else
   Dx = SolvePLU(stif,RHS);
#endif
   Magnitude1 = Dx.Norm();
   Magnitude2 = Magnitude1;

   cout << "\tForceNorm = " << ForceNorm << " \tDeltaNorm = " << Magnitude2 << "\n";

   int Converged = 0;
   do
   {
      itr++;

      ArcLenUpdate(Dx);
      Difference_ += Dx;
      // get stiffness first for efficiency
      stif=ArcLenStiffness(Difference_,Aspect_);
      RHS = -ArcLenForce(CurrentDS_,Difference_,Aspect_);
      ForceNorm = RHS.Norm();
#ifdef SOLVE_SVD
      Dx = SolveSVD(
         stif,
         RHS,MAXCONDITION,Echo_);
#else
      Dx = SolvePLU(stif,RHS);
#endif
      Magnitude2 = Dx.Norm();
      Kappa = Magnitude2/(Magnitude1+Tolerance_*eta);
      Magnitude1 = Magnitude2;
      
      cout << "\tForceNorm = " << ForceNorm
           << " \tDeltaNorm = " << Magnitude2
           << " \tContraction = " << Kappa << "\n";

      switch (ConvergeType_)
      {
         case Both:
            if ((RHS.Norm() <= Tolerance_) && (Dx.Norm() <= Tolerance_))
            {
               Converged = 1;
            }
            break;
         case Force:
            if (RHS.Norm() <= Tolerance_)
            {
               Converged = 1;
            }
            break;
         case Displacement:
            if (Dx.Norm() <= Tolerance_)
            {
               Converged = 1;
            }
            break;
      }
   }
   while ((itr < MaxIter_) && (!Converged));
   
   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << "\n";
      good = 0;
   }
   else
   {
      cout << "Prediction 1 Corrector Iterations: " << itr << "\n";

      if (BifStartFlag_)
      {
         Vector diff = Difference_;
         diff[Dim-1] = 0.0;
         diff /= Dx.Norm();
         cout << "Projection on BifTangent = " << Restrict_->DOF()*BifTangent_
              << ",     Angle (deg.) with BifTangent = "
              << acos(diff*BifTangent_)*(57.2957795130823) << "\n";
      }

      cout << "Converged with ForceNorm = " << RHS.Norm()
           << ",     CorrectorNorm = " << Dx.Norm() << "\n";
      
      good = 1;
   }
}

void ArcLengthSolution::FindCriticalPoint(Lattice* const Lat,int& TotalNumCPCrossings,
                                          PerlInput const& Input,int const& Width,ostream& out)
{
   ++counter_[10];
   
   Vector OriginalDOF=ArcLenDef();
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   int TestValueDiff;
   int temp;
   int size = Lat->NumTestFunctions();
   TF_LHS_static.Resize(size);
   TF_RHS_static.Resize(size);
   CurrentTF_static.Resize(size);
   double fa,fb;
   double LHSLambda,RHSLambda;
   int Multiplicity;
   int track;
   int num;
   int CP;
   int spot;
   int CPorBif;
   int Bif;

   RHSLambda = ArcLenDef()[DOFS_-1];
   LHSLambda = RHSLambda - Difference_[DOFS_-1];
   
   TestValueDiff = Lat->TestFunctions(TF_LHS_static, Lattice::RHS, &TF_RHS_static);
   if (TestValueDiff < 0)
   {
      out << "Note: TestFunctions found a discrepancy between the\n"
          << "Note: difference in number of negative Test Functions\n"
          << "Note: and the number of Test Functions that change sign\n"
          << "Note: from LeftHandSide to RightHandSide.  This is usually\n"
          << "Note: caused by having too large of a path-following stepsize."
          << "\n";
      if (Echo_)
         cout << "Note: TestFunctions found a discrepancy between the\n"
              << "Note: difference in number of negative Test Functions\n"
              << "Note: and the number of Test Functions that change sign\n"
              << "Note: from LeftHandSide to RightHandSide.  This is usually\n"
              << "Note: caused by having too large of a path-following stepsize."
              << "\n";
      TestValueDiff = -TestValueDiff;
   }
   
   int* Index;
   Index = new int[TestValueDiff];
   Vector DSTrack(TestValueDiff);
   Vector* CPDOFs;
   double* CPLambdas;
   CPDOFs = new Vector[TestValueDiff];
   for (int i=0;i<TestValueDiff;++i) CPDOFs[i].Resize(Lat->DOF().Dim());
   CPLambdas = new double[TestValueDiff];
   int* CPorBifs;
   CPorBifs = new int[TestValueDiff];

   temp = 0;
   for (int i = 0; i< size; i++)
   {
      if ((TF_LHS_static[i]*TF_RHS_static[i]) < 0.0)
      {
         Index[temp] = i;
         temp++;
      }
   }
   
   cout << "TestValueDiff = "<< TestValueDiff << "\n";
   for (int i = 0; i < TestValueDiff; i++)
   {
      cout << "Index[" << setw(2) << i << "] = " << setw(6) << Index[i]
           << ",   TF_LHS[" << setw(6) << Index[i] << "] = "
           << setw(Width) << TF_LHS_static[Index[i]]
           << ",   TF_RHS[" << setw(6) << Index[i] << "] = "
           << setw(Width) << TF_RHS_static[Index[i]] << "\n";
   }
   
   num = 0;
   for (CP= 0; CP < TestValueDiff; CP++)
   {
      track = Index[CP];
      fa = TF_LHS_static[track];
      fb = TF_RHS_static[track];
      
      if(track>=0) //START OF IF STATEMENT
      {
         ZBrent(Lat,track,OriginalDiff,OriginalDS,fa,fb,CurrentTF_static);
         double lambda = ArcLenDef()[DOFS_-1];
         if ((Lat->UseEigenValTFs() == 0) || (track >= Lat->DOF().Dim()))
         {
            CPorBif = -1; // ExtraTF
         }
         else if (((lambda >= LHSLambda) && (lambda <= RHSLambda))
                  || ((lambda <= LHSLambda) && (lambda >= RHSLambda)))
         {
            CPorBif = 1; // bif point
         }
         else
         {
            CPorBif = 0; // turning point
         }
         
         Multiplicity = 1;
         for(int i=CP+1;i<TestValueDiff;i++)
         {
            temp = Index[i];
            if(fabs(CurrentTF_static[temp]) <= Tolerance_)
            {
               Index[i] = -1;
               Multiplicity++;
               //cout <<"i = " << i << "\n" <<  "CHECK POINT INDEX[CP] = " << Index[i] << "\n";
            }
         }
         
         // sort the critical points
         spot=num;
         while((spot!=0) && (DSTrack[spot-1] > CurrentDS_))
         {
            DSTrack[spot] = DSTrack[spot-1];
            CPDOFs[spot] = CPDOFs[spot - 1];
            CPLambdas[spot] = CPLambdas[spot - 1];
            CPorBifs[spot] = CPorBifs[spot - 1];
            spot = spot - 1;
         }
         DSTrack[spot] = CurrentDS_;
         CPDOFs[spot] = Lat->DOF();
         CPLambdas[spot] = ( (Lat->LoadParameter() == Lattice::Load)
                            ? Lat->Lambda() : Lat->Temp() );
         CPorBifs[spot] = CPorBif;

         num = num + 1;
      }//END OF IF STATEMENT
   }
   
   ////PRINT OUT CP DATA
   for (int i = 0; i < num; i++)
   {
      // reset to appropriate cp.
      Lat->SetDOF(CPDOFs[i]);
      if (Lat->LoadParameter() == Lattice::Load)
         Lat->SetLambda(CPLambdas[i]);
      else
         Lat->SetTemp(CPLambdas[i]);
      
      // Output Critical Point
      for (int j=0;j<70;j++)
      {
         if (Echo_) cout << "=";
         out << "=";
      }
      if (Echo_) cout << "\n";
      out << "\n";
      
      // Lattice takes care of echo
      out << setw(Width);
      if (0 == CPorBifs[i])
         Lat->Print(out,Lattice::PrintShort,Lattice::TurningPt);
      else if (1 == CPorBifs[i])
         Lat->Print(out,Lattice::PrintShort,Lattice::BifurcationPt);
      else
         Lat->Print(out,Lattice::PrintShort,Lattice::ExtraTFPt);
      
      for (int j=0;j<70;j++)
      {
         if (Echo_) cout << "=";
         out << "=";
      }
      if (Echo_) cout << "\n";
      out << "\n";
      
      // Call Lattice function to do any Lattice Specific things
      Bif=Lat->CriticalPointInfo(TotalNumCPCrossings,
                                 Restrict_->DrDt(Restrict_->DOF()-(OriginalDOF-OriginalDiff)),
                                 CPorBifs[i],Multiplicity,10.0*Tolerance_,Width,Input,out);
      
      if (Echo_) cout << "Success = 1" << "\n";
      out << "Success = 1" << "\n";

      ++TotalNumCPCrossings;
   }
   
   delete [] Index;
   delete [] CPDOFs;
   delete [] CPLambdas;
   delete [] CPorBifs;
   
   // Reset Lattice and ArcLengthSolution
   ArcLenSet(OriginalDOF);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   // Check to see if we should stop
   if ((StopAtCPCrossingNum_ > -1) && (TotalNumCPCrossings >= StopAtCPCrossingNum_))
      CurrentSolution_ = NumSolutions_;
}

int ArcLengthSolution::ZBrent(Lattice* const Lat,int const& track,Vector const& OriginalDiff,
                              double const& OriginalDS,double& fa,double& fb,Vector& CurrentTF)
{
   ++counter_[1];
   
   int retcode = 1;
   Vector LastDiff(Difference_.Dim(),0.0);
   double LastDS=CurrentDS_;
   double a,b,c,d,e,xm,p,fc, tol1,s,q,r,min1,min2;
   int good = 1;
   int loops = 0;
   double factor = 0.0;
   int oldprecision = cout.precision();
         
   b=OriginalDS;
   c=b;
   a=0.0;
   d=e=0.0; // arbitrary initial values.
   
   fc = fb;
   //cout << " a = " << a << "\n" << "fa = " << fa << "\n"
   //<<  "b = " << b << "\n" << "fb = " << fb << "\n"
   //<< "c = " << c << "\n" << "fc = " << fc << "\n";
   
   while (((fabs(fb) > Tolerance_)) && (loops < MaxIter_))
   {
      cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
      cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";

      ArcLenUpdate(-Difference_);
      
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
      {
         c=a;
         fc=fa;
         e=b-a;
         d=e;
      }
      
      if (fabs(fc) < fabs(fb))
      {
         a=b;
         b=c;
         c=a;
         fa=fb;
         fb=fc;
         fc=fa;
      }
      
      tol1 = 2.0*ARCLENEPS*fabs(b) + 0.5*Tolerance_;
      xm = 0.5 * (c-b);
      
      if(fabs(xm) <= tol1 || fb == 0.0)
      {
         cout << "Minimal Root found! XM too small! " << "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
         fb = CurrentTF[track];

         cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";
         break;
      }
      
      if((fabs(e) >= tol1) && (fabs(fa) > fabs(fb)))
      {
         s=fb/fa;
         if(a == c)
         {
            p=2.0*xm*s;
            q=1.0-s;
         }
         else
         {
            q=fa/fc;
            r=fb/fc;
            p=s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
            q=(q-1.0)*(r-1.0)*(s-1.0);
         }
         
         if(p > 0.0)
         {
            q=-q;
         }
         
         p=fabs(p);
         min1=3.0*xm*q - fabs(tol1*q);
         min2=fabs(e*q);
         
         if(2.0*p < (min1 < min2 ? min1 : min2))
         {
            e=d;
            d=p/q;
         }
         else
         {
            d=xm;
            e=d;
         }
      }
      else
      {
         d=xm;
         e=d;
      }
      
      a=b;
      fa=fb;
      if(fabs(d) > tol1)
      {
         b += d;
      }
      else
      {
         b += (xm >=0.0 ? fabs(tol1) : -fabs(tol1));
      }
      LastDiff = Difference_;
      LastDS = CurrentDS_;
      CurrentDS_ = b;
      factor = OriginalDS/b;
      Difference_ = OriginalDiff/factor;
      ArcLengthNewton(good);
      if (!good)
      {
         // set back to best solution
         cout << "ZBrent stopped due to ArcLengthNewton() convergence failure! " << "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         fb = CurrentTF[track];
         cout << setprecision(30) <<"CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";
         retcode = 0;
         break;
      }
      Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
      
      fb=CurrentTF[track];
      loops++;
      
      //cout << "CurrentTF = " << "\n" << setw(Width) << CurrentTF << "\n" << "\n";
      //cout << "CurrentTF[track] = " << "\n" << CurrentTF[track]<< "\n"<< "\n";
   }
   
   if (loops >= MaxIter_)
   {
      cout << "Error: ZBrent reached Maximum number of iterations before it converged and exited." << "\n";
   }
   
   return retcode;
}

//
//
//// routines to find bif point --- not very robust --- not used anywhere
//void ArcLengthSolution::FindSimpleBif(Lattice* const Lat,Vector const& OriginalDiff,
//                                      double const& OriginalDS,double& fa, double& fb,
//                                      Vector& CurrentTF)
//{
//   //Only works with NoRestriction
//   if (strcmp(Restrict_->Name(),"NoRestriction"))
//   {
//      cerr << "Error: FindSimpleCP only works with NoRestriction object.\n";
//      exit(-1);
//   }
//   
//   // check for turning point
//   int N = Restrict_->DOF().Dim() - 1;
//   int sgn1 = 1;
//   int sgn2 = 1;
//   Matrix const& stiff = Restrict_->Stiffness();
//   Matrix QQ(N+1,N+1),RR(N+1,N);
//   QR(stiff,QQ,RR,1);
//   for (int i=0;i<N;++i)
//   {
//      sgn1 *= int(RR[i][i]/fabs(RR[i][i]));
//   }
//   ArcLenUpdate(-Difference_);
//   Matrix const& stiff2 = Restrict_->Stiffness();
//   QR(stiff2,QQ,RR,1);
//   for (int i=0;i<N;++i)
//   {
//      sgn2 *= int(RR[i][i]/fabs(RR[i][i]));
//   }
//   // set back to RHS of bif pt.
//   ArcLenUpdate(Difference_);
//   // done determining CP type
//   
//   // set back to LHS of bif pt.
//   ArcLenUpdate(-Difference_);
//   // store DOFs for LHS
//   Vector LHSDef = ArcLenDef();
//   
//   Vector currentdef(N+1);
//   // set up augmented equation vector
//   Vector h(2*N+2);
//   Vector oldh(2*N+2);
//   Vector w(2*N+2);
//   Vector dw(2*N+2);
//   
//   // create storage for force and stiffness and eigenvectors and values
//   Vector force(N);
//   Matrix eigvecs(N,N);
//   Matrix eigvals(1,N);
//   
//   // set to initial guess for bif point
//   Vector initialguess = (-fa/(fb-fa))*OriginalDiff;
//   eigvals=SymEigVal(Lat->E2(),&eigvecs);
//   int minindex = 0;
//   double minval = fabs(eigvals[0][0]);
//   for (int i=1;i<N;++i)
//   {
//      if (fabs(eigvals[0][i]) < minval)
//      {
//         minindex = i;
//         minval = fabs(eigvals[0][i]);
//      }
//   }
//   
//   // fill w = (x,z,alpha)
//   for (int i=0;i<N+1;++i)
//   {
//      // x
//      w[i] = LHSDef[i] + initialguess[i];
//   }
//   for (int i=0;i<N;++i)
//   {
//      // z
//      w[N+1+i] = eigvecs[i][minindex];
//   }
//   // alpha
//   w[2*N+1] = 0.0;
//   
//   // update dofs
//   for (int i=0;i<N+1;++i)
//   {
//      currentdef[i]=w[i];
//   }
//   ArcLenSet(currentdef);
//   
//   
//   // update h
//   SetH(h,w);
//   // set oldh
//   oldh = h;
//   
//   cout << "h.Norm()=" << setw(20) << h.Norm() << "\n";
//   
//   // Only perform calculation if CP is a bif pt
//   if (sgn1 == -sgn2) // bif pt
//   {
//      // start Newton loops.
//      // setup initial guess for jacobian
//      Matrix J(2*N+2,2*N+2,0.0);
//      // identity as guess for D2f^T*z
//      Matrix const& e3 = Lat->E3();
//      for (int i=0;i<N;++i)
//      {
//         for (int j=0;j<N;++j)
//         {
//            J[i][j] = 0.0;
//            for (int k=0;k<N;++k)
//            {
//               J[i][j] += e3[N*i+j][k]*w[N+1+k];
//            }
//         }
//      }
//      Matrix const& stiffdl = (Lat->LoadParameter() == Lattice::Load) ?
//         Lat->StiffnessDL() : Lat->StiffnessDT() ;
//      for (int i=0;i<N;++i)
//      {
//         J[N][i] = J[i][N] = 0.0;
//         for (int j=0;j<N;++j)
//         {
//            J[i][N] = J[N][i] += stiffdl[i][j]*w[N+1+j];
//         }
//      }
//      // Df terms
//      for (int i=0;i<N+1;++i)
//      {
//         for (int j=0;j<N;++j)
//         {
//            J[N+1+j][i] = stiff[j][i];
//            J[i][N+1+j] = stiff[j][i];
//         }
//      }
//      // zero terms
//      for (int i=0;i<N+1;++i)
//      {
//         J[i][2*N+1] = 0.0;
//         J[2*N+1][i] = 0.0;
//      }
//      J[2*N+1][2*N+1] = 0.0;
//      //center -alpha*identity terms
//      for (int i=N+1;i<2*N+1;++i)
//      {
//         J[i][i] = -w[2*N+1];
//      }
//      // z terms
//      for (int i=0;i<N;++i)
//      {
//         J[N+1+i][2*N+1] = -w[N+1+i];
//         J[2*N+1][N+1+i] = w[N+1+i]/2.0;
//      }
//      
//      Matrix Q,R;
//      Q.SetIdentity(2*N+2);
//      R.SetIdentity(2*N+2);
//      QR(J,Q,R);
//      
//      int loops = 0;
//      do
//      {
//         SolveQR(Q,R,dw,h);
//         w -= dw;
//         // update dofs
//         for (int i=0;i<N+1;++i)
//         {
//            currentdef[i]=w[i];
//         }
//         ArcLenSet(currentdef);
//         
//         // calculate h
//         SetH(h,w);
//         // done calculating h
//         cout << "dw.Norm()=" << setw(20) << dw.Norm()
//              << ",  h.Norm()=" << setw(20) << h.Norm() << "\n";
//         
//         // Update Q and R
//         BroydenQRUpdate(Q,R,(h-oldh)/dw.Norm(),-dw/dw.Norm());
//         
//         // update oldh
//         oldh = h;
//         ++loops;
//      }
//      while ((h.Norm() > Tolerance_) && (loops < MaxIter_));
//      
//      if (loops >= MaxIter_)
//      {
//         cout << "Error: FindSimpleBif() did not converge\n";
//      }
//   }
//   else
//   {
//      cout << "Warning: FindSimpleCP identified a turningpoint CP\n"
//           << "         and has not accurately polished the point.\n";
//   }
//   
//   // set Difference_ and CurrentDS_ appropriately;
//   Difference_ = currentdef - LHSDef;
//   CurrentDS_ = 0.0;
//   for (int i=0;i<N;++i)
//   {
//      CurrentDS_ += Difference_[i]*Difference_[i];
//   }
//   CurrentDS_ += Difference_[N]*Difference_[N]/(Aspect_*Aspect_);
//   
//   // Find current testfunctions
//   Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
//}
//
//void ArcLengthSolution::SetH(Vector& h,Vector const& w)
//{
//   // h(x,z,alpha) = (Df(x)^T*z,f(x)-alpha*z,(z*z - 1)/2)
//   // size           (N+1,      N,            1) = 2N+2
//   // w = (x,z,alpha)
//   //size (N+1,N,1) = 2N+2
//   
//   Matrix const& stiff = Restrict_->Stiffness();
//   int N=stiff.Rows();
//   
//   // update Df(x)^T*z
//   for (int i=0;i<N+1;++i)
//   {
//      h[i] = 0.0;
//      for (int j=0;j<N;++j)
//      {
//         h[i] += stiff[j][i]*w[N+1+j];
//      }
//   }
//   
//   // update f(x)-alpha*z
//   Vector const& stress = Restrict_->Force();
//   for (int i=0;i<N;++i)
//   {
//      h[N+1+i] = stress[i] - w[2*N+1]*w[N+1+i];
//   }
//   
//   // update (z*z - 1)/2
//   h[2*N+1] = 0.0;
//   for (int i=0;i<N;++i)
//   {
//      h[2*N+1] += w[N+1+i]*w[N+1+i];
//   }
//   h[2*N+1] -= 1.0;
//   h[2*N+1] /= 2.0;
//}
//
//
