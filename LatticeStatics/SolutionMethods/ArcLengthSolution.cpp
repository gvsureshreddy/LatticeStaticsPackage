#include <math.h>
#include <sstream>
#include <string>
#include "ArcLengthSolution.h"

using namespace std;

#define ARCLENEPS 1.0e-15

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,Vector const& dofs,
                                     int const& MaxIter,double const& Tolerance,
                                     double const& BisectTolerance,double const& DSMax,
                                     double const& DSMin,double const& CurrentDS,
                                     double const& AngleCutoff,double const& AngleIncrease,
                                     double const& Aspect,int const& NumSolutions,
                                     int const& CurrentSolution,Vector const& FirstSolution,
                                     Vector const& Difference,int const& ClosedLoopStart,
                                     int const& StopAtCPNum,int const& Echo)
   : Echo_(Echo),
     Restrict_(Restrict),
     DOFS_(Restrict_->DOF().Dim()),
     MaxIter_(MaxIter),
     Tolerance_(Tolerance),
     BisectTolerance_(BisectTolerance),
     DSMax_(DSMax),
     DSMin_(DSMin),
     CurrentDS_(CurrentDS),
     AngleCutoff_(AngleCutoff),
     AngleIncrease_(AngleIncrease),
     Aspect_(Aspect),
     NumSolutions_(NumSolutions_),
     CurrentSolution_(CurrentSolution),
     ClosedLoopStart_(ClosedLoopStart),
     FirstSolution_(FirstSolution),
     StopAtCPNum_(StopAtCPNum),
     Difference_(Difference),
     force_static(DOFS_),
     mdfc_static(DOFS_-1),
     K_static(DOFS_,DOFS_),
     RestrictK_static(DOFS_-1,DOFS_)
{
   ArcLenSet(dofs);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,
                                     Vector const& one,Vector const& two,int const& Echo)
   : Echo_(Echo),
     Restrict_(Restrict),
     CurrentSolution_(0),
     Difference_(two-one)
{
   DOFS_=Restrict_->DOF().Dim();
   // initialize "static" members variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_-1);
   K_static.Resize(DOFS_,DOFS_);
   RestrictK_static.Resize(DOFS_-1,DOFS_);

   
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash,"MaxIterations");
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
   if (Input.ParameterOK(Hash,"StopAtCPNum"))
   {
      StopAtCPNum_ = Input.getInt(Hash,"StopAtCPNum");
   }
   else
   {
      StopAtCPNum_ = Input.useInt(-1,Hash,"StopAtCPNum"); // Default Value
   }
   Input.EndofInputSection();
   
   BisectTolerance_ = Tolerance_;
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   
   // Set Lattice to solution "two"
   ArcLenSet(two);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict,PerlInput const& Input,
                                     int const Echo)
   :  Echo_(Echo),
      Restrict_(Restrict),
      CurrentSolution_(0)
{
   DOFS_=Restrict_->DOF().Dim();
   // initialize "static" memver variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_-1);
   K_static.Resize(DOFS_,DOFS_);
   RestrictK_static.Resize(DOFS_-1,DOFS_);

   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash,"MaxIterations");
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
   if (Input.ParameterOK(Hash,"StopAtCPNum"))
   {
      StopAtCPNum_ = Input.getInt(Hash,"StopAtCPNum");
   }
   else
   {
      StopAtCPNum_ = Input.useInt(-1,Hash,"StopAtCPNum"); // Default Value
   }
   Input.EndofInputSection();
   
   BisectTolerance_ = Tolerance_;
   
   const char *starttype = Input.getString("StartType","Type");

   if (!strcmp("Bifurcation",starttype))
   {
      // Set Difference and Lattice state
      double eps = Input.getDouble("StartType","Epsilon");

      Difference_.Resize(DOFS_);
      Vector diff(Input.getArrayLength("StartType","Tangent"));
      Input.getVector(diff,"StartType","Tangent");
      Difference_ = Restrict_->TransformVector(diff);
      Difference_ *= eps;
      
      Vector stat(Input.getArrayLength("StartType","BifurcationPoint"));
      Input.getVector(stat,"StartType","BifurcationPoint");
      // Set Lattice state to the bifurcation point
      ArcLenSet(Restrict_->RestrictDOF(stat));
      
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
         FirstSolution_ = ArcLenDef();
         Input.useVector(FirstSolution_,"StartType","ClosedLoopFirstSolution"); // Default Value
      }
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
         FirstSolution_ = ArcLenDef();
         Input.useVector(FirstSolution_,"StartType","ClosedLoopFirstSolution"); // Default Value
      }
   }
   else if (!strcmp("ConsistenceCheck",starttype))
   {
      double ConsistencyEpsilon;
      int Width;
      Vector Solution1(DOFS_);
      Vector Solution2(DOFS_);

      Vector onetmp(Input.getArrayLength("StartType","Solution1"));
      Input.getVector(onetmp,"StartType","Solution1");
      Solution1 = Restrict_->RestrictDOF(onetmp);
      Vector twotmp(Input.getArrayLength("StartType","Solution2"));
      Input.getVector(twotmp,"StartType","Solution2");
      Solution2 = Restrict_->RestrictDOF(twotmp);
      // Get Epsilon and Width
      ConsistencyEpsilon = Input.getDouble("StartType","ConsistenceEpsilon");
      Width = Input.getPosInt("Main","FieldWidth");

      cout << "FIX UP CONSISTENCYCHECK!!!!!!\n";
      //ConsistencyCheck(Solution1,Solution2,ConsistencyEpsilon,Width,out);
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

void ArcLengthSolution::ConsistencyCheck(Vector const& Solution1,Vector const& Solution2,
                                         double const& ConsistencyEpsilon,int const& Width,
                                         fstream& out)
{
   double potential;
   int Dim=ArcLenDef().Dim();
   Matrix
      Stiff(Dim,Dim),
      PerturbedStiff(Dim,Dim);
   Vector
      Force(Dim),
      PerturbedForce(Dim),
      pert(Dim,0.0),
      RHS(Dim);
   
   Difference_.Resize(Dim);
   
   // Set Lattice state to Solution2
   ArcLenSet(Solution2);
   // Set Difference to Solution2 - Solution1
   Difference_ = Solution2 - Solution1;
   
   // Do Consistency check
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << "\n";
      cout << "Consistency Check." << "\n";
      cout << "F(U + DeltaU) * Epsilon" << "\n";
   }
   for (int i=0;i<70;i++) out << "="; out << "\n";
   out << "Consistency Check." << "\n";
   out << "F(U + DeltaU) * Epsilon" << "\n";
   ArcLenUpdate(Difference_);
   Force = ConsistencyEpsilon*ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
   if (Echo_)
   {
      cout << setw(Width) << Force << "\n";
      cout << "K(U + DeltaU) * Epsilon" << "\n";
   }
   out << setw(Width) << Force << "\n";
   out << "K(U + DeltaU) * Epsilon" << "\n";
   Stiff = ConsistencyEpsilon*ArcLenStiffness(Difference_,1.0);
   if (Echo_) cout << setw(Width) << Stiff << "\n";
   out << setw(Width) << Stiff << "\n";
   for (int i=0;i<Difference_.Dim();i++)
   {
      // Get RHS
      Difference_ = Solution2 - Solution1;
      ArcLenSet(Solution2 + Difference_);
      potential = Restrict_->Energy();
      RHS = ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
      
      // Perturb the lattice state
      pert=Vector(pert.Dim(),0.0);
      pert[i]=1.0;
      Difference_ = Solution2 - Solution1 + ConsistencyEpsilon*pert;
      ArcLenSet(Solution2 + Difference_);
      // Get Check
      potential = Restrict_->Energy() - potential;
      // fix-up the arclength equation part of PerturbedForce
      if (i == RHS.Dim()-1) potential = ConsistencyEpsilon*RHS[i];
      PerturbedForce[i] = potential;
      RHS = ArcLenForce(ConsistencyEpsilon,Difference_,1.0) - RHS;
      for (int j=0;j<Dim;j++)
         PerturbedStiff[j][i] = RHS[j];
   }
   
   // Print out the facts
   if (Echo_)
   {
      cout << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << "\n";
      cout << setw(Width) << PerturbedForce << "\n";
      cout << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << "\n";
      cout << setw(Width) << PerturbedStiff << "\n";
      cout << "Difference" << "\n";
      cout << setw(Width) << Force - PerturbedForce << "\n" << "\n";
      cout << setw(Width) << Stiff - PerturbedStiff << "\n";
   }
   
   out << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedForce << "\n";
   out << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedStiff << "\n";
   out << "Difference" << "\n";
   out << setw(Width) << Force - PerturbedForce << "\n" << "\n";
   out << setw(Width) << Stiff - PerturbedStiff << "\n";
   
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << "\n";
   }
   for (int i=0;i<70;i++) out << "="; out << "\n";
   
   // We are done -- set currentsolution to numsolutions
   CurrentSolution_ = NumSolutions_;
}

int ArcLengthSolution::AllSolutionsFound() const
{
   return (CurrentSolution_ >= NumSolutions_);
}

int ArcLengthSolution::FindNextSolution()
{
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
      
      if (Echo_)
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
      cerr << "Closed Loop detected at Solution # " << CurrentSolution_
           << " --- Terminating!" << "\n";
      
      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      CurrentSolution_++;
   }
   
   // Always have the current "solution" state printed as a solution point
   good = 1;
   
   return good;
}

void ArcLengthSolution::ArcLengthNewton(int& good)
{
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
   if (Echo_) cout << "Taking Predictor Step. CurrentDS = " << CurrentDS_ << "\n";
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

   if (Echo_) cout << "\tForceNorm = " << ForceNorm << " \tDeltaNorm = " << Magnitude2 << "\n";

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
      
      if (Echo_) cout << "\tForceNorm = " << ForceNorm
                      << " \tDeltaNorm = " << Magnitude2
                      << " \tContraction = " << Kappa << "\n";
   }
   while ((itr < MaxIter_) && ((RHS.Norm() > Tolerance_) || (Dx.Norm() > Tolerance_)));
   
   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << "\n";
      good = 0;
   }
   else
   {
      if (Echo_) cout << "Prediction 1 Corrector Iterations: " << itr << "\n";
      good = 1;
   }
}

void ArcLengthSolution::OldFindCriticalPoint(int const& LHN,double const& LHEV,int const& RHN,
                                             double const& RHEV,Lattice* const Lat,
                                             int& TotalNumCPs,PerlInput const& Input,
                                             int const& Width,fstream& out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   double CurrentMinEV=1.0, OldMinEV=1.0;
   double Delta_DS=0.0;
   int dummy = 1;
   int loops = 0;
   int RighthandTestValue = RHN;
   int CurrentTestValue = RHN;
   int LefthandTestValue = LHN;
   Vector EigenValues(Lat->DOF().Dim());
   
   Delta_DS = CurrentDS_;
   
   Vector Original_DOF((Restrict_->DOF()).Dim(),0);
   Original_DOF = Restrict_->DOF();
   
   if (Echo_) cout << "LHN = " << LHN << "\n" << "RHN = " << RHN << "\n";
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > BisectTolerance_)
           || (CurrentTestValue == RighthandTestValue))
          && (loops < MaxIter_))
   {
      if (Echo_) cout << "OldMinEV = " << OldMinEV << "\n"
                      << "CurrentTestValue = " << CurrentTestValue << "\n"
                      << "CurrentMinEV = "  << CurrentMinEV << "\n";
      
      //set to left hand point
      ArcLenUpdate(-Difference_);
      Delta_DS = Delta_DS/2.0;
      
      if(CurrentTestValue == RighthandTestValue)
      {
         CurrentDS_ = CurrentDS_ - Delta_DS;
         Difference_ = Difference_/2.0;
      }
      if(CurrentTestValue == LefthandTestValue)
      {
         CurrentDS_ = CurrentDS_ + Delta_DS;
         Difference_ = 1.5 * Difference_;
      }
      
      cout << "Current_DS = " << CurrentDS_ << "\n" << "\n";
      
      ArcLengthNewton(dummy);
      
      OldMinEV = CurrentMinEV;
      CurrentTestValue = Lat->TestFunctions(EigenValues);
      CurrentMinEV = EigenValues[1]; //not correct;
      
      loops++;
   }
   
   // Output Critical Point
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << "\n";
   out << "\n";

   // Lattice takes care of echo
   out << setw(Width) << *Lat << "\n";
      
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << "\n"; out << "\n";
   
   // Call Lattice function to do any Lattice Specific things
   //  abs(RighthandNulity - LefthandNulity) is the number of zero eigenvalues
   //  in a perfect situation. should check to see if this is found to be true.
   Lat->CriticalPointInfo(Restrict_->DrDt(Difference_),
                          abs(RighthandTestValue-LefthandTestValue),
                          BisectTolerance_,Width,Input,out,cout);
   
   if (Echo_) cout << "Success = 1" << "\n";
   out << "Success = 1" << "\n";
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   // Check to see if we should stop
   TotalNumCPs += 1;
   if ((StopAtCPNum_ > -1) && (TotalNumCPs >= StopAtCPNum_))
      CurrentSolution_ = NumSolutions_;
}

void ArcLengthSolution::FindCriticalPoint(Lattice* const Lat,int& TotalNumCPs,
                                          PerlInput const& Input,int const& Width,fstream& out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   int TestValueDiff;
   int temp;
   int size = Lat->DOF().Dim();
   TF_LHS_static.Resize(size);
   TF_RHS_static.Resize(size);
   CurrentTF_static.Resize(size);
   double fa,fb;
   int Multiplicity;
   int track;
   int num;
   int CP;
   int spot;
   ostringstream in_string;
   ostringstream in_newinput_string;

   // Setup in_string ios
   in_string << setiosflags(ios::fixed) << setprecision(out.precision());
   
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
   
   cout << "TF_LHS_static = " << setw(Width) << TF_LHS_static<< "\n";
   cout << "TF_RHS_static = " << setw(Width) << TF_RHS_static << "\n";
   
   int* Index;
   Index = new int[TestValueDiff];
   Vector DSTrack(TestValueDiff);
   string* out_string = new string[TestValueDiff];
   string* out_newinput_string = new string[TestValueDiff];
   int* Bif = new int[TestValueDiff];

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
   if (Echo_)
   {
      for (int i = 0; i < TestValueDiff; i++)
      {
         cout << "i = " << i<< "\n";
         cout << "Index[i] = " << Index[i]<< "\n";
      }
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
         Multiplicity = 1;
         for(int i=CP+1;i<TestValueDiff;i++)
         {
            temp = Index[i];
            if(fabs(CurrentTF_static[temp]) <= BisectTolerance_)
            {
               Index[i] = -1;
               Multiplicity++;
               //cout <<"i = " << i << "\n" <<  "CHECK POINT INDEX[CP] = " << Index[i] << "\n";
            }
         }
         
         // sort the critical points
         spot=num;
         while((spot!=0)&&(DSTrack[spot-1] > CurrentDS_) )
         {
            DSTrack[spot] = DSTrack[spot-1];
            out_string[spot] = out_string[spot - 1];
            out_newinput_string[spot] = out_newinput_string[spot - 1];
            Bif[spot] = Bif[spot - 1];
            spot = spot - 1;
         }         
         
         // Output Critical Point
         for (int i=0;i<70;i++)
         {
            if (Echo_) cout << "=";
            in_string << "=";
         }
         if (Echo_) cout << "\n";
         in_string << "\n";
         
         // Lattice takes care of echo
         in_string << setw(Width) << *Lat << "\n";         
         
         for (int i=0;i<70;i++)
         {
            if (Echo_) cout << "=";
            in_string << "=";
         }
         if (Echo_) cout << "\n";
         in_string << "\n";

         // Call Lattice function to do any Lattice Specific things
         Bif[spot]=Lat->CriticalPointInfo(Restrict_->DrDt(Difference_),Multiplicity,
                                          10.0*Tolerance_,Width,Input,in_string,
                                          in_newinput_string);
         if (Echo_) cout << "Success = 1" << "\n";
         in_string << "Success = 1" << "\n";
         
         DSTrack[spot] = CurrentDS_;
         out_string[spot] = in_string.str();
         out_newinput_string[spot] = in_newinput_string.str();
         num = num + 1;
         in_string.str("");
         in_newinput_string.str("");
      }//END OF IF STATEMENT
   }
   
   ////PRINT OUT CP DATA
   fstream cpfile;
   
   for (int i = 0; i < num; i++)
   {
      out << out_string[i];
      if (2 == Bif[i])
      {
         in_string.str("");
         in_string << Input.LastInputFileName() << ".CP."
                   << setw(2) << setfill('0') << TotalNumCPs;
         cpfile.open(in_string.str().c_str(),ios::out);
         cpfile << out_newinput_string[i];
         cpfile.close();
      }
      else if (1 == Bif[i])
      {
         in_string.str("");
         in_string << Input.LastInputFileName() << ".BP."
                   << setw(2) << setfill('0') << TotalNumCPs;
         cpfile.open(in_string.str().c_str(),ios::out);
         cpfile << out_newinput_string[i];
         cpfile.close();
      }
      else
      {
         in_string.str("");
         in_string << Input.LastInputFileName() << ".TP."
                   << setw(2) << setfill('0') << TotalNumCPs;
         cpfile.open(in_string.str().c_str(),ios::out);
         cpfile << out_newinput_string[i];
         cpfile.close();
      }
         ++TotalNumCPs;
   }
   
   delete [] Index;
   delete [] out_string;   //deletes memory allocated to out_string
   delete [] out_newinput_string;
   delete [] Bif;
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   // Check to see if we should stop
   if ((StopAtCPNum_ > -1) && (TotalNumCPs >= StopAtCPNum_))
      CurrentSolution_ = NumSolutions_;
}

void ArcLengthSolution::ZBrent(Lattice* const Lat,int const& track,Vector const& OriginalDiff,
                               double const& OriginalDS,double& fa,double& fb,Vector& CurrentTF)
{
   Vector LastDiff(Difference_.Dim(),0.0);
   double LastDS=CurrentDS_;
   double a,b,c,d,e,xm,p,fc, tol1,s,q,r,min1,min2;
   int dummy = 1;
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
   
   while (((fabs(fb) > BisectTolerance_)) && (loops < MaxIter_))
   {
      if (Echo_)
      {
         cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";
      }
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
      
      tol1 = 2.0*ARCLENEPS*fabs(b) + 0.5*BisectTolerance_;
      xm = 0.5 * (c-b);
      
      if(fabs(xm) <= tol1 || fb == 0.0)
      {
         cout << "Minimal Root found! XM too small! " << "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
         fb=CurrentTF[track];
         if(Echo_)
         {
            cout <<setprecision(30)<<"CurrentMinTF = " << fb << "\n";
            cout << "CurrentDS_ = " << CurrentDS_ <<setprecision(oldprecision)<< "\n";
         }
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
      ArcLengthNewton(dummy);
      Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
      
      fb=CurrentTF[track];
      loops++;
      
      //cout << "CurrentTF = " << "\n" << setw(Width) << CurrentTF << "\n" << "\n";
      //cout << "CurrentTF[track] = " << "\n" << CurrentTF[track]<< "\n"<< "\n";
   }
}
