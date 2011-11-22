#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include "ArcLengthSolution.h"

using namespace std;

#define ARCLENEPS 1.0e-15

extern "C" void qcbfb_output_(int& nfree, double* u, double& prop, int& nint, int* intdata, int& ndouble, double* doubledata);

ArcLengthSolution::~ArcLengthSolution()
{
   cout.width(0);
   cout << "ArcLengthSolution Stats:\n"
        << "\tCumulativeArcLength - " << CumulativeArcLength_ << "\n";
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

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict, Vector const& dofs,
                                     int const& MaxIter, double const& Tolerance,
                                     ConvergeType CnvrgTyp, double const& DSMax,
                                     double const& DSMin, double const& CurrentDS,
                                     double const& AngleCutoff, double const& AngleIncrease,
                                     double const& Aspect, double const& eig_angle_max,
                                     int const& NumSolutions, int const& CurrentSolution,
                                     Vector const& FirstSolution, Vector const& Difference,
                                     int const& BifStartFlag, Vector const& BifTangent,
                                     int const& ClosedLoopStart, int const& ClosedLoopUseAsFirst,
                                     double const& MaxCumulativeArcLength,
                                     int const& StopAtCPCrossingNum, int const& Echo) :
   Echo_(Echo),
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
   eig_angle_max_(eig_angle_max),
   NumSolutions_(NumSolutions_),
   CumulativeArcLength_(0.0),
   CurrentSolution_(CurrentSolution),
   BifStartFlag_(BifStartFlag),
   BifTangent_(BifTangent),
   ClosedLoopStart_(ClosedLoopStart),
   ClosedLoopUseAsFirst_(ClosedLoopUseAsFirst),
   FirstSolution_(FirstSolution),
   MaxCumulativeArcLength_(MaxCumulativeArcLength),
   StopAtCPCrossingNum_(StopAtCPCrossingNum),
   Difference_(Difference),
   force_static(DOFS_),
   mdfc_static(DOFS_ - 1),
   K_static(DOFS_, DOFS_),
   RestrictK_static(DOFS_ - 1, DOFS_)
{
   for (int i = 0; i < nocounters_; ++i)
   {
      counter_[i] = 0;
   }
   ArcLenSet(dofs);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict, PerlInput const& Input,
                                     Vector const& one, Vector const& two, int const& Echo) :
   Echo_(Echo),
   Restrict_(Restrict),
   CumulativeArcLength_(0.0),
   CurrentSolution_(0),
   BifStartFlag_(0),
   BifTangent_(),
   Difference_(two - one)
{
   for (int i = 0; i < nocounters_; ++i)
   {
      counter_[i] = 0;
   }

   DOFS_ = Restrict_->DOF().Dim();
   // initialize "static" members variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_ - 1);
   K_static.Resize(DOFS_, DOFS_);
   RestrictK_static.Resize(DOFS_ - 1, DOFS_);


   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash, "MaxIterations");
   if (Input.ParameterOK(Hash, "ConvergeType"))
   {
      char const* const cnvrgtyp = Input.getString(Hash, "ConvergeType");
      if (!strcmp("Both", cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force", cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement", cnvrgtyp))
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
      Input.useString("Both", Hash, "ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
   Tolerance_ = Input.getDouble(Hash, "Tolerance");
   DSMax_ = Input.getDouble(Hash, "DSMax");
   CurrentDS_ = Input.getDouble(Hash, "DSStart");
   DSMin_ = Input.getDouble(Hash, "DSMin");
   AngleCutoff_ = Input.getDouble(Hash, "AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash, "AngleIncrease");
   Aspect_ = Input.getDouble(Hash, "Aspect");
   if (Input.ParameterOK(Hash, "MaxEigVectAngle"))
   {
      eig_angle_max_ = Input.getDouble(Hash, "MaxEigVectAngle");
   }
   else
   {
      Input.useDouble(-1.0, Hash, "MaxEigVectAngle"); // Default Value
   }
   NumSolutions_ = Input.getPosInt(Hash, "NumSolutions");
   if (Input.ParameterOK(Hash, "ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash, "ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT, Hash, "ClosedLoopStart"); // Default Value
   }
   if (Input.ParameterOK(Hash, "MaxCumulativeArcLength"))
   {
      MaxCumulativeArcLength_ = Input.getDouble(Hash, "MaxCumulativeArcLength");
      if (MaxCumulativeArcLength_ < 0.0)
      {
         MaxCumulativeArcLength_ = -1.0;  // Negative values mean don't check
      }
   }
   else
   {
      MaxCumulativeArcLength_ = Input.useDouble(-1.0, Hash, "MaxCumulativeArcLength_"); // Default Value
   }
   if (Input.ParameterOK(Hash, "StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash, "StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1, Hash, "StopAtCPCrossingNum"); // Default Value
   }


   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   if (Input.ParameterOK("StartType", "ClosedLoopFirstSolution"))
   {
      Vector clfs(Input.getArrayLength("StartType", "ClosedLoopFirstSolution"));
      Input.getVector(clfs, "StartType", "ClosedLoopFirstSolution");
      FirstSolution_ = Restrict_->RestrictDOF(clfs);
   }
   else
   {
      if (Input.ParameterOK("StartType", "ClosedLoopUseAsFirst"))
      {
         ClosedLoopUseAsFirst_ = Input.getPosInt("StartType", "ClosedLoopUseAsFirst");
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
         ClosedLoopUseAsFirst_ = Input.usePosInt(0, "StartType", "ClosedLoopUseAsFirst");
         FirstSolution_ = ArcLenDef();
      }
   }
   Input.EndofInputSection();

   // Set Lattice to solution "two"
   ArcLenSet(two);
}

ArcLengthSolution::ArcLengthSolution(Restriction* const Restrict, PerlInput const& Input,
                                     int const Echo) :
   Echo_(Echo),
   Restrict_(Restrict),
   CumulativeArcLength_(0.0),
   CurrentSolution_(0),
   BifStartFlag_(0),
   BifTangent_()
{
   for (int i = 0; i < nocounters_; ++i)
   {
      counter_[i] = 0;
   }

   DOFS_ = Restrict_->DOF().Dim();
   // initialize "static" memver variables
   force_static.Resize(DOFS_);
   mdfc_static.Resize(DOFS_ - 1);
   K_static.Resize(DOFS_, DOFS_);
   RestrictK_static.Resize(DOFS_ - 1, DOFS_);

   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "ArcLengthSolution");
   MaxIter_ = Input.getPosInt(Hash, "MaxIterations");
   if (Input.ParameterOK(Hash, "ConvergeType"))
   {
      char const* const cnvrgtyp = Input.getString(Hash, "ConvergeType");
      if (!strcmp("Both", cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force", cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement", cnvrgtyp))
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
      Input.useString("Both", Hash, "ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
   Tolerance_ = Input.getDouble(Hash, "Tolerance");
   DSMax_ = Input.getDouble(Hash, "DSMax");
   CurrentDS_ = Input.getDouble(Hash, "DSStart");
   DSMin_ = Input.getDouble(Hash, "DSMin");
   AngleCutoff_ = Input.getDouble(Hash, "AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash, "AngleIncrease");
   Aspect_ = Input.getDouble(Hash, "Aspect");
   if (Input.ParameterOK(Hash, "MaxEigVectAngle"))
   {
      eig_angle_max_ = Input.getDouble(Hash, "MaxEigVectAngle");
   }
   else
   {
      Input.useDouble(-1.0, Hash, "MaxEigVectAngle"); // Default Value
   }
   NumSolutions_ = Input.getPosInt(Hash, "NumSolutions");
   if (Input.ParameterOK(Hash, "ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash, "ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT, Hash, "ClosedLoopStart"); // Default Value
   }
   if (Input.ParameterOK(Hash, "MaxCumulativeArcLength"))
   {
      MaxCumulativeArcLength_ = Input.getDouble(Hash, "MaxCumulativeArcLength");
      if (MaxCumulativeArcLength_ < 0.0)
      {
         MaxCumulativeArcLength_ = -1.0;  // Negative values mean don't check
      }
   }
   else
   {
      MaxCumulativeArcLength_ = Input.useDouble(-1.0, Hash, "MaxCumulativeArcLength_"); // Default Value
   }
   if (Input.ParameterOK(Hash, "StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash, "StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1, Hash, "StopAtCPCrossingNum"); // Default Value
   }
   Input.EndofInputSection();

   const char* starttype = Input.getString("StartType", "Type");

   if (!strcmp("Bifurcation", starttype))
   {
      // Bifurcation
      BifStartFlag_ = 1;
      BifTangent_.Resize(DOFS_);

      // Set Difference and Lattice state
      double eps = Input.getDouble("StartType", "Epsilon");

      Difference_.Resize(DOFS_);
      Vector diff(Input.getArrayLength("StartType", "Tangent"));
      Input.getVector(diff, "StartType", "Tangent");
      Difference_ = Restrict_->TransformVector(diff);
      Difference_ *= eps;
      // set BifTangent_ for projection output
      BifTangent_ = Difference_ / Difference_.Norm();
      BifTangent_[DOFS_ - 1] = 0.0;

      Vector stat(Input.getArrayLength("StartType", "BifurcationPoint"));
      Input.getVector(stat, "StartType", "BifurcationPoint");
      // Set Lattice state to the bifurcation point
      Vector RestrictedStat = Restrict_->RestrictDOF(stat);
      ArcLenSet(RestrictedStat);

      // Set FirstSolution
      FirstSolution_.Resize(DOFS_);
      if (Input.ParameterOK("StartType", "ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType", "ClosedLoopFirstSolution"));
         Input.getVector(clfs, "StartType", "ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType", "ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType", "ClosedLoopUseAsFirst");
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
            ClosedLoopUseAsFirst_ = Input.useInt(0, "StartType", "ClosedLoopUseAsFirst");
            FirstSolution_ = ArcLenDef();
         }
      }

      cout << "Projection on BifTangent of BifurcationPoint = " << RestrictedStat * BifTangent_ << "\n";
   }
   else if (!strcmp("Continuation", starttype))
   {
      // Get solution1
      Vector one(DOFS_);
      Vector onetmp(Input.getArrayLength("StartType", "Solution1"));
      Input.getVector(onetmp, "StartType", "Solution1");
      one = Restrict_->RestrictDOF(onetmp);

      // Set Lattice state to Solution2
      Vector two(DOFS_);
      Vector twotmp(Input.getArrayLength("StartType", "Solution2"));
      Input.getVector(twotmp, "StartType", "Solution2");
      two = Restrict_->RestrictDOF(twotmp);
      ArcLenSet(two);

      // Set Difference_ to   two - one
      Difference_.Resize(DOFS_);
      Difference_ = two - one;

      // Set FirstSolution
      FirstSolution_.Resize(DOFS_);
      if (Input.ParameterOK("StartType", "ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType", "ClosedLoopFirstSolution"));
         Input.getVector(clfs, "StartType", "ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType", "ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType", "ClosedLoopUseAsFirst");
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
            ClosedLoopUseAsFirst_ = Input.useInt(0, "StartType", "ClosedLoopUseAsFirst");
            FirstSolution_ = ArcLenDef();
         }
      }
   }
   else if (!strcmp("ConsistencyCheck", starttype))
   {
      double ConsistencyEpsilon;
      int Width;
      Vector Solution(DOFS_);

      Vector onetmp(Input.getArrayLength("StartType", "Solution"));
      Input.getVector(onetmp, "StartType", "Solution");
      Solution = Restrict_->RestrictDOF(onetmp);
      // Get Epsilon and Width
      ConsistencyEpsilon = Input.getDouble("StartType", "Epsilon");
      Width = Input.getPosInt("Main", "FieldWidth");

      ostream::fmtflags oldflags = cout.flags();
      cout << scientific;
      Restrict_->ConsistencyCheck(Solution, ConsistencyEpsilon, Width, cout);
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

Vector const& ArcLengthSolution::ArcLenForce(double const& DS, Vector const& Diff,
                                             double const& Aspect) const
{
   ++counter_[2];

   mdfc_static = Restrict_->Force();

   force_static[DOFS_ - 1] = DS * DS - Diff[DOFS_ - 1] * Diff[DOFS_ - 1] / (Aspect * Aspect);
   for (int i = 0; i < DOFS_ - 1; ++i)
   {
      force_static[i] = mdfc_static[i];
      force_static[DOFS_ - 1] -= Diff[i] * Diff[i];
   }

   return force_static;
}

Matrix const& ArcLengthSolution::ArcLenStiffness(Vector const& Diff, double const& Aspect)
const
{
   ++counter_[7];

   RestrictK_static = Restrict_->Stiffness();

   for (int i = 0; i < DOFS_ - 1; ++i)
   {
      for (int j = 0; j <= DOFS_ - 1; ++j)
      {
         K_static[i][j] = RestrictK_static[i][j];
      }
      K_static[DOFS_ - 1][i] = -2.0 * Diff[i];
   }
   K_static[DOFS_ - 1][DOFS_ - 1] = -2.0 * Diff[DOFS_ - 1] / (Aspect * Aspect);

   return K_static;
}

double ArcLengthSolution::ArcLenAngle(Vector const& Old, Vector const& New, double const& Aspect)
const
{
   ++counter_[6];

   double angle = 0.0;
   double NewNorm = 0.0;
   double OldNorm = 0.0;

   for (int i = 0; i < DOFS_ - 1; ++i)
   {
      angle += Old[i] * New[i];
      NewNorm += New[i] * New[i];
      OldNorm += Old[i] * Old[i];
   }

   angle += Old[DOFS_ - 1] * New[DOFS_ - 1] / (Aspect * Aspect);
   NewNorm += New[DOFS_ - 1] * New[DOFS_ - 1] / (Aspect * Aspect);
   OldNorm += Old[DOFS_ - 1] * Old[DOFS_ - 1] / (Aspect * Aspect);
   NewNorm = sqrt(NewNorm);
   OldNorm = sqrt(OldNorm);

   return fabs(acos(angle / (NewNorm * OldNorm)));
}

int ArcLengthSolution::AllSolutionsFound() const
{
   ++counter_[8];

   return (CurrentSolution_ >= NumSolutions_);
}

int ArcLengthSolution::FindNextSolution()
{
   ++counter_[9];

   int good = 0;
   double AngleTest;
   // Assume that first solution should not be strictly restricted to the
   // adaptive steping angle constraint
   double AngleFactor = CurrentSolution_ ? 1.0 : 5.0;

   Vector OldDiff = Difference_;

   double forcenorm = 0.0;
   double dxnorm = 0.0;
   int Prediction = 1;
   int rotations = 1;
   int itr = 0;
   do
   {
      ArcLengthNewton(good, itr, forcenorm, dxnorm);

      AngleTest = ArcLenAngle(OldDiff, Difference_, Aspect_);

      cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << "\n";

      if (eig_angle_max_ >= 0.0)
      {
         cout << "Eigen-vector rotations (1=PASS, 0=FAIL) = "
              << (rotations = RelativeEigVectsOK())
              << "\n";
      }

      cout << "Prediction " << Prediction << " Corrector Iterations: " << itr << "\n";

      ++Prediction;
   }
   while (((AngleTest >= AngleFactor * AngleCutoff_) || (rotations == 0) || !good)
          && (CurrentDS_ >= DSMin_)
          && (ArcLenUpdate(-Difference_), // back to previous solution
              Difference_ = OldDiff,
              CurrentDS_ = CurrentDS_ / 2.0));

   CumulativeArcLength_ += CurrentDS_;

   cout << "Converged to solution at CumulativeArcLength = " << CumulativeArcLength_
        << " with CorrectorNorm = " << dxnorm
        << ",     ForceNorm = " << forcenorm << "\n";

   if (!good)
   {
      cerr << "ArcLenghtSolution did not converge properly" << "\n";
   }

   if ((MaxCumulativeArcLength_ >= 0.0) && (CumulativeArcLength_ > MaxCumulativeArcLength_))
   {
      // We are done -- set currentsolution to numsolutions-1
      // it will then be incremented to numsolutions below.
      cout << "NOTE: MaxCumulativeArcLength reached at Solution # " << CurrentSolution_
           << " --- Terminating!" << "\n";

      CurrentSolution_ = NumSolutions_ - 1;
   }

   if ((AngleTest <= AngleIncrease_) && (CurrentDS_ < DSMax_))
   {
      CurrentDS_ *= 2.0;
      if (CurrentDS_ > DSMax_)
      {
         CurrentDS_ = DSMax_;
      }
   }

   if (BifStartFlag_)
   {
      Vector diff = Difference_;
      diff[ArcLenDef().Dim() - 1] = 0.0;
      diff /= diff.Norm();
      cout << "Projection on BifTangent = " << Restrict_->DOF() * BifTangent_
           << ",     Angle (deg.) with BifTangent = "
           << acos(diff * BifTangent_) * (57.2957795130823) << "\n";
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

void ArcLengthSolution::ArcLengthNewton(int& good, int& itr, double& forcenorm, double& dxnorm)
{
   ++counter_[0];

   itr = 0;
   int Dim = ArcLenDef().Dim();
   double ForceNorm = 0.0;
   double Magnitude1 = 0.0;
   double Magnitude2 = 0.0;
   double Kappa = 0.0;
   double const eta = 0.1;

   Vector Dx(Dim),
   RHS(Dim);
   Matrix stif(Dim, Dim);

   // Predictor step
   cout << "Taking Predictor Step. CurrentDS = " << CurrentDS_ << "\n";
   ArcLenUpdate(Difference_);

   // Iterate until convergence

   itr++;
   // get stiffness first for efficiency
   stif = ArcLenStiffness(Difference_, Aspect_);
   RHS = -ArcLenForce(CurrentDS_, Difference_, Aspect_);
   ForceNorm = RHS.Norm();
#ifdef SOLVE_SVD
   Dx = SolveSVD(
      stif,
      RHS, MAXCONDITION, Echo_);
#else
   Dx = SolvePLU(stif, RHS);
#endif
   Magnitude1 = Dx.Norm();
   Magnitude2 = Magnitude1;

   cout << "\tCorrectorNorm = " << Magnitude2
        << " \tForceNorm = " << ForceNorm
        << "\n";

   int Converged = 0;
   do
   {
      itr++;

      ArcLenUpdate(Dx);
      Difference_ += Dx;
      // get stiffness first for efficiency
      stif = ArcLenStiffness(Difference_, Aspect_);
      RHS = -ArcLenForce(CurrentDS_, Difference_, Aspect_);
      ForceNorm = RHS.Norm();
#ifdef SOLVE_SVD
      Dx = SolveSVD(
         stif,
         RHS, MAXCONDITION, Echo_);
#else
      Dx = SolvePLU(stif, RHS);
#endif
      Magnitude2 = Dx.Norm();
      Kappa = Magnitude2 / (Magnitude1 + Tolerance_ * eta);
      Magnitude1 = Magnitude2;

      cout << "\tCorrectorNorm = " << Magnitude2
           << " \tForceNorm = " << ForceNorm
           << " \tContraction = " << Kappa << "\n";

      switch (ConvergeType_)
      {
         case Both:
            if ((RHS.Norm() <= RHS.Dim()*Tolerance_) && (Dx.Norm() <= Dx.Dim()*Tolerance_))
            {
               Converged = 1;
            }
            break;
         case Force:
            if (RHS.Norm() <= RHS.Dim()*Tolerance_)
            {
               Converged = 1;
            }
            break;
         case Displacement:
            if (Dx.Norm() <= Dx.Dim()*Tolerance_)
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

   forcenorm = RHS.Norm();
   dxnorm = Dx.Norm();

   good = 1;
}

void ArcLengthSolution::FindCriticalPoint(Lattice* const Lat, int* const TotalNumCPCrossings,
                                          PerlInput const& Input, int const& Width, ostream& out)
{
   ++counter_[10];

   Vector OriginalDOF = ArcLenDef();
   Vector OriginalDiff = Difference_;
   double OriginalDS = CurrentDS_;
   int TestValueDiff;
   int temp;
   int size = Lat->NumTestFunctions();
   TF_LHS_static.Resize(size);
   TF_RHS_static.Resize(size);
   CurrentTF_static.Resize(size);
   double fa, fb;
   double LHSLambda, RHSLambda;
   int Multiplicity;
   int track;
   int num;
   int CP;
   int spot;
   int CPorBif;
   int Bif;

   RHSLambda = ArcLenDef()[DOFS_ - 1];
   LHSLambda = RHSLambda - Difference_[DOFS_ - 1];

   TestValueDiff = Lat->TestFunctions(TF_LHS_static, Lattice::RHS, &TF_RHS_static);


   int* Index;
   Index = new int[TestValueDiff];
   temp = 0;
   for (int i = 0; i < size; i++)
   {
      if ((TF_LHS_static[i] * TF_RHS_static[i]) < 0.0)
      {
         Index[temp] = i;
         temp++;
      }
   }

   int* ModeType;
   ModeType = new int[TestValueDiff];
   Matrix const& LHSEigVect = Lat->LHSEigVect();
   Vector Mode(LHSEigVect.Rows() + 1);
   for (int i = 0; i < TestValueDiff; ++i)
   {
      if ((Lat->UseEigenValTFs() == 0) || (Index[i] >= Lat->DOF().Dim()))
      {
         ModeType[i] = -1;
      }
      else
      {
         for (int j = 0; j < Mode.Dim() - 1; ++j)
         {
            Mode[j] = LHSEigVect[j][Index[i]];
         }
         Mode[Mode.Dim() - 1] = 0.0;
         // 1 = bif, 0 = turning pt;
         double nrm = (Restrict_->TransformVector(Mode)).Norm();
         if (nrm < 1.0e-9)
         {
            ModeType[i] = 1;
         }
         else
         {
            ModeType[i] = 0;
         }
      }
   }

   Vector DSTrack(TestValueDiff);
   Vector* CPDOFs;
   CPDOFs = new Vector[TestValueDiff];
   for (int i = 0; i < TestValueDiff; ++i)
   {
      CPDOFs[i].Resize(Lat->DOF().Dim());
   }
   double* CPLambdas;
   CPLambdas = new double[TestValueDiff];
   int* CPIndex;
   CPIndex = new int[TestValueDiff];
   int* CPType;
   CPType = new int[TestValueDiff];
   int* CPMultiplicity;
   CPMultiplicity = new int[TestValueDiff];

   cout << "TestValueDiff = " << TestValueDiff << "\n";
   for (int i = 0; i < TestValueDiff; i++)
   {
      cout << "Index[" << setw(2) << i << "] = " << setw(6) << Index[i]
           << ",   TF_LHS[" << setw(6) << Index[i] << "] = "
           << setw(Width) << TF_LHS_static[Index[i]]
           << ",   TF_RHS[" << setw(6) << Index[i] << "] = "
           << setw(Width) << TF_RHS_static[Index[i]] << "\n";
   }

   num = 0;
   for (CP = 0; CP < TestValueDiff; CP++)
   {
      track = Index[CP];
      fa = TF_LHS_static[track];
      fb = TF_RHS_static[track];

      if (track >= 0) // START OF IF STATEMENT
      {
         ZBrent(Lat, track, OriginalDiff, OriginalDS, fa, fb, CurrentTF_static);

         Multiplicity = 1;
         for (int i = CP + 1; i < TestValueDiff; i++)
         {
            temp = Index[i];
            if ((fabs(CurrentTF_static[temp]) <= Tolerance_) && (temp < Lat->DOF().Dim())) // only count multiplicity if NOT an ExtraTF
            {
               Index[i] = -1;
               Multiplicity++;
               // cout <<"i = " << i << "\n" <<  "CHECK POINT INDEX[CP] = " << Index[i] << "\n";
            }
         }

         // sort the critical points
         spot = num;
         while ((spot != 0) && (DSTrack[spot - 1] > CurrentDS_))
         {
            DSTrack[spot] = DSTrack[spot - 1];
            CPIndex[spot] = CPIndex[spot - 1];
            CPType[spot] = CPType[spot - 1];
            CPMultiplicity[spot] = CPMultiplicity[spot - 1];
            CPDOFs[spot] = CPDOFs[spot - 1];
            CPLambdas[spot] = CPLambdas[spot - 1];
            spot = spot - 1;
         }
         DSTrack[spot] = CurrentDS_;
         CPIndex[spot] = track;
         CPType[spot] = ModeType[CP];
         CPMultiplicity[spot] = Multiplicity;
         CPDOFs[spot] = Lat->DOF();
         CPLambdas[spot] = ((Lat->LoadParameter() == Lattice::Load)
                            ? Lat->Lambda() : Lat->Temp());

         num = num + 1;
      } // END OF IF STATEMENT
   }

   // //PRINT OUT CP DATA
   for (int i = 0; i < num; i++)
   {
      // reset to appropriate cp.
      Lat->SetDOF(CPDOFs[i]);
      if (Lat->LoadParameter() == Lattice::Load)
      {
         Lat->SetLambda(CPLambdas[i]);
      }
      else
      {
         Lat->SetTemp(CPLambdas[i]);
      }

      // Output Critical Point
      for (int j = 0; j < 70; j++)
      {
         if (Echo_)
         {
            cout << "=";
         }
         out << "=";
      }
      if (Echo_)
      {
         cout << "\n";
      }
      out << "\n";

      // determine CP type
      if ((Lat->UseEigenValTFs() == 0) || (CPIndex[i] >= Lat->DOF().Dim()))
      {
         CPorBif = -1; // ExtraTF
      }
      else
      {
         if (((CPLambdas[i] >= LHSLambda) && (CPLambdas[i] <= RHSLambda))
             || ((CPLambdas[i] <= LHSLambda) && (CPLambdas[i] >= RHSLambda)))
         {
            CPorBif = 1; // bif point
         }
         else
         {
            CPorBif = 0; // turning point
         }

         // Assume it suffices to check the first mode and don't bother with the multiplicity
         if (CPType[i] != CPorBif)
         {
            out << "NOTE: Conflict between critical point identification methods in ArcLengthSolution.\n"
                << "      Using characterization determined by Group Theory." << "\n";
            if (Echo_)
            {
               cout << "NOTE: Conflict between critical point identification methods in ArcLengthSolution.\n"
                    << "      Using characterization determined by Group Theory." << "\n";
            }
            CPorBif = CPType[i];
         }
      }

      // Lattice takes care of echo
      out << setw(Width);
      if (0 == CPorBif)
      {
         Lat->Print(out, Lattice::PrintShort, Lattice::TurningPt);
      }
      else if (1 == CPorBif)
      {
         Lat->Print(out, Lattice::PrintShort, Lattice::BifurcationPt);
      }
      else
      {
         Lat->Print(out, Lattice::PrintShort, Lattice::ExtraTFPt);
      }

      for (int j = 0; j < 70; j++)
      {
         if (Echo_)
         {
            cout << "=";
         }
         out << "=";
      }
      if (Echo_)
      {
         cout << "\n";
      }
      out << "\n";

      // Call Lattice function to do any Lattice Specific things
      Bif = Lat->CriticalPointInfo(TotalNumCPCrossings, CPIndex[i],
                                   Restrict_->DrDt(Restrict_->DOF() - (OriginalDOF - OriginalDiff)),
                                   CPorBif, CPMultiplicity[i], 10.0 * Tolerance_, Width, Input, out);

      if (Echo_)
      {
         cout << "Success = 1" << "\n";
      }
      out << "Success = 1" << "\n";

      ++TotalNumCPCrossings[CPIndex[i]];
   }

   delete[] Index;
   delete[] ModeType;
   delete[] CPDOFs;
   delete[] CPIndex;
   delete[] CPType;
   delete[] CPMultiplicity;
   delete[] CPLambdas;

   // Reset Lattice and ArcLengthSolution
   ArcLenSet(OriginalDOF);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   // Check to see if we should stop
   int cumulative = 0;
   for (int i = 0; i < Lat->NumTestFunctions(); ++i)
   {
      cumulative += TotalNumCPCrossings[i];
   }
   if ((StopAtCPCrossingNum_ > -1) && (cumulative >= StopAtCPCrossingNum_))
   {
      CurrentSolution_ = NumSolutions_;
   }
}

int ArcLengthSolution::ZBrent(Lattice* const Lat, int const& track, Vector const& OriginalDiff,
                              double const& OriginalDS, double& fa, double& fb, Vector& CurrentTF)
{
   ++counter_[1];

   int retcode = 1;
   Vector LastDiff(Difference_.Dim(), 0.0);
   double LastDS = CurrentDS_;
   double a, b, c, d, e, xm, p, fc, tol1, s, q, r, min1, min2;
   int good = 1;
   int loops = 0;
   double factor = 0.0;
   int oldprecision = cout.precision();
   int itr = 0;
   double forcenorm = 0.0;
   double dxnorm = 0.0;

   b = OriginalDS;
   c = b;
   a = 0.0;
   d = e = 0.0; // arbitrary initial values.

   fc = fb;
   // cout << " a = " << a << "\n" << "fa = " << fa << "\n"
   // <<  "b = " << b << "\n" << "fb = " << fb << "\n"
   // << "c = " << c << "\n" << "fc = " << fc << "\n";

   while (((fabs(fb) > Tolerance_)) && (loops < MaxIter_))
   {
      cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
      cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";

      ArcLenUpdate(-Difference_);

      if (((fb > 0.0) && (fc > 0.0)) || ((fb < 0.0) && (fc < 0.0)))
      {
         c = a;
         fc = fa;
         e = b - a;
         d = e;
      }

      if (fabs(fc) < fabs(fb))
      {
         a = b;
         b = c;
         c = a;
         fa = fb;
         fb = fc;
         fc = fa;
      }

      tol1 = 2.0* ARCLENEPS* fabs(b) + 0.5 * Tolerance_;
      xm = 0.5 * (c - b);

      if ((fabs(xm) <= tol1) || (fb == 0.0))
      {
         cout << "Minimal Root found! XM too small! " << "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         Lat->TestFunctions(CurrentTF, Lattice::INTERMED);
         fb = CurrentTF[track];

         cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";
         break;
      }

      if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb)))
      {
         s = fb / fa;
         if (a == c)
         {
            p = 2.0 * xm * s;
            q = 1.0 - s;
         }
         else
         {
            q = fa / fc;
            r = fb / fc;
            p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
            q = (q - 1.0) * (r - 1.0) * (s - 1.0);
         }

         if (p > 0.0)
         {
            q = -q;
         }

         p = fabs(p);
         min1 = 3.0 * xm * q - fabs(tol1 * q);
         min2 = fabs(e * q);

         if (2.0 * p < ((min1 < min2) ? min1 : min2))
         {
            e = d;
            d = p / q;
         }
         else
         {
            d = xm;
            e = d;
         }
      }
      else
      {
         d = xm;
         e = d;
      }

      a = b;
      fa = fb;
      if (fabs(d) > tol1)
      {
         b += d;
      }
      else
      {
         b += (xm >= 0.0 ? fabs(tol1) : -fabs(tol1));
      }
      LastDiff = Difference_;
      LastDS = CurrentDS_;
      CurrentDS_ = b;
      factor = OriginalDS / b;
      Difference_ = OriginalDiff / factor;
      ArcLengthNewton(good, itr, forcenorm, dxnorm);
      cout << "Converged with CorrectorNorm = " << dxnorm
           << ",     ForceNorm = " << forcenorm << "\n";

      if (!good)
      {
         // set back to best solution
         cout << "ZBrent stopped due to ArcLengthNewton() convergence failure! " << "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         fb = CurrentTF[track];
         cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ = " << CurrentDS_ << setprecision(oldprecision) << "\n";
         retcode = 0;
         break;
      }
      Lat->TestFunctions(CurrentTF, Lattice::INTERMED);

      fb = CurrentTF[track];
      loops++;

      // cout << "CurrentTF = " << "\n" << setw(Width) << CurrentTF << "\n" << "\n";
      // cout << "CurrentTF[track] = " << "\n" << CurrentTF[track]<< "\n"<< "\n";
   }

   if (loops >= MaxIter_)
   {
      cout << "Error: ZBrent reached Maximum number of iterations before it converged and exited." << "\n";
   }
   else
   {
      cout << "ZBrent finished with the " << track << "th TestFunction value of " << fb
           << "\n\n";
   }

   return retcode;
}

int ArcLengthSolution::RelativeEigVectsOK() const
{
   int retval = 1;

   if ((CurrentSolution_ > 0) && (eig_angle_max_ > 0.0)) // check enabled
   {
      double proj = cos(eig_angle_max_);
      Matrix RelEigVects = Restrict_->RelativeEigVects();
      int size = RelEigVects.Rows();
      
      for (int i = 0; i < size; ++i)
      {
         double maxval = fabs(RelEigVects[0][i]);
         int row = 0;
         for (int j = 0; j < size; ++j)
         {
            if (fabs(RelEigVects[j][i]) > maxval)
            {
               maxval = fabs(RelEigVects[j][i]);
               row = j;
            }
         }
         if ((row != i) || (maxval < proj))
         {
            retval = 0;
            break;
         }
      }
   }
   
   return retval;
}
