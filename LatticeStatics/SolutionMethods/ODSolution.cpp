#include <fstream>
#include "ODSolution.h"
#include "Matrix.h"

using namespace std;

ODSolution::ODSolution(Restriction* const Restrict, PerlInput const& Input,
                       int const& Echo) :
   Restrict_(Restrict),
   Echo_(Echo),
   SolutionFound_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "ODSolution");
   Converge_ = Input.getDouble(Hash, "ConvergeCriteria");

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

   Width_ = Input.getInt("Main", "FieldWidth");

   DimDOFS_ = (Restrict_->DOF()).Dim();

   X0_.Resize(DimDOFS_);
   X_.Resize(DimDOFS_);
   K1_.Resize(DimDOFS_);
   K2_.Resize(DimDOFS_);
   K3_.Resize(DimDOFS_);
   K4_.Resize(DimDOFS_);
   Delta_.Resize(DimDOFS_, 100);
   CurrentSolution_.Resize(DimDOFS_);

   if (Input.ParameterOK(Hash, "InitialState"))
   {
      Input.getVector(X0_, Hash, "InitialState");
   }
   Restrict_->SetDOF(X0_);
   Potential0_ = Restrict_->Energy();
   Tangent_.Resize(DimDOFS_);
   if (Input.ParameterOK(Hash, "InitialTangent"))
   {
      Input.getVector(Tangent_, Hash, "InitialTangent");
   }

   Epsilon_ = Input.getDouble(Hash, "Epsilon");
   dt_ = Input.getDouble(Hash, "dt");

   X_ = X0_ + Epsilon_ * Tangent_;
   CurrentSolution_ = X_;
   Restrict_->SetDOF(CurrentSolution_);
   Potential1_ = Restrict_->Energy();
   ForceNorm_ = (Restrict_->Force()).Norm();
   t_ = 0.0 + dt_;
   // Solution at t_ = 0.0 and t_ = dt_ is X0_ and X_ respectively;
}

int ODSolution::FindNextSolution(PerlInput const& Input, int const& Width, ostream& out)
{
   Vector yn(DimDOFS_);
   Vector yTemp(DimDOFS_);
   Vector Force(DimDOFS_ - 1);

   yn = CurrentSolution_;

   // Calculate K1_;
   Force = dt_ * Restrict_->Force();
   for (int i = 0; i < (DimDOFS_ - 1); ++i)
   {
      K1_[i] = -Force[i];
   }
   K1_[DimDOFS_-1] = 0.0; //X0_[DimDOFS_-1];

   // Calculate K2_;
   yTemp = yn + (0.5 * K1_);
   Restrict_->SetDOF(yTemp);
   Force = dt_ * Restrict_->Force();
   for (int i = 0; i < (DimDOFS_ - 1); ++i)
   {
      K2_[i] = -Force[i];
   }
   K2_[DimDOFS_-1] = 0.0; //X0_[DimDOFS_-1];

   // Calculate K3_;
   yTemp = yn + (0.5 * K2_);
   Restrict_->SetDOF(yTemp);
   Force = dt_ * Restrict_->Force();
   for (int i = 0; i < (DimDOFS_ - 1); ++i)
   {
      K3_[i] = -Force[i];
   }
   K3_[DimDOFS_-1] = 0.0; //X0_[DimDOFS_-1];

   // Calculate K4_;
   yTemp = yn + K3_;
   Restrict_->SetDOF(yTemp);
   Force = dt_ * Restrict_->Force();
   for (int i = 0; i < (DimDOFS_ - 1); ++i)
   {
      K4_[i] = -Force[i];
   }
   K4_[DimDOFS_-1] = 0.0; //X0_[DimDOFS_-1];

   Delta_ = (1.0 / 6.0) * (K1_ + 2.0 * K2_ + 2.0 * K3_ + K4_);

   CurrentSolution_ = yn + Delta_;
   t_ = t_ + dt_;

   Restrict_->SetDOF(CurrentSolution_);
   Potential_ = Restrict_->Energy();
   ForceNorm_ = (Restrict_->Force()).Norm();

   cout << "ForceNorm_ = " << ForceNorm_ << "\n";

   SolutionFound_++;

   if (SolutionFound_ == 1)
   {
      if (Echo_)
      {
         cout << "t = " << 0.0 << "\n";
         cout << "Potential = " << Potential0_ << "\n";
         cout << "t = " << dt_ << "\n";
         cout << "Potential = " << Potential1_ << "\n";
         //cout << "DOF = " << setw(Width) << CurrentSolution_ << "\n";
      }
      out << "t = " << 0.0 << "\n";
      out << "Potential = " << Potential0_ << "\n";
      out << "DOF = " << setw(Width) << X0_ << "\n";  
      out << "t = " << dt_ << "\n";
      out << "Potential = " << Potential1_ << "\n";
      out << "DOF = " << setw(Width) << X_ << "\n";
   }
   
   if (Echo_)
   {
      cout << "t = " << t_ << "\n";
      cout << "Potential = " << Potential_ << "\n";
      //cout << "DOF = " << setw(Width) << CurrentSolution_ << "\n";
   }
   out << "t = " << t_ << "\n";
   out << "Potential = " << Potential_ << "\n";		   
   out << "DOF = " << setw(Width) << CurrentSolution_ << "\n";
   
   return 1;
}

int ODSolution::AllSolutionsFound() const
{
   if ((ForceNorm_ > Converge_) || (t_ < (2 * dt_)))
   {
      return 0;
   }
   if ((ForceNorm_ <= Converge_) || ((Delta_.Norm()) <= Converge_))
   {
      return 1;
   }
}

double ODSolution::dt() const
{
   return dt_;
}

double ODSolution::Time() const
{
   return t_;
}

double ODSolution::Energy() const
{
   return Potential_;
}

double ODSolution::Energy0() const
{
   return Potential0_;
}

double ODSolution::Energy1() const
{
   return Potential1_;
}

Vector const& ODSolution::Solution()
{
   return CurrentSolution_;
}

Vector const& ODSolution::X0()
{
   return X0_;
}

Vector const& ODSolution::X1()
{
   return X_;
}
