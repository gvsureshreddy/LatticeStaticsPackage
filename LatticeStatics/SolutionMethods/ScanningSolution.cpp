#include <cmath>
#include "ScanningSolution.h"

using namespace std;

ScanningSolution::ScanningSolution(Restriction* const Restrict, int const& MaxIter,
                                   double const& Tolerance, double const& NewtonTolerance,
                                   YN const& ScanFullField, Vector const& InitialDef,
                                   int const& ScnDefParam, ScanDir const& Direction,
                                   double const& ScanStart, double const& ScanEnd,
                                   double const& ScanStep, double const& LineStart,
                                   double const& LineEnd, double const& LineStep,
                                   YN const& OnSolution, int const& Echo) :
   Echo_(Echo),
   Restrict_(Restrict),
   DOFS_(Restrict_->DOF().Dim()),
   DOF_(Restrict_->DOF()),
   MaxIter_(MaxIter),
   Tolerance_(Tolerance),
   NewtonTolerance_(NewtonTolerance),
   ScanFullField_(ScanFullField),
   OnSolution_(OnSolution),
   InitialDef_(InitialDef),
   ScnDefParam_(ScnDefParam),
   Direction_(Direction),
   ScanStart_(ScanStart),
   ScanEnd_(ScanEnd),
   ScanStep_(ScanStep),
   LineStart_(LineStart),
   LineEnd_(LineEnd),
   LineStep_(LineStep),
   CurrentScanLine_(ScanStart),
   stress_static(DOFS_ - 1),
   RestrictK_static(DOFS_ - 1, DOFS_)
{
   InitializeLine();
}

ScanningSolution::ScanningSolution(Restriction* const Restrict, PerlInput const& Input,
                                   int const& Echo) :
   Echo_(Echo),
   Restrict_(Restrict)
{
   DOFS_ = Restrict_->DOF().Dim();
   DOF_.Resize(DOFS_);
   DOF_ = Restrict_->DOF();
   // Get parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "ScanningSolution");
   MaxIter_ = Input.getPosInt(Hash, "MaxIterations");
   Tolerance_ = Input.getDouble(Hash, "Tolerance");
   NewtonTolerance_ = Input.getDouble(Hash, "NewtonTolerance");
   char const* const fullfld = Input.getString(Hash, "FullField");
   if (!strcmp("Yes", fullfld))
   {
      ScanFullField_ = Yes;
   }
   else if (!strcmp("No", fullfld))
   {
      ScanFullField_ = No;
   }
   else
   {
      cerr << "Error ScanningSolution: Unknown FullField string!" << "\n";
      exit(-1);
   }

   InitialDef_.Resize(DOFS_);
   Input.getVector(InitialDef_, Hash, "InitialDeformation");

   char const* const dir = Input.getString(Hash, "Direction");
   if (!strcmp("Loading", dir))
   {
      Direction_ = Loading;
   }
   else if (!strcmp("Deformation", dir))
   {
      Direction_ = Deformation;
   }
   else
   {
      cerr << "Error ScanningSolution: Unknown Direction string!" << "\n";
      exit(-1);
   }

   ScnDefParam_ = Input.getPosInt(Hash, "DefParam");

   if (ScnDefParam_ >= (int) DOFS_ - 1)
   {
      cerr << "ScanningSolution: ScanningDefParam too big!" << "\n";
      exit(-1);
   }

   ScanStart_ = Input.getDouble(Hash, "Start");
   ScanEnd_ = Input.getDouble(Hash, "End");
   ScanStep_ = Input.getDouble(Hash, "Step");
   LineStart_ = Input.getDouble(Hash, "LineStart");
   LineEnd_ = Input.getDouble(Hash, "LineEnd");
   LineStep_ = Input.getDouble(Hash, "LineStep");
   Input.EndofInputSection();

   // Initialize various data storage space
   stress_static.Resize(DOFS_ - 1);
   RestrictK_static.Resize(DOFS_ - 1, DOFS_);

   // Initialize Lattice to be ready to
   // find a solution
   CurrentScanLine_ = ScanStart_;
   OnSolution_ = No;
   InitializeLine();
}

void ScanningSolution::InitializeLine()
{
   if (Direction_ == Loading)
   {
      ScanningLoadParamSet(LineStart_);

      ScanningSet(InitialDef_);

      ScanningDefParamSet(CurrentScanLine_);
   }
   else
   {
      ScanningLoadParamSet(CurrentScanLine_);

      ScanningSet(InitialDef_);

      ScanningDefParamSet(LineStart_);
   }
}

// ----------------------------------------------------------------
double const& ScanningSolution::ScanningDefParameter() const
{
   return DOF_[ScnDefParam_];
}

void ScanningSolution::ScanningDefParamSet(double const& val)
{
   DOF_[ScnDefParam_] = val;
   Restrict_->SetDOF(DOF_);
}

void ScanningSolution::ScanningDefParamUpdate(double const& newval)
{
   DOF_[ScnDefParam_] += newval;
   Restrict_->SetDOF(DOF_);
}

double const& ScanningSolution::ScanningLoadParameter() const
{
   return DOF_[DOFS_ - 1];
}

void ScanningSolution::ScanningLoadParamSet(double const& val)
{
   DOF_[DOFS_ - 1] = val;
   Restrict_->SetDOF(DOF_);
}

void ScanningSolution::ScanningLoadParamUpdate(double const& newval)
{
   DOF_[DOFS_ - 1] += newval;
   Restrict_->SetDOF(DOF_);
}

double const& ScanningSolution::ScanningStressParameter() const
{
   return (Restrict_->Force())[ScnDefParam_];
}

Vector const& ScanningSolution::ScanningForce() const
{
   stress_static = Restrict_->Force();
   force_static.Resize(DOFS_ - 2, 0.0);
   int a = 0;

   if (DOFS_ != 2)
   {
      for (int i = 0; i < DOFS_ - 1; ++i)
      {
         if (i != ScnDefParam_)
         {
            force_static[a] = stress_static[i];
            ++a;
         }
      }
   }
   else
   {
      force_static.Resize(1, 0.0);
   }

   return force_static;
}

Vector const& ScanningSolution::ScanningDef() const
{
   DEF_static.Resize(DOFS_ - 2, 0.0);

   int a = 0;

   if (DOFS_ != 2)
   {
      for (int i = 0; i < DOFS_ - 1; ++i)
      {
         if (i != ScnDefParam_)
         {
            DEF_static[a] = DOF_[i];
            ++a;
         }
      }
   }
   else
   {
      DEF_static.Resize(1, 0.0);
   }

   return DEF_static;
}

void ScanningSolution::ScanningSet(Vector const& val)
{
   for (int i = 0; i < DOFS_ - 1; ++i)
   {
      if (i != ScnDefParam_)
      {
         DOF_[i] = val[i];
      }
   }

   Restrict_->SetDOF(DOF_);
}

void ScanningSolution::ScanningUpdate(Vector const& newval)
{
   for (int i = 0; i < DOFS_ - 1; ++i)
   {
      if (i != ScnDefParam_)
      {
         DOF_[i] += newval[(i > ScnDefParam_) ? i - 1 : i];
      }
   }

   Restrict_->SetDOF(DOF_);
}

Matrix const& ScanningSolution::ScanningStiffness() const
{
   RestrictK_static = Restrict_->Stiffness();
   K_static.Resize(DOFS_ - 2, DOFS_ - 2);

   int a = 0, b = 0;

   if (DOFS_ != 2)
   {
      for (int i = 0; i < DOFS_ - 1; ++i)
      {
         if (i != ScnDefParam_)
         {
            b = 0;
            for (int j = 0; j < DOFS_ - 1; ++j)
            {
               if (j != ScnDefParam_)
               {
                  K_static[a][b] = RestrictK_static[i][j];
                  ++b;
               }
            }
            ++a;
         }
      }
   }
   else
   {
      K_static.Resize(1, 1, 1.0);
   }

   return K_static;
}
// ----------------------------------------------------------------

int ScanningSolution::AllSolutionsFound() const
{
   return CurrentScanLine_ * (ScanStep_ / fabs(ScanStep_))
          > ScanEnd_ * (ScanStep_ / fabs(ScanStep_));
}

int ScanningSolution::FindNextSolution()
{
   int good = 0;
   int iteration = 0;

   if (OnSolution_ == Yes)
   {
      if (ScanFullField_ == No)
      {
         InitializeLine();
      }
      else
      {
         if (Direction_ == Loading)
         {
            ScanningLoadParamUpdate(LineStep_);
         }
         else
         {
            ScanningDefParamUpdate(LineStep_);
         }
      }
   }

   ScanningNewton(good);

   double stepsize;
   double val = ScanningStressParameter(),
          oldval = val,
          sign = fabs(val) / val,
          newsign = 0;

   if (fabs(val) < Tolerance_)
   {
      good = 1;
      OnSolution_ = Yes;
      return good;
   }

   // March along CurrentScanLine_ until
   // ScanningStressParameter changes sign
   while (sign * newsign >= 0)
   {
      if (Direction_ == Loading)
      {
         if (ScanningLoadParameter() * LineStep_ / fabs(LineStep_)
             > LineEnd_ * LineStep_ / fabs(LineStep_))
         {
            CurrentScanLine_ += ScanStep_;
            OnSolution_ = No;
            InitializeLine();
            good = 0;
            return good;
         }

         if (Echo_)
         {
            cout << ScanningLoadParameter();
         }
      }
      else
      {
         if (ScanningDefParameter() * LineStep_ / fabs(LineStep_)
             > LineEnd_ * LineStep_ / fabs(LineStep_))
         {
            CurrentScanLine_ += ScanStep_;
            OnSolution_ = No;
            InitializeLine();
            good = 0;
            return good;
         }

         if (Echo_)
         {
            cout << ScanningDefParameter();
         }
      }
      if (Echo_)
      {
         cout << "\t" << ScanningStressParameter() << "\n";
      }

      if (Direction_ == Loading)
      {
         ScanningLoadParamUpdate(LineStep_);
      }
      else
      {
         ScanningDefParamUpdate(LineStep_);
      }

      ScanningNewton(good);

      oldval = val;
      val = ScanningStressParameter();
      newsign = val / fabs(val);
   }

   // Iterate onto solution
   stepsize = LineStep_;
   while ((fabs(ScanningStressParameter()) > Tolerance_)
          && (fabs(val - oldval) > Tolerance_)
          && (iteration < MaxIter_))
   {
      if (Echo_)
      {
         if (Direction_ == Loading)
         {
            cout << ScanningLoadParameter();
         }
         else
         {
            cout << ScanningDefParameter();
         }
         cout << "\t" << ScanningStressParameter() << "\n";
      }

      iteration++;

      if (Direction_ == Loading)
      {
         // Bisection Method
         // ScanningLoadParamUpdate(
         //   -LineStep_*sign*newsign/(pow(2.0,iteration)));

         // Secant Method
         stepsize /= -(1.0 - oldval / val);
         ScanningLoadParamUpdate(stepsize);
      }
      else
      {
         // Bisection Method
         // ScanningDefParamUpdate(
         //   -LineStep_*sign*newsign/(pow(2.0,iteration)));

         // Secant Method
         stepsize /= -(1.0 - oldval / val);
         ScanningDefParamUpdate(stepsize);
      }

      ScanningNewton(good);

      oldval = val;
      val = ScanningStressParameter();
      newsign = val / fabs(val);
   }

   if (Echo_)
   {
      if (Direction_ == Loading)
      {
         cout << ScanningLoadParameter();
      }
      else
      {
         cout << ScanningDefParameter();
      }
      cout << "\t" << ScanningStressParameter() << "\n";
   }

   if (iteration >= MaxIter_)
   {
      good = 0;
      cerr << "Final Convergence Not Reached -- ScanningSolution" << "\n";
   }

   if (!good)
   {
      cerr << "ScanningNewton did not converge -- ScanningSolution" << "\n";
   }

   OnSolution_ = Yes;
   if (ScanFullField_ == No)
   {
      CurrentScanLine_ += ScanStep_;
   }

   return good;
}

void ScanningSolution::ScanningNewton(int& good)
{
   int itr = 0;

   Vector RHS = -ScanningForce();
   Vector dx(RHS.Dim());

   if (Direction_ == Loading)
   {
      ScanningLoadParamUpdate(-LineStep_);
   }
   else
   {
      ScanningDefParamUpdate(-LineStep_);
   }

   Matrix Stiff = ScanningStiffness();
#ifdef SOLVE_SVD
   dx = SolveSVD(Stiff,
                 RHS,
                 MAXCONDITION, Echo_);
#else
   dx = SolvePLU(Stiff, RHS);
#endif

   ScanningUpdate(dx);

   if (Direction_ == Loading)
   {
      ScanningLoadParamUpdate(LineStep_);
   }
   else
   {
      ScanningDefParamUpdate(LineStep_);
   }

   // Iterate until convergence
   // 1/05 changed from relative criterion to absolute stopping criterion
   while ((dx.Norm() > NewtonTolerance_) && (itr < MaxIter_))
   {
      itr++;

      // get stiffness first for efficiency
      Stiff = ScanningStiffness();
      RHS = -ScanningForce();

      if (Echo_)
      {
         cout << "ScanningNewton(dx) = " << setw(20) << dx
              << ", RHS = " << setw(20) << RHS << "\n";
      }

#ifdef SOLVE_SVD
      dx = SolveSVD(Stiff, RHS, MAXCONDITION, Echo_);
#else
      dx = SolvePLU(Stiff, RHS);
#endif

      ScanningUpdate(dx);
   }

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached -- ScanningNewton" << "\n";
      good = 0;
   }
   else
   {
      good = 1;
   }
}

