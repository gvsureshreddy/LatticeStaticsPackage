#include <cmath>
#include "ScanningSolution.h"

using namespace std;

ScanningSolution::ScanningSolution(LatticeMode *Mode,
                                   int MaxIter,double Tolerance,double NewtonTolerance,
                                   YN ScanFullField,Vector &InitialDef,
                                   unsigned ScnDefParam,ScanDir Direction,
                                   double ScanStart,double ScanEnd,double ScanStep,
                                   double LineStart,double LineEnd,double LineStep,
                                   YN OnSolution,int Echo)
   : Echo_(Echo),
     Mode_(Mode),
     ModeDOFS_(Mode_->ModeDOF().Dim()),
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
     CurrentScanLine_(ScanStart)
{
   InitializeLine();
}

ScanningSolution::ScanningSolution(LatticeMode *Mode,PerlInput &Input,int Echo)
   : Echo_(Echo),
     Mode_(Mode)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   // Get parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ScanningSolution");
   MaxIter_ = Input.getUnsigned(Hash,"MaxIterations");
   Tolerance_ = Input.getDouble(Hash,"Tolerance");
   NewtonTolerance_ = Input.getDouble(Hash,"NewtonTolerance");
   const char *fullfld = Input.getString(Hash,"FullField");
   if (!strcmp("Yes",fullfld))
      ScanFullField_ = Yes;
   else if (!strcmp("No",fullfld))
      ScanFullField_ = No;
   else
   {
      cerr << "Error ScanningSolution: Unknown FullField string!" << "\n";
      exit(-1);
   }
   
   InitialDef_.Resize(ModeDOFS_);
   Input.getVector(InitialDef_,Hash,"InitialDeformation");
   
   const char *dir = Input.getString(Hash,"Direction");
   if (!strcmp("Loading",dir))
      Direction_ = Loading;
   else if (!strcmp("Deformation",dir))
      Direction_ = Deformation;
   else
   {
      cerr << "Error ScanningSolution: Unknown Direction string!" << "\n";
      exit(-1);
   }

   ScnDefParam_ = Input.getUnsigned(Hash,"DefParam");

   if (ScnDefParam_ >= (unsigned) ModeDOFS_-1)
   {
      cerr << "ScanningSolution: ScanningDefParam too big!" << "\n";
      exit(-1);
   }

   ScanStart_ = Input.getDouble(Hash,"Start");
   ScanEnd_ = Input.getDouble(Hash,"End");
   ScanStep_ = Input.getDouble(Hash,"Step");
   LineStart_ = Input.getDouble(Hash,"LineStart");
   LineEnd_ = Input.getDouble(Hash,"LineEnd");
   LineStep_ = Input.getDouble(Hash,"LineStep");
      
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

//----------------------------------------------------------------
double ScanningSolution::ScanningDefParameter()
{
   return Mode_->ModeDOF()[ScnDefParam_];
}

void ScanningSolution::ScanningDefParamSet(const double val)
{
   static Vector ModeDOF(ModeDOFS_);
   ModeDOF = Mode_->ModeDOF();
   ModeDOF[ScnDefParam_] = val;
   Mode_->SetModeDOF(ModeDOF);
}

void ScanningSolution::ScanningDefParamUpdate(const double newval)
{
   static Vector ModeDOF(ModeDOFS_);
   ModeDOF = Mode_->ModeDOF();
   ModeDOF[ScnDefParam_] += newval;
   Mode_->SetModeDOF(ModeDOF);
}

double ScanningSolution::ScanningLoadParameter()
{
   return Mode_->ModeDOF()[ModeDOFS_-1];
}

void ScanningSolution::ScanningLoadParamSet(const double val)
{
   static Vector ModeDOF(ModeDOFS_);
   ModeDOF = Mode_->ModeDOF();
   ModeDOF[ModeDOFS_-1] = val;
   Mode_->SetModeDOF(ModeDOF);
}

void ScanningSolution::ScanningLoadParamUpdate(const double newval)
{
   static Vector ModeDOF(ModeDOFS_);
   ModeDOF = Mode_->ModeDOF();
   ModeDOF[ModeDOFS_-1] += newval;
   Mode_->SetModeDOF(ModeDOF);
}

double ScanningSolution::ScanningStressParameter()
{
   return (Mode_->ModeForce())[ScnDefParam_];
}

Vector ScanningSolution::ScanningForce()
{
   Vector stress = Mode_->ModeForce();
   Vector force(ModeDOFS_-2,0.0);
   int a=0;
   
   if (ModeDOFS_ != 2)
   {
      for (unsigned i=0;i<ModeDOFS_-1;++i)
      {
         if (i != ScnDefParam_)
         {
            force[a] = stress[i];
            ++a;
         }
      }
      
      return force;
   }
   else
      return Vector(1,0.0);
}

Vector ScanningSolution::ScanningDef()
{
   static Vector ModeDOF(ModeDOFS_);
   static Vector DEF(ModeDOFS_-2,0.0);
   
   ModeDOF = Mode_->ModeDOF();
   int a=0;
   
   if (ModeDOFS_ != 2)
   {
      for (unsigned i=0;i<ModeDOFS_-1;++i)
      {
         if (i != ScnDefParam_)
         {
            DEF[a] = ModeDOF[i];
            ++a;
         }
      }
      
      return DEF;
   }
   else
      return Vector(1,0.0);
}

void ScanningSolution::ScanningSet(const Vector &val)
{
   static Vector ModeDOF(ModeDOFS_);
   
   ModeDOF = Mode_->ModeDOF();
   
   for (unsigned i=0;i<ModeDOFS_-1;++i)
   {
      if (i != ScnDefParam_)
      {
         ModeDOF[i] = val[i];
      }
   }
   
   Mode_->SetModeDOF(ModeDOF);
}

void ScanningSolution::ScanningUpdate(const Vector &newval)
{
   static Vector ModeDOF(ModeDOFS_);
   
   ModeDOF = Mode_->ModeDOF();
   
   for (unsigned i=0;i<ModeDOFS_-1;++i)
   {
      if (i != ScnDefParam_)
      {
         ModeDOF[i] += newval[(i>ScnDefParam_)?i-1:i];
      }
   }
   
   Mode_->SetModeDOF(ModeDOF);
}

Matrix ScanningSolution::ScanningStiffness()
{
   static Matrix ModeK(ModeDOFS_,ModeDOFS_+1);
   static Matrix K(ModeDOFS_-2,ModeDOFS_-2);
   
   int a=0,b=0;
   
   if (ModeDOFS_ != 2)
   {
      for (unsigned i=0;i<ModeDOFS_-1;++i)
      {
         if (i != ScnDefParam_)
         {
            for (unsigned j=0;j<ModeDOFS_-1;++j)
            {
               if (j != ScnDefParam_)
               {
                  K[a][b] = ModeK[i][j];
               }
               ++b;
            }
            ++a;
         }
      }
      return K;
   }
   else
      return Matrix(1,1,1.0);
}
//----------------------------------------------------------------

int ScanningSolution::AllSolutionsFound()
{
   return CurrentScanLine_*(ScanStep_/fabs(ScanStep_))
      > ScanEnd_*(ScanStep_/fabs(ScanStep_));
}

int ScanningSolution::FindNextSolution()
{
   int good=0;
   unsigned iteration=0;
   
   if (OnSolution_ == Yes)
   {
      if (ScanFullField_ == No)
      {
         InitializeLine();
      }
      else
      {
         if (Direction_ == Loading)
            ScanningLoadParamUpdate(LineStep_);
         else
            ScanningDefParamUpdate(LineStep_);
      }
   }
   
   ScanningNewton(good);
   
   double stepsize;
   double val = ScanningStressParameter(),
      oldval = val,
      sign = fabs(val)/val,
      newsign = 0;
   
   if (fabs(val) < Tolerance_)
   {
      good = 1;
      OnSolution_ = Yes;
      return good;
   }
   
   // March along CurrentScanLine_ until
   // ScanningStressParameter changes sign
   while (sign*newsign >= 0)
   {
      if (Direction_ == Loading)
      {
         if (ScanningLoadParameter()*LineStep_/fabs(LineStep_)
             > LineEnd_*LineStep_/fabs(LineStep_))
         {
            CurrentScanLine_ += ScanStep_;
            OnSolution_ = No;
            InitializeLine();
            good = 0;
            return good;
         }
         
         if (Echo_) cout << ScanningLoadParameter();
      }
      else
      {
         if (ScanningDefParameter()*LineStep_/fabs(LineStep_)
             > LineEnd_*LineStep_/fabs(LineStep_))
         {
            CurrentScanLine_ += ScanStep_;
            OnSolution_ = No;
            InitializeLine();
            good = 0;
            return good;
         }
         
         if (Echo_) cout << ScanningDefParameter();
      }
      if (Echo_) cout << "\t" << ScanningStressParameter() << "\n";
      
      if (Direction_ == Loading)
         ScanningLoadParamUpdate(LineStep_);
      else
         ScanningDefParamUpdate(LineStep_);
      
      ScanningNewton(good);
      
      oldval = val;
      val = ScanningStressParameter();
      newsign = val/fabs(val);
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
            cout << ScanningLoadParameter();
         else
            cout << ScanningDefParameter();
         cout << "\t" << ScanningStressParameter() << "\n";
      }
      
      iteration++;
      
      if (Direction_ == Loading)
      {
         // Bisection Method
         //ScanningLoadParamUpdate(
         //   -LineStep_*sign*newsign/(pow(2.0,iteration)));
         
         // Secant Method
         stepsize /= -(1.0 - oldval/val);
         ScanningLoadParamUpdate(stepsize);
      }
      else
      {
         // Bisection Method
         //ScanningDefParamUpdate(
         //   -LineStep_*sign*newsign/(pow(2.0,iteration)));
         
         // Secant Method
         stepsize /= -(1.0 - oldval/val);
         ScanningDefParamUpdate(stepsize);
      }
      
      ScanningNewton(good);
      
      oldval = val;
      val = ScanningStressParameter();
      newsign =val/fabs(val);
   }
   
   if (Echo_)
   {
      if (Direction_ == Loading)
         cout << ScanningLoadParameter();
      else
         cout << ScanningDefParameter();
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

void ScanningSolution::ScanningNewton(int &good)
{
   unsigned itr=0;
   
   Vector RHS=-ScanningForce();
   Vector dx(RHS.Dim());
   
   if (Direction_ == Loading)
   {
      ScanningLoadParamUpdate(-LineStep_);
   }
   else
   {
      ScanningDefParamUpdate(-LineStep_);
   }
   
#ifdef SOLVE_SVD
   dx = SolveSVD(ScanningStiffness(),
                 RHS,
                 MAXCONDITION,Echo_);
#else
   dx = SolvePLU(ScanningStiffness(),RHS);
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
      
#ifdef SOLVE_SVD
      dx=SolveSVD(ScanningStiffness(),
                  -ScanningForce(),
                  MAXCONDITION,Echo_);
#else
      dx=SolvePLU(ScanningStiffness(),-ScanningForce());
#endif
      
      ScanningUpdate(dx);
      if (Echo_) cout << "ScanningNewton(dx) = " << setw(20) << dx
                      << ", RHS = " << setw(20) << ScanningForce() << "\n";
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
