#include "ScanningSolution.h"
#include <cmath>

#include "UtilityFunctions.h"

using namespace std;

ScanningSolution::ScanningSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   int Echo)
   : Mode_(Mode), Echo_(Echo)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   const char *yn[]={"No","Yes"};
   int ans;
   // Get parameters
   if(!GetParameter(prefix,"ScanningMaxIterations",datafile,'u',&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ScanningTolerance",datafile,'l',&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ScanningNewtonTolerance",datafile,'l',&NewtonTolerance_)) exit(-1);
   ans=GetStringParameter(prefix,"ScanningFullField",datafile,yn,2);
   if (ans == 1)
      ScanFullField_ = Yes;
   else if (ans == 0)
      ScanFullField_ = No;
   else
   {
      cerr << "Error: Unknown ScanningFullField value!" << endl;
      exit(-1);
   }
   
   InitialDef_.Resize(ModeDOFS_);
   if(!GetVectorParameter(prefix,"ScanningInitialDeformation",datafile,&InitialDef_)) exit(-1);
   
   const char *dir[]={"Loading","Deformation"};
   ans=GetStringParameter(prefix,"ScanningDirection",datafile,dir,2);
   if (ans == 0)
      Direction_ = Loading;
   else if (ans == 1)
      Direction_ = Deformation;
   else
   {
      cerr << "Unknown Scanning direction" << endl;
      exit(-1);
   }
   
   if (!GetParameter(prefix,"ScanningDefParam",datafile,'u',&ScnDefParam_))
      exit(-1);
   if ((ScnDefParam_ < 0) || (ScnDefParam_ >= ModeDOFS_-1))
   {
      cerr << "ScanningSolution: ScanningDefParam too small or too big!" << endl;
      exit(-1);
   }
   
   if(!GetParameter(prefix,"ScanningStart",datafile,'l',&ScanStart_)) exit(-1);
   if(!GetParameter(prefix,"ScanningEnd",datafile,'l',&ScanEnd_)) exit(-1);
   if(!GetParameter(prefix,"ScanningStep",datafile,'l',&ScanStep_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineStart",datafile,'l',&LineStart_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineStep",datafile,'l',&LineStep_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineEnd",datafile,'l',&LineEnd_)) exit(-1);
   
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
      for (int i=0;i<ModeDOFS_-1;++i)
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
      for (int i=0;i<ModeDOFS_-1;++i)
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
   
   for (int i=0;i<ModeDOFS_-1;++i)
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
   
   for (int i=0;i<ModeDOFS_-1;++i)
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
      for (int i=0;i<ModeDOFS_-1;++i)
      {
	 if (i != ScnDefParam_)
	 {
	    for (int j=0;j<ModeDOFS_-1;++j)
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

double ScanningSolution::FindNextSolution(int &good)
{
   double uncertainty;
   int iteration=0;
   
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
   
   uncertainty = ScanningNewton(good);
   
   double stepsize;
   double val = ScanningStressParameter(),
      oldval = val,
      sign = fabs(val)/val,
      newsign = 0;
   
   if (fabs(val) < Tolerance_)
   {
      good = 1;
      OnSolution_ = Yes;
      return uncertainty;
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
	    return uncertainty;
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
	    return uncertainty;
	 }
	 
	 if (Echo_) cout << ScanningDefParameter();
      }
      if (Echo_) cout << "\t" << ScanningStressParameter() << endl;
      
      if (Direction_ == Loading)
	 ScanningLoadParamUpdate(LineStep_);
      else
	 ScanningDefParamUpdate(LineStep_);
      
      uncertainty = ScanningNewton(good);
      
      oldval = val;
      val = ScanningStressParameter();
      newsign = val/fabs(val);
   }
   
   // Iterate onto solution
   stepsize = LineStep_;
   while ((fabs(ScanningStressParameter()) > Tolerance_)
	  && (iteration < MaxIter_))
   {
      if (Echo_)
      {
	 if (Direction_ == Loading)
	    cout << ScanningLoadParameter();
	 else
	    cout << ScanningDefParameter();
	 cout << "\t" << ScanningStressParameter() << endl;
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
      
      uncertainty = ScanningNewton(good);
      
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
      cout << "\t" << ScanningStressParameter() << endl;
   }
   
   if (iteration >= MaxIter_)
   {
      good = 0;
      cerr << "Final Convergence Not Reached -- ScanningSolution" << endl;
   }
   
   if (!good)
   {
      cerr << "ScanningNewton did not converge -- ScanningSolution" << endl;
   }
   
   OnSolution_ = Yes;
   if (ScanFullField_ == No)
   {
      CurrentScanLine_ += ScanStep_;
   }
   
   return (stepsize > uncertainty)?stepsize:uncertainty;
}

double ScanningSolution::ScanningNewton(int &good)
{
   int itr=0;
   
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
		      << ", RHS = " << setw(20) << ScanningForce() << endl;
   }
   
   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached -- ScanningNewton" << endl;
      good = 0;
   }
   else
   {
      good = 1;
   }
   
   return dx.Norm();
}
