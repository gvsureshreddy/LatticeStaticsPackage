#include "ScanningSolution.h"
#include <math.h>

#include "UtilityFunctions.h"

ScanningSolution::ScanningSolution(LatticeMode *Mode,char *datafile)
   : Mode_(Mode)
{
   FILE *pipe;

   const char *yn[]={"No","Yes"};
   int ans;
   // Get parameters
   GetParameter("^ScanningMaxIterations",datafile,"%u",&MaxIter_);
   GetParameter("^ScanningTolerance",datafile,"%lf",&Tolerance_);
   GetParameter("^ScanningNewtonTolerance",datafile,"%lf",&NewtonTolerance_);
   ans=GetStringParameter("^ScanningFullField",datafile,yn,2);
   if (ans == 1)
      ScanFullField_ = Yes;
   else if (ans == 0)
      ScanFullField_ = No;
   else
   {
      cerr << "Error: Unknown ScanningFullField value!" << endl;
      exit(-1);
   }

   char tmp[LINELENGTH];
   char inidef[]="^ScanningInitialDeformation";
   SetPerlCommand(tmp,datafile,inidef);
   pipe=OpenPipe(tmp,"r");
   InitialDef_.Resize(Mode->ScanningRHS().Dim());
   InitialDef_=Mode->ScanningRHS();
   for (int i=0;i<InitialDef_.Dim();i++)
   {
      fscanf(pipe,"%lf",&InitialDef_[i]);
   }
   if (pclose(pipe)) Errfun(inidef);

   const char *dir[]={"Loading","Deformation"};
   ans=GetStringParameter("^ScanningDirection",datafile,dir,2);
   if (ans == 0)
      Direction_ = Loading;
   else if (ans == 1)
      Direction_ = Deformation;
   else
   {
      cerr << "Unknown Scanning direction" << endl;
      exit(-1);
   }

   GetParameter("^ScanningStart",datafile,"%lf",&ScanStart_);
   GetParameter("^ScanningEnd",datafile,"%lf",&ScanEnd_);
   GetParameter("^ScanningStep",datafile,"%lf",&ScanStep_);
   GetParameter("^ScanningLineStart",datafile,"%lf",&LineStart_);
   GetParameter("^ScanningLineStep",datafile,"%lf",&LineStep_);
   GetParameter("^ScanningLineEnd",datafile,"%lf",&LineEnd_);

   // Initialize Lattice to be ready to
   // find a solution
   CurrentScanLine_ = ScanStart_;
   OnSolution_ = No;
   InitializeLine();
   
};

void ScanningSolution::InitializeLine()
{
   if (Direction_ == Loading)
   {
      Mode_->ScanningLoadParamUpdate(
	 LineStart_ - Mode_->ScanningLoadParameter());
      
      Mode_->ScanningDefParamUpdate(
	 CurrentScanLine_ - Mode_->ScanningDefParameter());

      Mode_->ScanningUpdate(
	 Mode_->ScanningDef() - InitialDef_);
   }
   else
   {
      Mode_->ScanningLoadParamUpdate(
	 CurrentScanLine_ - Mode_->ScanningLoadParameter());

      Mode_->ScanningDefParamUpdate(
	 LineStart_ - Mode_->ScanningDefParameter());

      Mode_->ScanningUpdate(
	 Mode_->ScanningDef() - InitialDef_);
   }
}
   

int ScanningSolution::AllSolutionsFound()
{
   return CurrentScanLine_*(ScanStep_/fabs(ScanStep_))
      > ScanEnd_*(ScanStep_/fabs(ScanStep_));
}

int ScanningSolution::FindNextSolution()
{
   int iteration=0,
      good = 1;

   if (OnSolution_ == Yes)
   {
      if (ScanFullField_ == No)
      {
	 InitializeLine();
      }
      else
      {
	 if (Direction_ == Loading)
	    Mode_->ScanningLoadParamUpdate(LineStep_);
	 else
	    Mode_->ScanningDefParamUpdate(LineStep_);
      }
   }

   good = ScanningNewton();

   double val = Mode_->ScanningStressParameter(),
      sign = fabs(val)/val,
      newsign = 0;
   val = fabs(val);

   if (val < Tolerance_)
      return 1;

   // March along CurrentScanLine_ until
   // ScanningStressParameter changes sign
   while (sign*newsign >= 0)
   {
      if (Direction_ == Loading)
      {
	 if (Mode_->ScanningLoadParameter()*LineStep_/fabs(LineStep_)
	     > LineEnd_*LineStep_/fabs(LineStep_))
	 {
	    CurrentScanLine_ += ScanStep_;
	    OnSolution_ = No;
	    InitializeLine();
	    return 0;
	 }
	 
	 cout << Mode_->ScanningLoadParameter();
      }
      else
      {
	 if (Mode_->ScanningDefParameter()*LineStep_/fabs(LineStep_)
	     > LineEnd_*LineStep_/fabs(LineStep_))
	 {
	    CurrentScanLine_ += ScanStep_;
	    OnSolution_ = No;
	    InitializeLine();
	    return 0;
	 }
	 
	 cout << Mode_->ScanningDefParameter();
      }
      cout << "\t" << Mode_->ScanningStressParameter() << endl;

      if (Direction_ == Loading)
	 Mode_->ScanningLoadParamUpdate(LineStep_);
      else
	 Mode_->ScanningDefParamUpdate(LineStep_);

      good = (ScanningNewton() && good);

      val = Mode_->ScanningStressParameter();
      newsign = val/fabs(val);
   }

   // Bisection method
   while ((fabs(Mode_->ScanningStressParameter()) > Tolerance_)
	  && (iteration < MaxIter_))
   {
      if (Direction_ == Loading)
	 cout << Mode_->ScanningLoadParameter();
      else
	 cout << Mode_->ScanningDefParameter();
      cout << "\t" << Mode_->ScanningStressParameter() << endl;

      iteration++;

      if (Direction_ == Loading)
      {
	 Mode_->ScanningLoadParamUpdate(
	    LineStep_*sign*newsign/(pow(2.0,iteration)));
      }
      else
      {
	 Mode_->ScanningDefParamUpdate(
	    LineStep_*sign*newsign/(pow(2.0,iteration)));
      }

      good = (ScanningNewton() && good);

      val = Mode_->ScanningStressParameter();
      newsign =val/fabs(val);
   }

   if (iteration >= MaxIter_)
   {
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
   
   return 1;
}

int ScanningSolution::ScanningNewton()
{
   int good=1;
   int itr=0;
   Vector RHS=Mode_->ScanningRHS();
   Vector dx(RHS.Dim());

   if (Direction_ == Loading)
   {
      Mode_->ScanningLoadParamUpdate(-LineStep_);
   }
   else
   {
      Mode_->ScanningDefParamUpdate(-LineStep_);
   }

   dx = SolveSVD(Mode_->ScanningStiffness(),
		 RHS,
		 MAXCONDITION,1);
   Mode_->ScanningUpdate(dx);

   if (Direction_ == Loading)
   {
      Mode_->ScanningLoadParamUpdate(LineStep_);
   }
   else
   {
      Mode_->ScanningDefParamUpdate(LineStep_);
   }

   // Iterate until convergence
   while ((dx.Norm() > NewtonTolerance_*Mode_->ScanningDef().Norm()) &&
	  (itr < MaxIter_))
   {
      itr++;

      dx=SolveSVD(Mode_->ScanningStiffness(),
		  Mode_->ScanningRHS(),
		  MAXCONDITION,1);

      Mode_->ScanningUpdate(dx);
   }

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached -- ScanningNewton" << endl;
      good =0;
   }

   return good;
}
