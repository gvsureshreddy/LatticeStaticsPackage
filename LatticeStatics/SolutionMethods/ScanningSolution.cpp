#include "ScanningSolution.h"
#include <math.h>

#include "UtilityFunctions.h"

ScanningSolution::ScanningSolution(LatticeMode *Mode,char *datafile,const char *prefix)
   : Mode_(Mode)
{
   FILE *pipe;

   const char *yn[]={"No","Yes"};
   int ans;
   // Get parameters
   if(!GetParameter(prefix,"ScanningMaxIterations",datafile,"%u",&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ScanningTolerance",datafile,"%lf",&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ScanningNewtonTolerance",datafile,"%lf",&NewtonTolerance_)) exit(-1);
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

   char tmp[LINELENGTH];
   char inidef[]="ScanningInitialDeformation";
   SetPerlCommand(tmp,datafile,prefix,inidef);
   pipe=OpenPipe(tmp,"r");
   InitialDef_.Resize(Mode->ScanningRHS().Dim());
   InitialDef_=Mode->ScanningRHS();
   for (int i=0;i<InitialDef_.Dim();i++)
   {
      fscanf(pipe,"%lf",&InitialDef_[i]);
   }
   if (pclose(pipe)) Errfun(inidef);

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

   if(!GetParameter(prefix,"ScanningStart",datafile,"%lf",&ScanStart_)) exit(-1);
   if(!GetParameter(prefix,"ScanningEnd",datafile,"%lf",&ScanEnd_)) exit(-1);
   if(!GetParameter(prefix,"ScanningStep",datafile,"%lf",&ScanStep_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineStart",datafile,"%lf",&LineStart_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineStep",datafile,"%lf",&LineStep_)) exit(-1);
   if(!GetParameter(prefix,"ScanningLineEnd",datafile,"%lf",&LineEnd_)) exit(-1);

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
	 Mode_->ScanningLoadParameter() - LineStart_);
      
      Mode_->ScanningDefParamUpdate(
	 Mode_->ScanningDefParameter() - CurrentScanLine_);

      Mode_->ScanningUpdate(
	 Mode_->ScanningDef() - InitialDef_);
   }
   else
   {
      Mode_->ScanningLoadParamUpdate(
	 Mode_->ScanningLoadParameter() - CurrentScanLine_);

      Mode_->ScanningDefParamUpdate(
	 Mode_->ScanningDefParameter() - LineStart_);

      Mode_->ScanningUpdate(
	 Mode_->ScanningDef() - InitialDef_);
   }
}
   

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
	    Mode_->ScanningLoadParamUpdate(-LineStep_);
	 else
	    Mode_->ScanningDefParamUpdate(-LineStep_);
      }
   }

   uncertainty = ScanningNewton(good);

   double stepsize;
   double val = Mode_->ScanningStressParameter(),
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
	 if (Mode_->ScanningLoadParameter()*LineStep_/fabs(LineStep_)
	     > LineEnd_*LineStep_/fabs(LineStep_))
	 {
	    CurrentScanLine_ += ScanStep_;
	    OnSolution_ = No;
	    InitializeLine();
	    good = 0;
	    return uncertainty;
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
	    good = 0;
	    return uncertainty;
	 }
	 
	 cout << Mode_->ScanningDefParameter();
      }
      cout << "\t" << Mode_->ScanningStressParameter() << endl;

      if (Direction_ == Loading)
	 Mode_->ScanningLoadParamUpdate(-LineStep_);
      else
	 Mode_->ScanningDefParamUpdate(-LineStep_);

      uncertainty = ScanningNewton(good);

      oldval = val;
      val = Mode_->ScanningStressParameter();
      newsign = val/fabs(val);
   }

   // Iterate onto solution
   stepsize = LineStep_;
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
	 // Bisection Method
	 //Mode_->ScanningLoadParamUpdate(
	 //   -LineStep_*sign*newsign/(pow(2.0,iteration)));

	 // Secant Method
	 stepsize /= -(1.0 - oldval/val);
	 Mode_->ScanningLoadParamUpdate(-stepsize);
      }
      else
      {
	 // Bisection Method
	 //Mode_->ScanningDefParamUpdate(
	 //   -LineStep_*sign*newsign/(pow(2.0,iteration)));

	 // Secant Method
	 stepsize /= -(1.0 - oldval/val);
	 Mode_->ScanningDefParamUpdate(-stepsize);
      }

      uncertainty = ScanningNewton(good);

      oldval = val;
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
   
   return (stepsize > uncertainty)?stepsize:uncertainty;
}

double ScanningSolution::ScanningNewton(int &good)
{
   int itr=0;
   Vector RHS=Mode_->ScanningRHS();
   Vector dx(RHS.Dim());

   if (Direction_ == Loading)
   {
      Mode_->ScanningLoadParamUpdate(LineStep_);
   }
   else
   {
      Mode_->ScanningDefParamUpdate(LineStep_);
   }
   
#ifdef SOLVE_SVD
   dx = SolveSVD(Mode_->ScanningStiffness(),
		 RHS,
		 MAXCONDITION,1);
#else
   dx = SolvePLU(Mode_->ScanningStiffness(),RHS);
#endif
   
   Mode_->ScanningUpdate(dx);

   if (Direction_ == Loading)
   {
      Mode_->ScanningLoadParamUpdate(-LineStep_);
   }
   else
   {
      Mode_->ScanningDefParamUpdate(-LineStep_);
   }

   // Iterate until convergence
   while ((dx.Norm() > NewtonTolerance_*Mode_->ScanningDef().Norm()) &&
	  (itr < MaxIter_))
   {
      itr++;

#ifdef SOLVE_SVD
      dx=SolveSVD(Mode_->ScanningStiffness(),
		  Mode_->ScanningRHS(),
		  MAXCONDITION,1);
#else
      dx=SolvePLU(Mode_->ScanningStiffness(),Mode_->ScanningRHS());
#endif
      
      Mode_->ScanningUpdate(dx);
      cout << "ScanningNewton(dx.Norm) = " << dx.Norm() << endl;
   }

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached -- ScanningNewton" << endl;
      good = 0;
   }

   return dx.Norm();
}
