#include "ArcLengthSolution.h"


#include "UtilityFunctions.h"

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,
				     const Vector &one,const Vector &two)
   : Mode_(Mode),Difference_(two-one), CurrentSolution_(0)
{
   FILE *pipe;
   char command[LINELENGTH];

   char itr[]="^ArcLenMaxIterations";
   SetPerlCommand(command,datafile,itr);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%u",&MaxIter_);
   if (pclose(pipe)) Errfun(itr);

   char tol[]="^ArcLenTolerance";
   SetPerlCommand(command,datafile,tol);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&Tolerance_);
   if (pclose(pipe)) Errfun(tol);

   char dsmax[]="^ArcLenDSMax";
   SetPerlCommand(command,datafile,dsmax);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&DSMax_);
   if (pclose(pipe)) Errfun(dsmax);
   CurrentDS_ = DSMax_;

   char dsmin[]="^ArcLenDSMin";
   SetPerlCommand(command,datafile,dsmin);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&DSMin_);
   if (pclose(pipe)) Errfun(dsmin);

   char angle[]="^ArcLenAngleCutoff";
   SetPerlCommand(command,datafile,angle);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&AngleCutoff_);
   if (pclose(pipe)) Errfun(angle);

   char angleinc[]="^ArcLenAngleIncrease";
   SetPerlCommand(command,datafile,angleinc);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&AngleIncrease_);
   if (pclose(pipe)) Errfun(angleinc);

   char aspect[]="^ArcLenAspect";
   SetPerlCommand(command,datafile,aspect);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%lf",&Aspect_);
   if (pclose(pipe)) Errfun(aspect);

   char nosol[]="^ArcLenNumSolutions";
   SetPerlCommand(command,datafile,nosol);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%u",&NumSolutions_);
   if (pclose(pipe)) Errfun(nosol);

   // Set Lattice to solution "two"
   Mode_->ArcLenUpdate(
      two - Mode_->ArcLenDef());
}

int ArcLengthSolution::AllSolutionsFound()
{
   return (CurrentSolution_ == NumSolutions_);
}

int ArcLengthSolution::FindNextSolution()
{
   int good=1;

   double AngleTest;

   Vector OldDiff = Difference_;

   do
   {
      cout << "DS= " << CurrentDS_ << endl;

      good = (ArcLengthNewton() && good);

      AngleTest = Mode_->ArcLenAngle(OldDiff,Difference_,Aspect_);
   }
   while ((AngleTest >= AngleCutoff_)
	  && (CurrentDS_ >= DSMin_)
	  && (CurrentDS_ = CurrentDS_/2.0));

   if ((AngleTest >= AngleIncrease_) && (CurrentDS_ < DSMax_))
   {
      CurrentDS_ *= 2.0;
      cout << "DS= " << CurrentDS_ << endl;
   }

   if (!good)
   {
      cerr << "ArcLenghtSolution did not converge properly" << endl;
   }

   CurrentSolution_++;

   return 1;
}

int ArcLengthSolution::ArcLengthNewton()
{
   int good = 1;
   int itr = 0;

   Vector Dx(Mode_->ArcLenDef().Dim());
   
   Mode_->ArcLenUpdate(-Difference_);

   // Iterate until convergence
   cout << "ArcLenNewton: Number of Iterations --\n";
   do
   {
      itr++;

      Dx = SolveSVD(
	 Mode_->ArcLenStiffness(Difference_,Aspect_),
	 Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_),
	 MAXCONDITION,1);

      Mode_->ArcLenUpdate(Dx);
      Difference_ -= Dx;

      cout << itr << "(" << setw(20)
	   << Mode_->ScanningStressParameter() << "), ";
   }
   while ((itr < MaxIter_)
	  && (fabs(Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_).Norm())
	      > Tolerance_)
	  && (fabs(Dx.Norm()) > Tolerance_));

   cout << endl;

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << endl;
      good = 0;
   }

   return good;
}

int ArcLengthSolution::BisectAlert(Lattice *Lat,int Width,fstream &out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   double CurrentMinEV,OldMinEV;
   int loops = 0;
   int CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);
   int OriginalNulity = CurrentNulity;
   int OldNulity;

   // Set Lattice back to previous solution
   Mode_->ArcLenUpdate(Difference_);
   OldNulity = Lat->StiffnessNulity(&OldMinEV);
   
   cout << "\t" << setw(Width) << OldNulity << setw(Width) << OldMinEV << endl;

   while ((fabs(CurrentMinEV - OldMinEV) > Tolerance_)
	  && (loops < MaxIter_))
   {
      cout << setw(Width) << CurrentNulity
	   <<setw(Width) << CurrentMinEV << endl;

      loops++;

      if ((OldNulity - CurrentNulity) != 0)
	 CurrentDS_ = -CurrentDS_/2.0;
      else
	 CurrentDS_ = CurrentDS_/2.0;

      ArcLengthNewton();
      OldMinEV = CurrentMinEV;
      OldNulity = CurrentNulity;
      CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);
   }

   // Output Critical Point
   for (int i=0;i<70;i++)
   {
      cout << "=";
      out << "=";
   }
   cout << endl; out << endl;

   cout << setw(Width) << Lat;
   out << setw(Width) << Lat;
      
   for (int i=0;i<70;i++)
   {
      cout << "=";
      out << "=";
   }
   cout << endl; out << endl;

   // Reset Lattice and ArcLengthSolution
   Mode_->ArcLenUpdate(-(OriginalDiff - Difference_));
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   return 1;
}
