#include "ArcLengthSolution.h"


#include "UtilityFunctions.h"

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,
				     const Vector &one,const Vector &two)
   : Mode_(Mode),Difference_(two-one), CurrentSolution_(0)
{
   GetParameter("^ArcLenMaxIterations",datafile,"%u",&MaxIter_);
   GetParameter("^ArcLenTolerance",datafile,"%lf",&Tolerance_);
   GetParameter("^ArcLenDSMax",datafile,"%lf",&DSMax_);
   CurrentDS_ = DSMax_;
   GetParameter("^ArcLenDSMin",datafile,"%lf",&DSMin_);
   GetParameter("^ArcLenAngleCutoff",datafile,"%lf",&AngleCutoff_);
   GetParameter("^ArcLenAngleIncrease",datafile,"%lf",&AngleIncrease_);
   GetParameter("^ArcLenAspect",datafile,"%lf",&Aspect_);
   GetParameter("^ArcLenNumSolutions",datafile,"%u",&NumSolutions_);

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
	  && (Mode_->ArcLenUpdate(Difference_),CurrentDS_ = CurrentDS_/2.0));

   if ((AngleTest <= AngleIncrease_) && (CurrentDS_ < DSMax_))
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
	  && ((fabs(Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_).Norm())
	       > Tolerance_)
	      || (fabs(Dx.Norm()) > Tolerance_)));

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
   Vector IntermediateDiff(Difference_.Dim(),0.0);
   double OriginalDS = CurrentDS_;
   double CurrentMinEV,OldMinEV;
   int loops = 0;
   int OldNulity = Lat->StiffnessNulity(&OldMinEV);
   int CurrentNulity;

   // Set Lattice back to previous solution
   Mode_->ArcLenUpdate(Difference_);
   CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);
   
   cout << "\t" << setw(Width) << OldNulity << setw(Width) << OldMinEV
	<< " DS " << setw(Width) << CurrentDS_ << endl;

   while ((fabs(CurrentMinEV - OldMinEV) > Tolerance_)
	  && (loops < MaxIter_))
   {
      cout << setw(Width) << CurrentNulity
	   << setw(Width) << CurrentMinEV
	   << " DS " << setw(Width) << CurrentDS_ << endl;

      //CurrentDS_ /= 2.0; // Bisection Method
      CurrentDS_ /= (1- (OldMinEV/CurrentMinEV)); // Secant Method
      if (((OldNulity - CurrentNulity) != 0)
	  && (loops != 0))
	 Difference_ = -Difference_;

      ArcLengthNewton();
      IntermediateDiff += Difference_;
      OldMinEV = CurrentMinEV;
      OldNulity = CurrentNulity;
      CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);

      loops++;
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
   Mode_->ArcLenUpdate(-(OriginalDiff - IntermediateDiff));
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   return 1;
}
