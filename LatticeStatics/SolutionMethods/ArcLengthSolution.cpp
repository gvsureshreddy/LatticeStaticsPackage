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
      Mode_->ArcLenDef() - two);
}

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,
				     char *startfile,fstream &out)
   : Mode_(Mode), CurrentSolution_(0)
{
   FILE *pipe;
   
   GetParameter("^ArcLenMaxIterations",datafile,"%u",&MaxIter_);
   GetParameter("^ArcLenTolerance",datafile,"%lf",&Tolerance_);
   GetParameter("^ArcLenDSMax",datafile,"%lf",&DSMax_);
   CurrentDS_ = DSMax_;
   GetParameter("^ArcLenDSMin",datafile,"%lf",&DSMin_);
   GetParameter("^ArcLenAngleCutoff",datafile,"%lf",&AngleCutoff_);
   GetParameter("^ArcLenAngleIncrease",datafile,"%lf",&AngleIncrease_);
   GetParameter("^ArcLenAspect",datafile,"%lf",&Aspect_);
   GetParameter("^ArcLenNumSolutions",datafile,"%u",&NumSolutions_);


   const char *StartType[] = {"Bifurcation","Continuation","ConsistencyCheck"};
   switch (GetStringParameter("^StartType",startfile,StartType,3))
   {
      case 0:
      {
	 // Set Difference and Lattice state
	 double eps;
	 GetParameter("^Epsilon",startfile,"%lf",&eps);

	 char tmp[LINELENGTH];
	 char tang[]="^Tangent";
	 SetPerlCommand(tmp,startfile,tang);
	 pipe=OpenPipe(tmp,"r");
	 Difference_.Resize(Mode_->ArcLenDef().Dim());
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&Difference_[i]);
	    Difference_[i] *= eps;
	    
	 }
	 if (pclose(pipe)) Errfun(tang);

	 char bifpt[]="^BifurcationPoint";
	 SetPerlCommand(tmp,startfile,bifpt);
	 pipe=OpenPipe(tmp,"r");
	 Vector stat(Difference_.Dim());
	 for (int i=0;i<stat.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&stat[i]);
	 }
	 if (pclose(pipe)) Errfun(bifpt);
	 // Set Lattice state to the bifurcation point
	 Mode_->ArcLenUpdate(Mode_->ArcLenDef() - stat);

	 break;
      }
      case 1:
      {
	 // Set Difference_ to   two - one
	 char tmp[LINELENGTH];
	 char sol2[]="^Solution2";
	 SetPerlCommand(tmp,startfile,sol2);
	 pipe=OpenPipe(tmp,"r");
	 Difference_.Resize(Mode_->ArcLenDef().Dim());
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&Difference_[i]);
	 }
	 if (pclose(pipe)) Errfun(sol2);
	 
         // Set Lattice state to Solution2
	 Mode_->ArcLenUpdate(Mode_->ArcLenDef() - Difference_);
	 
	 // Get solution1 and set Difference
	 char sol1[]="^Solution1";
	 double tmpval;
	 SetPerlCommand(tmp,startfile,sol1);
	 pipe=OpenPipe(tmp,"r");
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&tmpval);
	    Difference_[i] -= tmpval;
	 }
	 if (pclose(pipe)) Errfun(sol1);
	 break;
      }
      case 2:
      {
	 int Width;
	 int Dim=Mode_->ArcLenDef().Dim();
	 Matrix
	    Stiff(Dim,Dim),
	    Perturbed(Dim,Dim);
	 Vector
	    Solution1(Dim),
	    Solution2(Dim),
	    pert(Dim,0.0),
	    RHS(Dim);
	 
	 Difference_.Resize(Dim);
	 char tmp[LINELENGTH];
	 char sol2[]="^Solution2";
	 SetPerlCommand(tmp,startfile,sol2);
	 pipe=OpenPipe(tmp,"r");
	 for (int i=0;i<Dim;i++)
	 {
	    fscanf(pipe,"%lf",&Solution2[i]);
	 }
	 if (pclose(pipe)) Errfun(sol2);

	 // Get solution1
	 char sol1[]="^Solution1";
	 SetPerlCommand(tmp,startfile,sol1);
	 pipe=OpenPipe(tmp,"r");
	 for (int i=0;i<Dim;i++)
	 {
	    fscanf(pipe,"%lf",&Solution1[i]);
	 }
	 if (pclose(pipe)) Errfun(sol1);

	 // Set Lattice state to Solution2
	 Mode_->ArcLenUpdate(Mode_->ArcLenDef() - Solution2);
	 // Set Difference to Solution2 - Solution1
	 Difference_ = Solution2 - Solution1;

	 // Get Epsilon and Width
	 GetParameter("^ConsistencyEpsilon",startfile,"%lf",&ConsistencyEpsilon_);
	 GetParameter("^MainFieldWidth",datafile,"%i",&Width);

	 // Do Consistency check
	 for (int i=0;i<70;i++) cout << "="; cout << endl;
	 for (int i=0;i<70;i++) out << "="; out << endl;
	 cout << "Consistency Check." << endl
	      << "K(U + DeltaU) * Epsilon" << endl;
	 out << "Consistency Check." << endl
	     << "K(U + DeltaU) * Epsilon" << endl;
	 Mode_->ArcLenUpdate(-Difference_);
	 Stiff = ConsistencyEpsilon_*Mode_->ArcLenStiffness(Difference_,1.0);
	 cout << setw(Width) << Stiff << endl;
	 out << setw(Width) << Stiff << endl;
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    // Get RHS
	    Difference_ = Solution2 - Solution1;
	    Mode_->ArcLenUpdate(Mode_->ArcLenDef() - (Solution2 + Difference_));
	    RHS = Mode_->ArcLenRHS(ConsistencyEpsilon_,Difference_,1.0);
	    
	    // Perturb the lattice state
	    pert=Vector(pert.Dim(),0.0);
	    pert[i]=1.0;
	    Difference_ = Solution2 - Solution1 + ConsistencyEpsilon_*pert;
	    Mode_->ArcLenUpdate(Mode_->ArcLenDef() - (Solution2 + Difference_));
	    // Get Check
	    RHS = RHS - Mode_->ArcLenRHS(ConsistencyEpsilon_,Difference_,1.0);
	    for (int j=0;j<Dim;j++)
	       Perturbed[j][i] = -RHS[j];
	 }

	 // Print out the facts
	 cout << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << endl;
	 cout << setw(Width) << Perturbed << endl;
	 cout << "Difference" << endl;
	 cout << setw(Width) << Stiff - Perturbed << endl;

	 out << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << endl;
	 out << setw(Width) << Perturbed << endl;
	 out << "Difference" << endl;
	 out << setw(Width) << Stiff - Perturbed << endl;
	 
	 for (int i=0;i<70;i++) cout << "="; cout << endl;
	 for (int i=0;i<70;i++) out << "="; out << endl;

	 // We are done -- set currentsolution to numsolutions
	 CurrentSolution_ = NumSolutions_;
      }
   }
}   

int ArcLengthSolution::AllSolutionsFound()
{
   return (CurrentSolution_ == NumSolutions_);
}

double ArcLengthSolution::FindNextSolution(int &good)
{
   double uncertainty;
   double AngleTest;

   Vector OldDiff = Difference_;

   do
   {
      cout << "DS= " << CurrentDS_ << endl;

      uncertainty = ArcLengthNewton(good);

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

   good = 1;
   
   return uncertainty;
}

double ArcLengthSolution::ArcLengthNewton(int &good)
{
   int loc_good = 1;
   double uncertainty;
   
   int itr = 0;
   int Dim=Mode_->ArcLenDef().Dim();

   Vector Dx(Dim),
      RHS(Dim);
   
   Mode_->ArcLenUpdate(-Difference_);

   // Iterate until convergence
   cout << setiosflags(ios::scientific)
	<< "ArcLenNewton: Number of Iterations --\n";

   RHS = Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_);
   do
   {
      itr++;

#ifdef SOLVE_SVD
      Dx = SolveSVD(
	 Mode_->ArcLenStiffness(Difference_,Aspect_),
	 RHS,MAXCONDITION,1);
#else
      Dx = SolvePLU(Mode_->ArcLenStiffness(Difference_,Aspect_),RHS);
#endif

      Mode_->ArcLenUpdate(Dx);
      Difference_ -= Dx;
      RHS = Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_);

      cout << itr << "(" << setw(20)
	   << Mode_->ScanningStressParameter() << ","
	   << setw(20) << RHS.Norm() << ","
	   << setw(20) << Dx.Norm() << "), ";
#ifndef SOLVE_SVD
      cout << endl;
#endif
   }
   while ((itr < MaxIter_)
	  && ((fabs(RHS.Norm()) > Tolerance_) || (fabs(Dx.Norm()) > Tolerance_)));

   cout << resetiosflags(ios::scientific) << endl;
   uncertainty = Dx.Norm();

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << endl;
      loc_good = 0;
   }

   good = good && loc_good;
   return uncertainty;
}

int ArcLengthSolution::BisectAlert(Lattice *Lat,int Width,fstream &out)
{
   double NewtonTolFactor = 10.0,
      ConvergenceFactor = 100;
   
   Vector OriginalDiff=Difference_;
   Vector IntermediateDiff(Difference_.Dim(),0.0);
   double OriginalDS = CurrentDS_;
   double CurrentMinEV,OldMinEV;
   double uncertainty;
   int dummy;
   int loops = 0;
   int OldNulity = Lat->StiffnessNulity(&OldMinEV);
   int CurrentNulity;

   // Set Lattice back to previous solution
   Mode_->ArcLenUpdate(Difference_);
   CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);

   // Set Tolerance_ tighter
   Tolerance_ /= NewtonTolFactor;
   
   cout << "\t" << setw(Width) << OldNulity << setw(Width) << OldMinEV
	<< " DS " << setw(Width) << CurrentDS_ << endl;

   while ((fabs(CurrentMinEV) > ConvergenceFactor*Tolerance_)
	  && (loops < MaxIter_))
   {
      cout << setw(Width) << CurrentNulity
	   << setw(Width) << CurrentMinEV
	   << " DS " << setw(Width) << CurrentDS_ << endl;

      //CurrentDS_ /= 2.0; // Bisection Method
      CurrentDS_ /= (1.0 - (OldMinEV/CurrentMinEV)); // Secant Method
      if (((OldNulity - CurrentNulity) != 0)
	  && (loops != 0))
	 Difference_ = -Difference_;

      uncertainty = ArcLengthNewton(dummy);
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

   cout << setw(Width) << Lat
	<< "Uncertainty = " << setw(Width) << uncertainty << endl
	<< "Success = 1" << endl;
   out << setw(Width) << Lat
       << "Uncertainty = " << setw(Width) << uncertainty << endl
       << "Success = 1" << endl;
      
   for (int i=0;i<70;i++)
   {
      cout << "=";
      out << "=";
   }
   cout << endl; out << endl;

   // Call Lattice function to do any Lattice Specific things
   Lat->CriticalPointInfo(Width,out);
   
   // Reset Lattice and ArcLengthSolution
   Mode_->ArcLenUpdate(-(OriginalDiff - IntermediateDiff));
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;

   // Reste Tolerance_
   Tolerance_ *= NewtonTolFactor;

   return 1;
}
