#include "ArcLengthSolution.h"


#include "UtilityFunctions.h"

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				     const Vector &one,const Vector &two,int Echo)
   : Mode_(Mode),Difference_(two-one), CurrentSolution_(0), Echo_(Echo)
{
   if(!GetParameter(prefix,"ArcLenMaxIterations",datafile,"%u",&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenTolerance",datafile,"%lf",&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMax",datafile,"%lf",&DSMax_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSStart",datafile,"%lf",&CurrentDS_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMin",datafile,"%lf",&DSMin_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleCutoff",datafile,"%lf",&AngleCutoff_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleIncrease",datafile,"%lf",&AngleIncrease_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAspect",datafile,"%lf",&Aspect_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenNumSolutions",datafile,"%u",&NumSolutions_)) exit(-1);

   // Set Lattice to solution "two"
   Mode_->ArcLenUpdate(
      Mode_->ArcLenDef() - two);
}

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				     char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   FILE *pipe;
   
   if(!GetParameter(prefix,"ArcLenMaxIterations",datafile,"%u",&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenTolerance",datafile,"%lf",&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMax",datafile,"%lf",&DSMax_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSStart",datafile,"%lf",&CurrentDS_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMin",datafile,"%lf",&DSMin_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleCutoff",datafile,"%lf",&AngleCutoff_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleIncrease",datafile,"%lf",&AngleIncrease_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAspect",datafile,"%lf",&Aspect_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenNumSolutions",datafile,"%u",&NumSolutions_)) exit(-1);


   const char *StartType[] = {"Bifurcation","Continuation","ConsistencyCheck"};
   switch (GetStringParameter(prefix,"StartType",startfile,StartType,3))
   {
      case -1:
      {
	 cerr << "Unknown StartType!" << endl;
	 exit(-1);
	 break;
      }
      case 0:
      {
	 // Set Difference and Lattice state
	 double eps;
	 if(!GetParameter(prefix,"Epsilon",startfile,"%lf",&eps)) exit(-1);

	 char tmp[LINELENGTH];
	 char tang[]="Tangent";
	 SetPerlCommand(tmp,startfile,prefix,tang);
	 pipe=OpenPipe(tmp,"r");
	 Difference_.Resize(Mode_->ArcLenDef().Dim());
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&Difference_[i]);
	    Difference_[i] *= eps;
	    
	 }
	 if (pclose(pipe)) Errfun(tang);

	 char bifpt[]="BifurcationPoint";
	 SetPerlCommand(tmp,startfile,prefix,bifpt);
	 pipe=OpenPipe(tmp,"r");
	 Vector stat(Difference_.Dim());
	 for (int i=0;i<stat.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&stat[i]);
	 }
	 if (pclose(pipe)) Errfun(bifpt);
	 // Set Lattice state to the bifurcation point
	 Mode_->ArcLenSet(stat);

	 break;
      }
      case 1:
      {
	 // Set Difference_ to   two - one
	 char tmp[LINELENGTH];
	 char sol2[]="Solution2";
	 SetPerlCommand(tmp,startfile,prefix,sol2);
	 pipe=OpenPipe(tmp,"r");
	 Difference_.Resize(Mode_->ArcLenDef().Dim());
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    fscanf(pipe,"%lf",&Difference_[i]);
	 }
	 if (pclose(pipe)) Errfun(sol2);
	 
         // Set Lattice state to Solution2
	 Mode_->ArcLenSet(Difference_);
	 
	 // Get solution1 and set Difference
	 char sol1[]="Solution1";
	 double tmpval;
	 SetPerlCommand(tmp,startfile,prefix,sol1);
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
	 double potential;
	 int Width;
	 int Dim=Mode_->ArcLenDef().Dim();
	 Matrix
	    Stiff(Dim,Dim),
	    PerturbedStiff(Dim,Dim);
	 Vector
	    Force(Dim),
	    PerturbedForce(Dim),
	    Solution1(Dim),
	    Solution2(Dim),
	    pert(Dim,0.0),
	    RHS(Dim);
	 
	 Difference_.Resize(Dim);
	 char tmp[LINELENGTH];
	 char sol2[]="Solution2";
	 SetPerlCommand(tmp,startfile,prefix,sol2);
	 pipe=OpenPipe(tmp,"r");
	 for (int i=0;i<Dim;i++)
	 {
	    fscanf(pipe,"%lf",&Solution2[i]);
	 }
	 if (pclose(pipe)) Errfun(sol2);

	 // Get solution1
	 char sol1[]="Solution1";
	 SetPerlCommand(tmp,startfile,prefix,sol1);
	 pipe=OpenPipe(tmp,"r");
	 for (int i=0;i<Dim;i++)
	 {
	    fscanf(pipe,"%lf",&Solution1[i]);
	 }
	 if (pclose(pipe)) Errfun(sol1);

	 // Set Lattice state to Solution2
	 Mode_->ArcLenSet(Solution2);
	 // Set Difference to Solution2 - Solution1
	 Difference_ = Solution2 - Solution1;
	 
	 // Get Epsilon and Width
	 if(!GetParameter(prefix,"ConsistencyEpsilon",startfile,"%lf",&ConsistencyEpsilon_)) exit(-1);
	 if(!GetParameter(prefix,"MainFieldWidth",datafile,"%i",&Width)) exit(-1);

	 // Do Consistency check
	 if (Echo_)
	 {
	    for (int i=0;i<70;i++) cout << "="; cout << endl;
	    cout << "Consistency Check." << endl;
	    cout << "F(U + DeltaU) * Epsilon" << endl;
	 }
	 for (int i=0;i<70;i++) out << "="; out << endl;
	 out << "Consistency Check." << endl;
	 out << "F(U + DeltaU) * Epsilon" << endl;
	 Mode_->ArcLenUpdate(-Difference_);
	 Force = ConsistencyEpsilon_*Mode_->ArcLenRHS(ConsistencyEpsilon_,Difference_,1.0);
	 if (Echo_)
	 {
	    cout << setw(Width) << Force << endl;
	    cout << "K(U + DeltaU) * Epsilon" << endl;
	 }
	 out << setw(Width) << Force << endl;
	 out << "K(U + DeltaU) * Epsilon" << endl;
	 Stiff = ConsistencyEpsilon_*Mode_->ArcLenStiffness(Difference_,1.0);
	 if (Echo_) cout << setw(Width) << Stiff << endl;
	 out << setw(Width) << Stiff << endl;
	 for (int i=0;i<Difference_.Dim();i++)
	 {
	    // Get RHS
	    Difference_ = Solution2 - Solution1;
	    Mode_->ArcLenUpdate(Mode_->ArcLenDef() - (Solution2 + Difference_));
	    potential = Mode_->ModeEnergy();
	    RHS = Mode_->ArcLenRHS(ConsistencyEpsilon_,Difference_,1.0);
	    
	    // Perturb the lattice state
	    pert=Vector(pert.Dim(),0.0);
	    pert[i]=1.0;
	    Difference_ = Solution2 - Solution1 + ConsistencyEpsilon_*pert;
	    Mode_->ArcLenUpdate(Mode_->ArcLenDef() - (Solution2 + Difference_));
	    // Get Check
	    potential = potential - Mode_->ModeEnergy();
	    PerturbedForce[i] = -potential;
	    RHS = RHS - Mode_->ArcLenRHS(ConsistencyEpsilon_,Difference_,1.0);
	    for (int j=0;j<Dim;j++)
	       PerturbedStiff[j][i] = -RHS[j];
	 }

	 // Print out the facts
	 if (Echo_)
	 {
	    cout << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << endl;
	    cout << setw(Width) << PerturbedForce << endl;
	    cout << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << endl;
	    cout << setw(Width) << PerturbedStiff << endl;
	    cout << "Difference" << endl;
	    cout << setw(Width) << Force - PerturbedForce << endl << endl;
	    cout << setw(Width) << Stiff - PerturbedStiff << endl;
	 }

	 out << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << endl;
	 out << setw(Width) << PerturbedForce << endl;
	 out << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << endl;
	 out << setw(Width) << PerturbedStiff << endl;
	 out << "Difference" << endl;
	 out << setw(Width) << Force - PerturbedForce << endl << endl;
	 out << setw(Width) << Stiff - PerturbedStiff << endl;

	 if (Echo_)
	 {
	    for (int i=0;i<70;i++) cout << "="; cout << endl;
	 }
	 for (int i=0;i<70;i++) out << "="; out << endl;

	 // We are done -- set currentsolution to numsolutions
	 CurrentSolution_ = NumSolutions_;
	 break;
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
      if (Echo_) cout << "DS= " << CurrentDS_ << endl;

      uncertainty = ArcLengthNewton(good);

      AngleTest = Mode_->ArcLenAngle(OldDiff,Difference_,Aspect_);

      if (Echo_)
	 cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << endl;
   }
   while (((AngleTest >= AngleCutoff_) || !good)
	  && (CurrentDS_ >= DSMin_)
	  && (Mode_->ArcLenUpdate(Difference_),
	      Difference_ = OldDiff,
	      CurrentDS_=CurrentDS_/2.0,
	      good=1));

   if ((AngleTest <= AngleIncrease_) && (CurrentDS_ < DSMax_))
   {
      CurrentDS_ *= 2.0;
      if (CurrentDS_ > DSMax_) CurrentDS_ = DSMax_;
      if (Echo_) cout << "DS= " << CurrentDS_ << endl;
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
   if (Echo_) cout << setiosflags(ios::scientific)
		  << "ArcLenNewton: Number of Iterations --\n";

   RHS = Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_);
   do
   {
      itr++;

#ifdef SOLVE_SVD
      Dx = SolveSVD(
	 Mode_->ArcLenStiffness(Difference_,Aspect_),
	 RHS,MAXCONDITION,Echo_);
#else
      Dx = SolvePLU(Mode_->ArcLenStiffness(Difference_,Aspect_),RHS);
#endif

      Mode_->ArcLenUpdate(Dx);
      Difference_ -= Dx;
      RHS = Mode_->ArcLenRHS(CurrentDS_,Difference_,Aspect_);

      if (Echo_) cout << itr << "(" << setw(20)
		      << Mode_->ScanningStressParameter() << ","
		      << setw(20) << RHS.Norm() << ","
		      << setw(20) << Dx.Norm() << "), ";
#ifndef SOLVE_SVD
      if (Echo_) cout << endl;
#endif
   }
   while ((itr < MaxIter_)
	  && ((fabs(RHS.Norm()) > Tolerance_) || (fabs(Dx.Norm()) > Tolerance_)));

   if (Echo_) cout << resetiosflags(ios::scientific) << endl;
   uncertainty = Dx.Norm();

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << endl;
      loc_good = 0;
   }

   good = good && loc_good;
   return uncertainty;
}

int ArcLengthSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
				   int Width,fstream &out)
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
   // OriginalNulity is the nulity on the front side of the path being traced
   int OriginalNulity = OldNulity; 
   int CurrentNulity;

   // Set Lattice back to previous solution
   Mode_->ArcLenUpdate(Difference_);
   CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);

   // Set Tolerance_ tighter
   Tolerance_ /= NewtonTolFactor;
   
   if (Echo_) cout << "\t" << setw(Width) << OldNulity << setw(Width) << OldMinEV
		  << " DS " << setw(Width) << CurrentDS_ << endl;
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > ConvergenceFactor*Tolerance_)
	   || (CurrentNulity == OriginalNulity))
	  && (loops < MaxIter_))
   {
      if (Echo_) cout << setw(Width) << CurrentNulity
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
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << endl;
   out << endl;

   if (Echo_) cout << setw(Width) << Lat
		  << "Uncertainty = " << setw(Width) << uncertainty << endl
		  << "Success = 1" << endl;
   out << setw(Width) << Lat
       << "Uncertainty = " << setw(Width) << uncertainty << endl
       << "Success = 1" << endl;
      
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << endl; out << endl;

   // Call Lattice function to do any Lattice Specific things
   Lat->CriticalPointInfo(Mode_->DrDt(Difference_),ConvergenceFactor*Tolerance_,datafile,
			  prefix,Width,out);
   
   // Reset Lattice and ArcLengthSolution
   Mode_->ArcLenUpdate(-(OriginalDiff - IntermediateDiff));
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   // Reset Tolerance_
   Tolerance_ *= NewtonTolFactor;

   return 1;
}
