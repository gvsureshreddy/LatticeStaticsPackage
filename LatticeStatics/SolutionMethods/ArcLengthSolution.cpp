#include "ArcLengthSolution.h"


#include "UtilityFunctions.h"

using namespace std;

#define CLOSEDDEFAULT 30

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
   if(!GetParameter(prefix,"ArcLenClosedLoopStart",datafile,"%i",&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   
   // Set Lattice to solution "two"
   Mode_->ArcLenSet(two);
}

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				     char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
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
   if(!GetParameter(prefix,"ArcLenClosedLoopStart",datafile,"%i",&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }

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

	 Difference_.Resize(Mode_->ArcLenDef().Dim());
	 if(!GetVectorParameter(prefix,"Tangent",startfile,&Difference_)) exit(-1);
	 Difference_ *= eps;

	 Vector stat(Difference_.Dim());
	 if(!GetVectorParameter(prefix,"BifurcationPoint",startfile,&stat)) exit(-1);
	 // Set Lattice state to the bifurcation point
	 Mode_->ArcLenSet(stat);

	 // Set FirstSolution
	 FirstSolution_.Resize(stat.Dim());
	 if(!GetVectorParameter(prefix,"ClosedLoopFirstPoint",startfile,&FirstSolution_,0))
	 {
	    FirstSolution_ = stat;
	 }

	 break;
      }
      case 1:
      {
         // Set Lattice state to Solution2
	 Vector two(Mode_->ArcLenDef().Dim());
	 if(!GetVectorParameter(prefix,"Solution2",startfile,&two)) exit(-1);
	 Mode_->ArcLenSet(two);
	 
	 // Get solution1
	 Vector one(two.Dim());
	 if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
	 // Set Difference_ to   two - one
	 Difference_.Resize(two.Dim());
	 Difference_ = two - one;

	 // Set FirstSolution
	 FirstSolution_.Resize(two.Dim());
	 if(!GetVectorParameter(prefix,"ClosedLoopFirstPoint",startfile,&FirstSolution_,0))
	 {
	    FirstSolution_ = two;
	 }
	 
	 
	 break;
      }
      case 2:
      {
	 double ConsistencyEpsilon;
	 int Width,
	    Dim=Mode_->ArcLenDef().Dim();
	 Vector Solution1(Dim),
	    Solution2(Dim);

	 if(!GetVectorParameter(prefix,"Solution2",startfile,&Solution2)) exit(-1);
	 if(!GetVectorParameter(prefix,"Solution1",startfile,&Solution1)) exit(-1);
	 // Get Epsilon and Width
	 if(!GetParameter(prefix,"ConsistencyEpsilon",startfile,"%lf",&ConsistencyEpsilon)) exit(-1);
	 if(!GetParameter(prefix,"MainFieldWidth",datafile,"%i",&Width)) exit(-1);


	 ConsistencyCheck(Solution1,Solution2,ConsistencyEpsilon,Width,out);
	 break;
      }
   }
}


void ArcLengthSolution::ConsistencyCheck(Vector &Solution1,Vector &Solution2,
					 double ConsistencyEpsilon,int Width,fstream &out)
{
   double potential;
   int Dim=Mode_->ArcLenDef().Dim();
   Matrix
      Stiff(Dim,Dim),
      PerturbedStiff(Dim,Dim);
   Vector
      Force(Dim),
      PerturbedForce(Dim),
      pert(Dim,0.0),
      RHS(Dim);
   
   Difference_.Resize(Dim);
   
   // Set Lattice state to Solution2
   Mode_->ArcLenSet(Solution2);
   // Set Difference to Solution2 - Solution1
   Difference_ = Solution2 - Solution1;
	 
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
   Mode_->ArcLenUpdate(Difference_);
   Force = ConsistencyEpsilon*Mode_->ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
   if (Echo_)
   {
      cout << setw(Width) << Force << endl;
      cout << "K(U + DeltaU) * Epsilon" << endl;
   }
   out << setw(Width) << Force << endl;
   out << "K(U + DeltaU) * Epsilon" << endl;
   Stiff = ConsistencyEpsilon*Mode_->ArcLenStiffness(Difference_,1.0);
   if (Echo_) cout << setw(Width) << Stiff << endl;
   out << setw(Width) << Stiff << endl;
   for (int i=0;i<Difference_.Dim();i++)
   {
      // Get RHS
      Difference_ = Solution2 - Solution1;
      Mode_->ArcLenSet(Solution2 + Difference_);
      potential = Mode_->ModeEnergy();
      RHS = Mode_->ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
	    
      // Perturb the lattice state
      pert=Vector(pert.Dim(),0.0);
      pert[i]=1.0;
      Difference_ = Solution2 - Solution1 + ConsistencyEpsilon*pert;
      Mode_->ArcLenSet(Solution2 + Difference_);
      // Get Check
      potential = Mode_->ModeEnergy() - potential;
      // fix-up the arclength equation part of PerturbedForce
      if (i == RHS.Dim()-1) potential = ConsistencyEpsilon*RHS[i];
      PerturbedForce[i] = potential;
      RHS = Mode_->ArcLenForce(ConsistencyEpsilon,Difference_,1.0) - RHS;
      for (int j=0;j<Dim;j++)
	 PerturbedStiff[j][i] = RHS[j];
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
}

int ArcLengthSolution::AllSolutionsFound()
{
   return (CurrentSolution_ >= NumSolutions_);
}

double ArcLengthSolution::FindNextSolution(int &good)
{
   double uncertainty;
   double AngleTest;
   // Assume that first solution should not be strictly restricted to the
   // adaptive steping angle constraint
   double AngleFactor = CurrentSolution_ ? 1.0 : 5.0;

   Vector OldDiff = Difference_;

   do
   {
      if (Echo_) cout << "DS= " << CurrentDS_ << endl;

      uncertainty = ArcLengthNewton(good);

      AngleTest = Mode_->ArcLenAngle(OldDiff,Difference_,Aspect_);

      if (Echo_)
	 cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << endl;
   }
   while (((AngleTest >= AngleFactor*AngleCutoff_) || !good)
	  && (CurrentDS_ >= DSMin_)
	  && (Mode_->ArcLenUpdate(-Difference_),// back to previous solution
	      Difference_ = OldDiff,
	      CurrentDS_=CurrentDS_/2.0));

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

   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Mode_->ArcLenDef() - FirstSolution_).Norm() < CurrentDS_))
   {
      // We are done -- set currentsolution to numsolutions
      cerr << "Closed Loop detected at Solution # " << CurrentSolution_
	   << " --- Terminating!" << endl;

      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      CurrentSolution_++;
   }

   // Always have the current "solution" state printed as a solution point
   good = 1;
   
   return uncertainty;
}

double ArcLengthSolution::ArcLengthNewton(int &good)
{
   double uncertainty;
   
   int itr = 0;
   int Dim=Mode_->ArcLenDef().Dim();

   Vector Dx(Dim),
      RHS(Dim);

   // Predictor step
   Mode_->ArcLenUpdate(Difference_);

   // Iterate until convergence
   if (Echo_) cout << setiosflags(ios::scientific)
		  << "ArcLenNewton: Number of Iterations --\n";

   RHS = -Mode_->ArcLenForce(CurrentDS_,Difference_,Aspect_);
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
      Difference_ += Dx;
      RHS = -Mode_->ArcLenForce(CurrentDS_,Difference_,Aspect_);

      if (Echo_) cout << itr << "(" << setw(20)
		      << Mode_->ScanningStressParameter() << ","
		      << setw(20) << RHS.Norm() << ","
		      << setw(20) << Dx.Norm()  << "), ";
#ifdef SOLVE_SVD
      if (Echo_) cout << endl;
#endif
   }
   while ((itr < MaxIter_) && ((RHS.Norm() > Tolerance_) || (Dx.Norm() > Tolerance_)));

   if (Echo_) cout << resetiosflags(ios::scientific) << endl;
   uncertainty = Dx.Norm();

   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << endl;
      good = 0;
   }
   else
   {
      good = 1;
   }

   return uncertainty;
}

int ArcLengthSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
				   int Width,fstream &out)
{
   double BisectTolerance;
   if(!GetParameter(prefix,"ArcLenBisectTolerance",datafile,"%lf",&BisectTolerance))
   {
      // Default to 10*Tolerance_
      BisectTolerance = 10*Tolerance_;
   }
   
   Vector OriginalDiff=Difference_;
   Vector IntermediateDiff(Difference_.Dim(),0.0);
   double OriginalDS = CurrentDS_;
   double CurrentMinEV,OldMinEV;
   double uncertainty;
   int dummy = 1;
   int loops = 0;
   int OldNulity = Lat->StiffnessNulity(&OldMinEV);
   // RigthhandNulity is the nulity on the front side of the path being traced
   int RighthandNulity = OldNulity; 
   int CurrentNulity;

   // Set Lattice back to previous solution
   Mode_->ArcLenUpdate(-Difference_);
   CurrentNulity = Lat->StiffnessNulity(&CurrentMinEV);
   // LefthandNulity is the nulity on the back side of the path being traced
   int LefthandNulity = CurrentNulity;

   if (Echo_) cout << "\t" << setw(Width) << OldNulity << setw(Width) << OldMinEV
		  << " DS " << setw(Width) << CurrentDS_ << endl;
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > BisectTolerance)
	   || (CurrentNulity == RighthandNulity))
	  && (loops < MaxIter_))
   {
      if (Echo_) cout << setw(Width) << CurrentNulity
		     << setw(Width) << CurrentMinEV
		     << " DS " << setw(Width) << CurrentDS_ << endl;

      // stick with bisection --- secant method proves problematic.
      CurrentDS_ /= 2.0; // Bisection Method
      // this takes care of the direction that needs to be searched (thus the positive DS)
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
		   << "Uncertainty = " << setw(Width) << uncertainty << endl;
   out << setw(Width) << Lat
       << "Uncertainty = " << setw(Width) << uncertainty << endl;
      
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << endl; out << endl;

   // Call Lattice function to do any Lattice Specific things
   //  abs(RighthandNulity - LefthandNulity) is the number of zero eigenvalues
   //  in a perfect situation. should check to see if this is found to be true.
   Lat->CriticalPointInfo(Mode_->DrDt(Difference_),abs(RighthandNulity-LefthandNulity),
			  BisectTolerance,datafile,prefix,Width,out);

   if (Echo_) cout << "Success = 1" << endl;
   out << "Success = 1" << endl;
   
   // Reset Lattice and ArcLengthSolution
   Mode_->ArcLenUpdate(OriginalDiff - IntermediateDiff);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}
