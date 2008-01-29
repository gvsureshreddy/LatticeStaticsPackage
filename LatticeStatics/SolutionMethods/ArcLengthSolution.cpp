#include <math.h>
#include "ArcLengthSolution.h"
#include "UtilityFunctions.h"

using namespace std;

#define CLOSEDDEFAULT 30
#define ARCLENEPS 1.0e-15

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				     const Vector &one,const Vector &two,int Echo)
   : Mode_(Mode),Difference_(two-one), CurrentSolution_(0), Echo_(Echo)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   if(!GetParameter(prefix,"ArcLenMaxIterations",datafile,'u',&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenTolerance",datafile,'l',&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMax",datafile,'l',&DSMax_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSStart",datafile,'l',&CurrentDS_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMin",datafile,'l',&DSMin_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleCutoff",datafile,'l',&AngleCutoff_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleIncrease",datafile,'l',&AngleIncrease_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAspect",datafile,'l',&Aspect_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenClosedLoopStart",datafile,'i',&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   
   BisectTolerance_ = Tolerance_;
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   
   // Set Lattice to solution "two"
   ArcLenSet(two);
}

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				     char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   if(!GetParameter(prefix,"ArcLenMaxIterations",datafile,'u',&MaxIter_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenTolerance",datafile,'l',&Tolerance_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMax",datafile,'l',&DSMax_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSStart",datafile,'l',&CurrentDS_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenDSMin",datafile,'l',&DSMin_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleCutoff",datafile,'l',&AngleCutoff_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAngleIncrease",datafile,'l',&AngleIncrease_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenAspect",datafile,'l',&Aspect_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
   if(!GetParameter(prefix,"ArcLenClosedLoopStart",datafile,'i',&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   
   BisectTolerance_ = Tolerance_;
   
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
	 if(!GetParameter(prefix,"Epsilon",startfile,'l',&eps)) exit(-1);
	 
	 Difference_.Resize(ArcLenDef().Dim());
	 if(!GetVectorParameter(prefix,"Tangent",startfile,&Difference_)) exit(-1);
	 Difference_ *= eps;
	 
	 Vector stat(Difference_.Dim());
	 if(!GetVectorParameter(prefix,"BifurcationPoint",startfile,&stat)) exit(-1);
	 // Set Lattice state to the bifurcation point
	 ArcLenSet(stat);
	 
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
	 Vector two(ArcLenDef().Dim());
	 if(!GetVectorParameter(prefix,"Solution2",startfile,&two)) exit(-1);
	 ArcLenSet(two);
	 
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
	    Dim=ArcLenDef().Dim();
	 Vector Solution1(Dim),
	    Solution2(Dim);
	 
	 if(!GetVectorParameter(prefix,"Solution2",startfile,&Solution2)) exit(-1);
	 if(!GetVectorParameter(prefix,"Solution1",startfile,&Solution1)) exit(-1);
	 // Get Epsilon and Width
	 if(!GetParameter(prefix,"ConsistencyEpsilon",startfile,'l',&ConsistencyEpsilon)) exit(-1);
	 if(!GetParameter(prefix,"MainFieldWidth",datafile,'i',&Width)) exit(-1);
	 
	 ConsistencyCheck(Solution1,Solution2,ConsistencyEpsilon,Width,out);
	 break;
      }
   }
}

Vector ArcLengthSolution::ArcLenForce(double DS,const Vector &Diff,
				      double Aspect)
{
   static Vector force(ModeDOFS_);
   static Vector mdfc(ModeDOFS_-1);
   mdfc = Mode_->ModeForce();
   
   force[ModeDOFS_-1] = DS*DS - Diff[ModeDOFS_-1]*Diff[ModeDOFS_-1]/(Aspect*Aspect);
   for (int i=0;i<ModeDOFS_-1;++i)
   {
      force[i] = mdfc[i];
      force[ModeDOFS_-1] -= Diff[i]*Diff[i];
   }
   
   return force;
}

Matrix ArcLengthSolution::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   static Matrix K(ModeDOFS_,ModeDOFS_);
   static Matrix ModeK(ModeDOFS_-1,ModeDOFS_);
   
   ModeK = Mode_->ModeStiffness();
   
   for (int i=0;i<ModeDOFS_-1;++i)
   {
      for (int j=0;j<=ModeDOFS_-1;++j)
      {
	 K[i][j] = ModeK[i][j];
      }
      K[ModeDOFS_-1][i] = -2.0*Diff[i];
   }
   K[ModeDOFS_-1][ModeDOFS_-1] = -2.0*Diff[ModeDOFS_-1]/(Aspect*Aspect);
   
   return K;
}

double ArcLengthSolution::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[ModeDOFS_-1] /= Aspect;
   New[ModeDOFS_-1] /= Aspect;
   
   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

void ArcLengthSolution::ConsistencyCheck(Vector &Solution1,Vector &Solution2,
					 double ConsistencyEpsilon,int Width,fstream &out)
{
   double potential;
   int Dim=ArcLenDef().Dim();
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
   ArcLenSet(Solution2);
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
   ArcLenUpdate(Difference_);
   Force = ConsistencyEpsilon*ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
   if (Echo_)
   {
      cout << setw(Width) << Force << endl;
      cout << "K(U + DeltaU) * Epsilon" << endl;
   }
   out << setw(Width) << Force << endl;
   out << "K(U + DeltaU) * Epsilon" << endl;
   Stiff = ConsistencyEpsilon*ArcLenStiffness(Difference_,1.0);
   if (Echo_) cout << setw(Width) << Stiff << endl;
   out << setw(Width) << Stiff << endl;
   for (int i=0;i<Difference_.Dim();i++)
   {
      // Get RHS
      Difference_ = Solution2 - Solution1;
      ArcLenSet(Solution2 + Difference_);
      potential = Mode_->ModeEnergy();
      RHS = ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
      
      // Perturb the lattice state
      pert=Vector(pert.Dim(),0.0);
      pert[i]=1.0;
      Difference_ = Solution2 - Solution1 + ConsistencyEpsilon*pert;
      ArcLenSet(Solution2 + Difference_);
      // Get Check
      potential = Mode_->ModeEnergy() - potential;
      // fix-up the arclength equation part of PerturbedForce
      if (i == RHS.Dim()-1) potential = ConsistencyEpsilon*RHS[i];
      PerturbedForce[i] = potential;
      RHS = ArcLenForce(ConsistencyEpsilon,Difference_,1.0) - RHS;
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
      
      AngleTest = ArcLenAngle(OldDiff,Difference_,Aspect_);
      
      if (Echo_)
	 cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << endl;
   }
   while (((AngleTest >= AngleFactor*AngleCutoff_) || !good)
	  && (CurrentDS_ >= DSMin_)
	  && (ArcLenUpdate(-Difference_),// back to previous solution
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
       ((ArcLenDef() - FirstSolution_).Norm() < CurrentDS_))
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
   int Dim=ArcLenDef().Dim();
   
   Vector Dx(Dim),
      RHS(Dim);
   Matrix stif(Dim,Dim);
   
   // Predictor step
   ArcLenUpdate(Difference_);
   
   // Iterate until convergence
   if (Echo_) cout << setiosflags(ios::scientific)
		   << "ArcLenNewton: Number of Iterations --\n";
   
   RHS = -ArcLenForce(CurrentDS_,Difference_,Aspect_);
   stif=ArcLenStiffness(Difference_,Aspect_);
   
   do
   {
      itr++;
      
#ifdef SOLVE_SVD
      Dx = SolveSVD(
	 stif,
	 RHS,MAXCONDITION,Echo_);
#else
      Dx = SolvePLU(stif,RHS);
#endif
      
      ArcLenUpdate(Dx);
      Difference_ += Dx;
      RHS = -ArcLenForce(CurrentDS_,Difference_,Aspect_);
      stif=ArcLenStiffness(Difference_,Aspect_);
      
      if (Echo_) cout << itr << "("
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

int ArcLengthSolution::OldBisectAlert(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				      char *datafile,const char *prefix,int Width,fstream &out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   double CurrentMinEV=1.0, OldMinEV;
   double uncertainty;
   double Delta_DS=0.0;
   int dummy = 1;
   int loops = 0;
   int RighthandNulity = RHN;
   int CurrentNulity= RHN;
   int LefthandNulity = LHN;
   
   Delta_DS = CurrentDS_;
   
   Vector Original_DOF((Mode_->ModeDOF()).Dim(),0);
   Original_DOF = Mode_->ModeDOF();
   
   if (Echo_) cout << "LHN = " << LHN << endl << "RHN = " << RHN << endl;
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > BisectTolerance_)
	   || (CurrentNulity == RighthandNulity))
	  && (loops < MaxIter_))
   {
      if (Echo_) cout << "OldMinEV = " << OldMinEV << endl
		      << "CurrentNulity = " << CurrentNulity << endl
		      << "CurrentMinEV = "  << CurrentMinEV << endl;
      
      //set to left hand point
      ArcLenUpdate(-Difference_);
      Delta_DS = Delta_DS/2.0;
      
      if(CurrentNulity == RighthandNulity)
      {
	 CurrentDS_ = CurrentDS_ - Delta_DS;
	 Difference_ = Difference_/2.0;
      }
      if(CurrentNulity == LefthandNulity)
      {
	 CurrentDS_ = CurrentDS_ + Delta_DS;
	 Difference_ = 1.5 * Difference_;
      }
      
      cout << "Current_DS = " << CurrentDS_ << endl << endl;
      
      uncertainty = ArcLengthNewton(dummy);
      
      OldMinEV = CurrentMinEV;
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
			  BisectTolerance_,datafile,prefix,Width,out);
   
   if (Echo_) cout << "Success = 1" << endl;
   out << "Success = 1" << endl;
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}

int ArcLengthSolution::BisectAlert(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
				   char *datafile,const char *prefix,int Width,fstream &out)
{
   Vector OriginalDiff=Difference_;
   Vector LastDiff(Difference_.Dim(),0.0);
   double OriginalDS = CurrentDS_;
   double LastDS;
   double uncertainty;
   double a,b,c,d,e,xm,p, fa,fb,fc, tol1,s,q,r,min1,min2;
   int dummy = 1;
   int loops = 0;
   int RighthandNulity = RHN;
   int CurrentNulity = RHN;
   int LefthandNulity = LHN;
   double factor = 0.0;
   
   b=OriginalDS;
   c=b;
   a=0.0;
   fa=LHEV;
   fb=RHEV;
   fc=fb;
   
   //CurrentNulity == RighthandNulity gives some problems
   while (((fabs(fb) > BisectTolerance_) || (CurrentNulity == RighthandNulity))
	  && (loops < MaxIter_))
   {
      if (Echo_)
      {
	 cout << setprecision(30) << "CurrentMinEV = " << fb << endl;
	 cout << "CurrentDS_ =" << CurrentDS_ << setprecision(10) << endl;
	 cout << "LHN = " << LHN << endl << "RHN = " << RHN << endl
	      << "CurrentNulity = " << CurrentNulity << endl << endl;
      }
      ArcLenUpdate(-Difference_);
      
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
      {
	 c=a;
	 fc=fa;
	 e=b-a;
	 d=e;
      }
      
      if (fabs(fc) < fabs(fb))
      {
	 a=b;
	 b=c;
	 c=a;
	 fa=fb;
	 fb=fc;
	 fc=fa;
      }
      
      tol1 = 2.0*ARCLENEPS*fabs(b) + 0.5*BisectTolerance_;
      xm = 0.5 * (c-b);
      
      if(fabs(xm) <= tol1 || fb == 0.0)
      {
	 //cout << "Minimal Root found! XM too small! " < endl;
	 CurrentDS_ = LastDS;
	 Difference_ = LastDiff;
	 ArcLenUpdate(Difference_);
	 CurrentNulity = Lat->StiffnessNulity(&fb);
	 
	 if(Echo_)
	 {
	    cout <<setprecision(30)<<"CurrentMinEV = " << fb << endl;
	    cout << "CurrentDS_ =" << CurrentDS_ <<setprecision(10)<< endl;
	    cout << "LHN = " << LHN << endl << "RHN = " << RHN << endl
		 << "CurrentNulity = " << CurrentNulity << endl << endl;
	 }
	 break;
      }
      
      if((fabs(e) >= tol1) && (fabs(fa) > fabs(fb)))
      {
	 s=fb/fa;
	 if(a == c)
	 {
	    p=2.0*xm*s;
	    q=1.0-s;
	 }
	 else
	 {
	    q=fa/fc;
	    r=fb/fc;
	    p=s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
	    q=(q-1.0)*(r-1.0)*(s-1.0);
	 }
	 
	 if(p > 0.0)
	 {
	    q=-q;
	 }
	 
	 p=fabs(p);
	 min1=3.0*xm*q - fabs(tol1*q);
	 min2=fabs(e*q);
	 
	 if(2.0*p < (min1 < min2 ? min1 : min2))
	 {
	    e=d;
	    d=p/q;
	 }
	 else
	 {
	    d=xm;
	    e=d;
            //cout << "BISECT ALERT POINT 1" << endl;
	 }
      }
      else
      {
	 d=xm;
	 e=d;
	 //cout << "BISECT ALERT POINT 2" << endl;
      }
      
      a=b;
      fa=fb;
      if(fabs(d) > tol1)
      {
	 b += d;
      }
      else
      {
	 b += (xm >=0.0 ? fabs(tol1) : -fabs(tol1));
      }
      LastDiff = Difference_;
      LastDS = CurrentDS_;
      CurrentDS_ = b;
      factor = OriginalDS/b;
      Difference_ = OriginalDiff/factor;
      uncertainty = ArcLengthNewton(dummy);
      CurrentNulity = Lat->StiffnessNulity(&fb);
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
			  BisectTolerance_,datafile,prefix,Width,out);
   
   if (Echo_) cout << "Success = 1" << endl;
   out << "Success = 1" << endl;
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}
