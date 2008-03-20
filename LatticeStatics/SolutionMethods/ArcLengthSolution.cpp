#include <math.h>
#include "ArcLengthSolution.h"
#include "UtilityFunctions.h"
#include <sstream>
#include <string>

using namespace std;

#define CLOSEDDEFAULT 30
#define ARCLENEPS 1.0e-15

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
                                     const Vector &one,const Vector &two,int Echo)
   : Echo_(Echo),
     Mode_(Mode),
     CurrentSolution_(0),
     Difference_(two-one)
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
   :  Echo_(Echo),
      Mode_(Mode),
      CurrentSolution_(0)
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
   for (unsigned i=0;i<Difference_.Dim();i++)
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

int ArcLengthSolution::FindNextSolution()
{
   int good=0;
   double AngleTest;
   // Assume that first solution should not be strictly restricted to the
   // adaptive steping angle constraint
   double AngleFactor = CurrentSolution_ ? 1.0 : 5.0;
   
   Vector OldDiff = Difference_;
   
   do
   {
      if (Echo_) cout << "DS= " << CurrentDS_ << endl;
      
      ArcLengthNewton(good);
      
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
   
   return good;
}

void ArcLengthSolution::ArcLengthNewton(int &good)
{
   unsigned itr = 0;
   unsigned Dim=ArcLenDef().Dim();
   
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
   
   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << endl;
      good = 0;
   }
   else
   {
      good = 1;
   }
}

int ArcLengthSolution::OldFindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,
                                            Lattice *Lat,char *datafile,const char *prefix,
                                            int Width,fstream &out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   double CurrentMinEV=1.0, OldMinEV=1.0;
   double Delta_DS=0.0;
   int dummy = 1;
   unsigned loops = 0;
   int RighthandTestValue = RHN;
   int CurrentTestValue = RHN;
   int LefthandTestValue = LHN;
   Vector EigenValues(Lat->DOF().Dim());
   
   Delta_DS = CurrentDS_;
   
   Vector Original_DOF((Mode_->ModeDOF()).Dim(),0);
   Original_DOF = Mode_->ModeDOF();
   
   if (Echo_) cout << "LHN = " << LHN << endl << "RHN = " << RHN << endl;
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > BisectTolerance_)
           || (CurrentTestValue == RighthandTestValue))
          && (loops < MaxIter_))
   {
      if (Echo_) cout << "OldMinEV = " << OldMinEV << endl
                      << "CurrentTestValue = " << CurrentTestValue << endl
                      << "CurrentMinEV = "  << CurrentMinEV << endl;
      
      //set to left hand point
      ArcLenUpdate(-Difference_);
      Delta_DS = Delta_DS/2.0;
      
      if(CurrentTestValue == RighthandTestValue)
      {
         CurrentDS_ = CurrentDS_ - Delta_DS;
         Difference_ = Difference_/2.0;
      }
      if(CurrentTestValue == LefthandTestValue)
      {
         CurrentDS_ = CurrentDS_ + Delta_DS;
         Difference_ = 1.5 * Difference_;
      }
      
      cout << "Current_DS = " << CurrentDS_ << endl << endl;
      
      ArcLengthNewton(dummy);
      
      OldMinEV = CurrentMinEV;
      CurrentTestValue = Lat->TestFunctions(EigenValues);
      CurrentMinEV = EigenValues[1]; //not correct;
      
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
   
   if (Echo_) cout << setw(Width) << Lat << endl;
   out << setw(Width) << Lat << endl;
      
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << endl; out << endl;
   
   // Call Lattice function to do any Lattice Specific things
   //  abs(RighthandNulity - LefthandNulity) is the number of zero eigenvalues
   //  in a perfect situation. should check to see if this is found to be true.
   Lat->CriticalPointInfo(Mode_->DrDt(Difference_),abs(RighthandTestValue-LefthandTestValue),
                          BisectTolerance_,datafile,prefix,Width,out);
   
   if (Echo_) cout << "Success = 1" << endl;
   out << "Success = 1" << endl;
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}

int ArcLengthSolution::FindCriticalPoint(Lattice *Lat,char *datafile,const char *prefix,
                                         int Width,fstream &out)
{
   Vector OriginalDiff=Difference_;
   double OriginalDS = CurrentDS_;
   int TestValueDiff;
   int temp;
   int size = Lat->DOF().Dim();
   static Vector TF_LHS(size);
   static Vector TF_RHS(size);
   static Vector CurrentTF(size);
   double fa,fb;
   int Multiplicity;
   int track;
   int num;
   int CP;
   int spot;
   ostringstream in_string;

   // Setup in_string ios
   in_string << setiosflags(ios::fixed) << setprecision(out.precision());
   
   TestValueDiff = Lat->TestFunctions(TF_LHS, Lattice::RHS, &TF_RHS);
   if (TestValueDiff < 0)
   {
      out << "Note: TestFunctions found a discrepancy between the\n"
          << "Note: difference in number of negative Test Functions\n"
          << "Note: and the number of Test Functions that change sign\n"
          << "Note: from LeftHandSide to RightHandSide.  This is usually\n"
          << "Note: caused by having too large of a path-following stepsize."
          << endl;
      if (Echo_)
         cout << "Note: TestFunctions found a discrepancy between the\n"
              << "Note: difference in number of negative Test Functions\n"
              << "Note: and the number of Test Functions that change sign\n"
              << "Note: from LeftHandSide to RightHandSide.  This is usually\n"
              << "Note: caused by having too large of a path-following stepsize."
              << endl;
      TestValueDiff = -TestValueDiff;
   }
   
   cout << "TF_LHS = " << setw(15) << TF_LHS<< endl;
   cout << "TF_RHS = " << setw(15) << TF_RHS << endl;
   
   int *Index;
   Index = new int[TestValueDiff];
   Vector DSTrack(TestValueDiff);
   string *out_string;
   out_string = new string[TestValueDiff];
   
   temp = 0;
   for (int i = 0; i< size; i++)
   {
      if ((TF_LHS[i]*TF_RHS[i]) < 0.0)
      {
         Index[temp] = i;
         temp++;
      }
   }
   
   cout << "TestValueDiff = "<< TestValueDiff << endl;
   if (Echo_)
   {
      for (int i = 0; i < TestValueDiff; i++)
      {
         cout << "i = " << i<< endl;
         cout << "Index[i] =" << Index[i]<< endl;
      }
   }
   
   num = 0;
   for (CP= 0; CP < TestValueDiff; CP++)
   {
      track = Index[CP];
      fa = TF_LHS[track];
      fb = TF_RHS[track];
      
      if(track>=0) //START OF IF STATEMENT
      {
         ZBrent(Lat, track,fa, fb, OriginalDiff, OriginalDS, CurrentTF);
         Multiplicity = 1;
         for(int i=CP+1;i<TestValueDiff;i++)
         {
            temp = Index[i];
            if(fabs(CurrentTF[temp]) <= BisectTolerance_)
            {
               Index[i] = -1;
               Multiplicity++;
               //cout <<"i = " << i << endl <<  "CHECK POINT INDEX[CP] = " << Index[i] << endl;
            }
         }
         
         // sort the critical points
         spot=num;
         while((spot!=0)&&(DSTrack[spot-1] > CurrentDS_) )
         {
            DSTrack[spot] = DSTrack[spot-1];
            out_string[spot] = out_string[spot - 1];
            spot = spot - 1;
         }
         
         
         // Output Critical Point
         for (int i=0;i<70;i++)
         {
            if (Echo_) cout << "=";
            in_string << "=";
         }
         if (Echo_) cout << endl;
         in_string << endl;
         
         if (Echo_) cout << setw(Width) << Lat << endl;
         in_string << setw(Width) << Lat << endl;         
         
         for (int i=0;i<70;i++)
         {
            if (Echo_) cout << "=";
            in_string << "=";
         }
         if (Echo_) cout << endl;
         in_string << endl;
         
         // Call Lattice function to do any Lattice Specific things
         Lat->CriticalPointInfo(Mode_->DrDt(Difference_),Multiplicity,
                                BisectTolerance_,datafile,prefix,Width,in_string);
         
         if (Echo_) cout << "Success = 1" << endl;
         in_string << "Success = 1" << endl;
         
         DSTrack[spot] = CurrentDS_;
         out_string[spot] = in_string.str();
         num = num + 1;
         in_string.str("");
      }//END OF IF STATEMENT
   }
   
   ////PRINT OUT CP DATA
   for (int i = 0; i < num; i++)
   {
      out <<out_string[i];
   }
   
   delete [] Index;
   delete [] out_string;   //deletes memory allocated to out_string
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}

void ArcLengthSolution::ZBrent(Lattice *Lat,int track,double fa,double fb,
                               const Vector &OriginalDiff,const double OriginalDS,
                               Vector &CurrentTF)
{
   Vector LastDiff(Difference_.Dim(),0.0);
   double LastDS=CurrentDS_;
   double a,b,c,d,e,xm,p,fc, tol1,s,q,r,min1,min2;
   int dummy = 1;
   unsigned loops = 0;
   double factor = 0.0;
   
   b=OriginalDS;
   c=b;
   a=0.0;
   
   fc = fb;
   //cout << " a = " << a << endl << "fa = " << fa << endl
   //<<  "b = " << b << endl << "fb = " << fb << endl
   //<< "c = " << c << endl << "fc = " << fc << endl;
   
   while (((fabs(fb) > BisectTolerance_) )&& (loops < MaxIter_))
   {
      if (Echo_)
      {
         cout << setprecision(30) << "CurrentMinTF = " << fb << endl;
         cout << "CurrentDS_ =" << CurrentDS_ << setprecision(10) << endl;
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
         Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
         fb=CurrentTF[track];
         if(Echo_)
         {
            cout <<setprecision(30)<<"CurrentMinTF = " << fb << endl;
            cout << "CurrentDS_ =" << CurrentDS_ <<setprecision(10)<< endl;
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
         }
      }
      else
      {
         d=xm;
         e=d;
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
      ArcLengthNewton(dummy);
      Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
      
      fb=CurrentTF[track];
      loops++;
      
      //cout << "CurrentTF = " << endl << setw(15) << CurrentTF << endl << endl;
      //cout << "CurrentTF[track] = " << endl << CurrentTF[track]<< endl<< endl;
   }
}
