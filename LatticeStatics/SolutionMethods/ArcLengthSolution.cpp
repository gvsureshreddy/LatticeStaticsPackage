#include <math.h>
#include <sstream>
#include <string>
#include "ArcLengthSolution.h"

using namespace std;

#define ARCLENEPS 1.0e-15

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,const Vector &dofs,
                                     unsigned MaxIter,double Tolerance,
                                     double BisectTolerance,double DSMax,double DSMin,
                                     double CurrentDS,double AngleCutoff,double AngleIncrease,
                                     double Aspect,unsigned NumSolutions,
                                     unsigned CurrentSolution,const Vector &FirstSolution,
                                     const Vector &Difference,unsigned ClosedLoopStart,int Echo)
   : Echo_(Echo),
     Mode_(Mode),
     ModeDOFS_(Mode_->ModeDOF().Dim()),
     MaxIter_(MaxIter),
     Tolerance_(Tolerance),
     BisectTolerance_(BisectTolerance),
     DSMax_(DSMax),
     DSMin_(DSMin),
     CurrentDS_(CurrentDS),
     AngleCutoff_(AngleCutoff),
     AngleIncrease_(AngleIncrease),
     Aspect_(Aspect),
     NumSolutions_(NumSolutions_),
     CurrentSolution_(CurrentSolution),
     ClosedLoopStart_(ClosedLoopStart),
     FirstSolution_(FirstSolution),
     Difference_(Difference)
{
   ArcLenSet(dofs);
}

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,PerlInput &Input,
                                     const Vector &one,const Vector &two,int Echo)
   : Echo_(Echo),
     Mode_(Mode),
     CurrentSolution_(0),
     Difference_(two-one)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getUnsigned(Hash,"MaxIterations");
   Tolerance_ = Input.getDouble(Hash,"Tolerance");
   DSMax_ = Input.getDouble(Hash,"DSMax");
   CurrentDS_ = Input.getDouble(Hash,"DSStart");
   DSMin_ = Input.getDouble(Hash,"DSMin");
   AngleCutoff_ = Input.getDouble(Hash,"AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash,"AngleIncrease");
   Aspect_ = Input.getDouble(Hash,"Aspect");
   NumSolutions_ = Input.getUnsigned(Hash,"NumSolutions");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getUnsigned(Hash,"ClosedLoopStart");
   }
   else
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

ArcLengthSolution::ArcLengthSolution(LatticeMode *Mode,PerlInput &Input,int Echo)
   :  Echo_(Echo),
      Mode_(Mode),
      CurrentSolution_(0)
{
   ModeDOFS_=Mode_->ModeDOF().Dim();
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","ArcLengthSolution");
   MaxIter_ = Input.getUnsigned(Hash,"MaxIterations");
   Tolerance_ = Input.getDouble(Hash,"Tolerance");
   DSMax_ = Input.getDouble(Hash,"DSMax");
   CurrentDS_ = Input.getDouble(Hash,"DSStart");
   DSMin_ = Input.getDouble(Hash,"DSMin");
   AngleCutoff_ = Input.getDouble(Hash,"AngleCutoff");
   AngleIncrease_ = Input.getDouble(Hash,"AngleIncrease");
   Aspect_ = Input.getDouble(Hash,"Aspect");
   NumSolutions_ = Input.getUnsigned(Hash,"NumSolutions");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getUnsigned(Hash,"ClosedLoopStart");
   }
   else
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   
   BisectTolerance_ = Tolerance_;
   
   const char *starttype = Input.getString("StartType","Type");

   if (!strcmp("Bifurcation",starttype))
   {
      // Set Difference and Lattice state
      double eps = Input.getDouble("StartType","Epsilon");
      
      Difference_.Resize(ArcLenDef().Dim());
      Input.getVector(Difference_,"StartType","Tangent");
      Difference_ *= eps;
      
      Vector stat(Difference_.Dim());
      Input.getVector(stat,"StartType","BifurcationPoint");
      // Set Lattice state to the bifurcation point
      ArcLenSet(stat);
      
      // Set FirstSolution
      FirstSolution_.Resize(stat.Dim());
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Input.getVector(FirstSolution_,"StartType","ClosedLoopFirstSolution");
      }
      else
      {
         FirstSolution_ = stat;
      }
   }
   else if (!strcmp("Continuation",starttype))
   {
      // Set Lattice state to Solution2
      Vector two(ArcLenDef().Dim());
      Input.getVector(two,"StartType","Solution2");
      ArcLenSet(two);
      
      // Get solution1
      Vector one(two.Dim());
      Input.getVector(one,"StartType","Solution1");
      // Set Difference_ to   two - one
      Difference_.Resize(two.Dim());
      Difference_ = two - one;
      
      // Set FirstSolution
      FirstSolution_.Resize(two.Dim());
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Input.getVector(FirstSolution_,"StartType","ClosedLoopFirstSolution");
      }
      else
      {
         FirstSolution_ = one;
      }
   }
   else if (!strcmp("ConsistenceCheck",starttype))
   {
      double ConsistencyEpsilon;
      int Width,
         Dim=ArcLenDef().Dim();
      Vector Solution1(Dim),
         Solution2(Dim);

      Input.getVector(Solution2,"StartType","Solution2");
      Input.getVector(Solution1,"StartType","Solution1");
      // Get Epsilon and Width
      ConsistencyEpsilon = Input.getDouble("StartType","ConsistenceEpsilon");
      Width = Input.getInt("Main","FieldWidth");

      cout << "FIX UP CONSISTENCYCHECK!!!!!!\n";
      //ConsistencyCheck(Solution1,Solution2,ConsistencyEpsilon,Width,out);
   }
   else
   {
      cerr << "Unknown StartType!" << "\n";
      exit(-1);
   }
}

Vector ArcLengthSolution::ArcLenForce(double DS,const Vector &Diff,
                                      double Aspect)
{
   static Vector force(ModeDOFS_);
   static Vector mdfc(ModeDOFS_-1);
   mdfc = Mode_->ModeForce();
   
   force[ModeDOFS_-1] = DS*DS - Diff[ModeDOFS_-1]*Diff[ModeDOFS_-1]/(Aspect*Aspect);
   for (unsigned i=0;i<ModeDOFS_-1;++i)
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
   
   for (unsigned i=0;i<ModeDOFS_-1;++i)
   {
      for (unsigned j=0;j<=ModeDOFS_-1;++j)
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
      for (int i=0;i<70;i++) cout << "="; cout << "\n";
      cout << "Consistency Check." << "\n";
      cout << "F(U + DeltaU) * Epsilon" << "\n";
   }
   for (int i=0;i<70;i++) out << "="; out << "\n";
   out << "Consistency Check." << "\n";
   out << "F(U + DeltaU) * Epsilon" << "\n";
   ArcLenUpdate(Difference_);
   Force = ConsistencyEpsilon*ArcLenForce(ConsistencyEpsilon,Difference_,1.0);
   if (Echo_)
   {
      cout << setw(Width) << Force << "\n";
      cout << "K(U + DeltaU) * Epsilon" << "\n";
   }
   out << setw(Width) << Force << "\n";
   out << "K(U + DeltaU) * Epsilon" << "\n";
   Stiff = ConsistencyEpsilon*ArcLenStiffness(Difference_,1.0);
   if (Echo_) cout << setw(Width) << Stiff << "\n";
   out << setw(Width) << Stiff << "\n";
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
      cout << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << "\n";
      cout << setw(Width) << PerturbedForce << "\n";
      cout << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << "\n";
      cout << setw(Width) << PerturbedStiff << "\n";
      cout << "Difference" << "\n";
      cout << setw(Width) << Force - PerturbedForce << "\n" << "\n";
      cout << setw(Width) << Stiff - PerturbedStiff << "\n";
   }
   
   out << "P(U + DeltaU) - P(U + DeltaU + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedForce << "\n";
   out << "Fi(U + DeltaU) - Fi(U + DeltaU + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedStiff << "\n";
   out << "Difference" << "\n";
   out << setw(Width) << Force - PerturbedForce << "\n" << "\n";
   out << setw(Width) << Stiff - PerturbedStiff << "\n";
   
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << "\n";
   }
   for (int i=0;i<70;i++) out << "="; out << "\n";
   
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
      if (Echo_) cout << "DS= " << CurrentDS_ << "\n";
      
      ArcLengthNewton(good);
      
      AngleTest = ArcLenAngle(OldDiff,Difference_,Aspect_);
      
      if (Echo_)
         cout << "AngleTest = " << AngleTest << "  Cutoff = " << AngleCutoff_ << "\n";
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
      if (Echo_) cout << "DS= " << CurrentDS_ << "\n";
   }
   
   if (!good)
   {
      cerr << "ArcLenghtSolution did not converge properly" << "\n";
   }
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((ArcLenDef() - FirstSolution_).Norm() < CurrentDS_))
   {
      // We are done -- set currentsolution to numsolutions
      cerr << "Closed Loop detected at Solution # " << CurrentSolution_
           << " --- Terminating!" << "\n";
      
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
      if (Echo_) cout << "\n";
#endif
   }
   while ((itr < MaxIter_) && ((RHS.Norm() > Tolerance_) || (Dx.Norm() > Tolerance_)));
   
   if (Echo_) cout << resetiosflags(ios::scientific) << "\n";
   
   if (itr >= MaxIter_)
   {
      cerr << "Convergence Not Reached!!! -- ArcLengthNewton" << "\n";
      good = 0;
   }
   else
   {
      good = 1;
   }
}

int ArcLengthSolution::OldFindCriticalPoint(int LHN,double LHEV,int RHN,double RHEV,
                                            Lattice *Lat,PerlInput &Input,
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
   
   if (Echo_) cout << "LHN = " << LHN << "\n" << "RHN = " << RHN << "\n";
   
   // Find bifurcation point and make sure we are on the back side edge
   while (((fabs(CurrentMinEV) > BisectTolerance_)
           || (CurrentTestValue == RighthandTestValue))
          && (loops < MaxIter_))
   {
      if (Echo_) cout << "OldMinEV = " << OldMinEV << "\n"
                      << "CurrentTestValue = " << CurrentTestValue << "\n"
                      << "CurrentMinEV = "  << CurrentMinEV << "\n";
      
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
      
      cout << "Current_DS = " << CurrentDS_ << "\n" << "\n";
      
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
   if (Echo_) cout << "\n";
   out << "\n";
   
   if (Echo_) cout << setw(Width) << Lat << "\n";
   out << setw(Width) << Lat << "\n";
      
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "=";
      out << "=";
   }
   if (Echo_) cout << "\n"; out << "\n";
   
   // Call Lattice function to do any Lattice Specific things
   //  abs(RighthandNulity - LefthandNulity) is the number of zero eigenvalues
   //  in a perfect situation. should check to see if this is found to be true.
   Lat->CriticalPointInfo(Mode_->DrDt(Difference_),abs(RighthandTestValue-LefthandTestValue),
                          BisectTolerance_,Width,out);
   
   if (Echo_) cout << "Success = 1" << "\n";
   out << "Success = 1" << "\n";
   
   // Reset Lattice and ArcLengthSolution
   ArcLenUpdate(OriginalDiff-Difference_);
   CurrentDS_ = OriginalDS;
   Difference_ = OriginalDiff;
   
   return 1;
}

int ArcLengthSolution::FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out)
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
          << "\n";
      if (Echo_)
         cout << "Note: TestFunctions found a discrepancy between the\n"
              << "Note: difference in number of negative Test Functions\n"
              << "Note: and the number of Test Functions that change sign\n"
              << "Note: from LeftHandSide to RightHandSide.  This is usually\n"
              << "Note: caused by having too large of a path-following stepsize."
              << "\n";
      TestValueDiff = -TestValueDiff;
   }
   
   cout << "TF_LHS = " << setw(15) << TF_LHS<< "\n";
   cout << "TF_RHS = " << setw(15) << TF_RHS << "\n";
   
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
   
   cout << "TestValueDiff = "<< TestValueDiff << "\n";
   if (Echo_)
   {
      for (int i = 0; i < TestValueDiff; i++)
      {
         cout << "i = " << i<< "\n";
         cout << "Index[i] =" << Index[i]<< "\n";
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
               //cout <<"i = " << i << "\n" <<  "CHECK POINT INDEX[CP] = " << Index[i] << "\n";
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
         if (Echo_) cout << "\n";
         in_string << "\n";
         
         if (Echo_) cout << setw(Width) << Lat << "\n";
         in_string << setw(Width) << Lat << "\n";         
         
         for (int i=0;i<70;i++)
         {
            if (Echo_) cout << "=";
            in_string << "=";
         }
         if (Echo_) cout << "\n";
         in_string << "\n";
         
         // Call Lattice function to do any Lattice Specific things
         Lat->CriticalPointInfo(Mode_->DrDt(Difference_),Multiplicity,
                                BisectTolerance_,Width,in_string);
         
         if (Echo_) cout << "Success = 1" << "\n";
         in_string << "Success = 1" << "\n";
         
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
   //cout << " a = " << a << "\n" << "fa = " << fa << "\n"
   //<<  "b = " << b << "\n" << "fb = " << fb << "\n"
   //<< "c = " << c << "\n" << "fc = " << fc << "\n";
   
   while (((fabs(fb) > BisectTolerance_) )&& (loops < MaxIter_))
   {
      if (Echo_)
      {
         cout << setprecision(30) << "CurrentMinTF = " << fb << "\n";
         cout << "CurrentDS_ =" << CurrentDS_ << setprecision(10) << "\n";
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
         //cout << "Minimal Root found! XM too small! " < "\n";
         CurrentDS_ = LastDS;
         Difference_ = LastDiff;
         ArcLenUpdate(Difference_);
         Lat->TestFunctions(CurrentTF,Lattice::CRITPT);
         fb=CurrentTF[track];
         if(Echo_)
         {
            cout <<setprecision(30)<<"CurrentMinTF = " << fb << "\n";
            cout << "CurrentDS_ =" << CurrentDS_ <<setprecision(10)<< "\n";
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
      
      //cout << "CurrentTF = " << "\n" << setw(15) << CurrentTF << "\n" << "\n";
      //cout << "CurrentTF[track] = " << "\n" << CurrentTF[track]<< "\n"<< "\n";
   }
}
