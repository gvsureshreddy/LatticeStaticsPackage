#include "NewtonUpdatePCSolution.h"
#include "Matrix.h"
#include "UtilityFunctions.h"
#include <cmath>
#include "ArcLengthSolution.h"

using namespace std;

#define CLOSEDDEFAULT 30

NewtonUpdatePCSolution::NewtonUpdatePCSolution(LatticeMode *Mode,char *datafile,
					       const char *prefix,const Vector &one,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   // get needed parameters
   if(!GetParameter(prefix,"PCNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
   if(!GetParameter(prefix,"PCStepLength",datafile,'l',&MaxDS_)) exit(-1);
   if(!GetParameter(prefix,"PCContraction",datafile,'l',&cont_rate_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCDistance",datafile,'l',&delta_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCAngle",datafile,'l',&alpha_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCConvergeCriteria",datafile,'l',&Converge_)) exit(-1);
   if(!GetParameter(prefix,"MinDSRatio",datafile,'l',&MinDSRatio_)) exit(-1);
   if(!GetParameter(prefix,"PCClosedLoopStart",datafile,'u',&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   
   CurrentDS_ = MaxDS_;
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   Mode_->SetModeDOF(one);
   
   int count = Mode_->ModeStiffness().Cols();
   int count_minus_one = Mode_->ModeStiffness().Rows();
   
   //QR Decomposition of Stiffness Matrix
   
   Matrix Q(count, count);
   Matrix R(count, count_minus_one);

   //Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
   QR(Mode->ModeStiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   
   for(int i=0;i< count;i++)
   {
      Tangent1_[i] = Q[i][count_minus_one];
      Tangent2_[i] = Q[i][count_minus_one];
   }
}

NewtonUpdatePCSolution::NewtonUpdatePCSolution(LatticeMode *Mode,char *datafile,
					       const char *prefix,char *startfile,
					       fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   // get needed parameters12
   if(!GetParameter(prefix,"PCNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
   if(!GetParameter(prefix,"PCStepLength",datafile,'l',&MaxDS_)) exit(-1);
   if(!GetParameter(prefix,"PCContraction",datafile,'l',&cont_rate_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCDistance",datafile,'l',&delta_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCAngle",datafile,'l',&alpha_nom_)) exit(-1);
   if(!GetParameter(prefix,"PCConvergeCriteria",datafile,'l',&Converge_)) exit(-1);
   if(!GetParameter(prefix,"MinDSRatio",datafile,'l',&MinDSRatio_)) exit(-1);
   if(!GetParameter(prefix,"PCClosedLoopStart",datafile,'u',&ClosedLoopStart_,0))
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }
   
   CurrentDS_ = MaxDS_;
   
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
	 // Bifurcation
	 // Get solution1
	 int count = (Mode_->ModeDOF()).Dim();
	 int i;
	 Vector one(count);
	 Previous_Solution_.Resize(count);
	 Tangent1_.Resize(count);
	 Tangent2_.Resize(count);
	 
	 if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
	 if(!GetVectorParameter(prefix,"Tangent",startfile,&Tangent1_)) exit(-1);
	 
	 FirstSolution_.Resize(one.Dim());
	 FirstSolution_ = one;
	 
	 Mode_->SetModeDOF(one);
	 
	 for(i =0; i< count; i++)
	 {
	    Tangent2_[i] = Tangent1_[i];
	 }
	 
         break;
      }
      case 1:
      {
	 // Continuation
	 
         // Get solution1
	 int count = (Mode_->ModeDOF()).Dim();
	 int count_minus_one = count -1;
	 int i;
	 
	 Vector one(count);
	 
	 Previous_Solution_.Resize(count);
	 Tangent1_.Resize(count);
	 Tangent2_.Resize(count);
	 
	 if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
	 
	 FirstSolution_.Resize(one.Dim());
	 FirstSolution_ = one;
	 Mode_->SetModeDOF(one);
	 
	 Matrix Q(count, count);
	 Matrix R(count, count_minus_one);
	 QR(Mode_->ModeStiffness(),Q,R,1);
	 
	 for( i=0;i< count;i++)
	 {
	    Tangent1_[i] = Q[i][count_minus_one];
	    Tangent2_[i] = Tangent1_[i];
	 }
	 
	 break;
      }
      case 2:
      {
	 // ConsistencyCheck
	 
	 // do nothing for now
         break;
      }
   }
}

int NewtonUpdatePCSolution::AllSolutionsFound()
{
   if (CurrentSolution_ < NumSolutions_)
   {
      return 0;
   }
   else
   {
      return 1;
   }
}

double NewtonUpdatePCSolution::FindNextSolution(int &good)
{
   //Finds the next solution
   //Stiffness: NxN+1
   //DOF: N+1
   //Force: N
   static int count = FirstSolution_.Dim();
   static int count_minus_one = count -1;
   static Vector v(count);
   static Vector w(count);
   //static Vector Force1(count_minus_one);
   static Vector Force2(count_minus_one);
   static Vector Corrector(count);
   static Vector difference(count);
   static Matrix Q(count, count);
   static Matrix R(count, count_minus_one);
   static Matrix Stiffness(count_minus_one, count);
   static Matrix Temporary(count_minus_one, count);
   int omega=1;
   int i, j, Converge_Test ;
   double Kappa, Alpha, Delta, Magnitude1, Magnitude2, temp, f;
   
   Previous_Solution_ = Mode_->ModeDOF();
   
   //Force1 = Mode_->ModeForce();
   
   temp=0;
   for(i=0;i<count;i++)
   {
      temp = temp + Tangent1_[i] * Tangent2_[i];
   }
   
   if (temp < 0)
   {
      omega = -omega;
      
      for(i=0;i<count;i++)
      {
	 Tangent1_[i] = Tangent2_[i] * omega;
      }
   }
   else
   {
      for(i=0;i<count;i++)
      {
	 Tangent1_[i] = Tangent2_[i];
      }
   }
   
   do
   {
      for (i=0;i< count;i++)
      {
	 v[i] = Previous_Solution_[i] + CurrentDS_ * Tangent1_[i];
      }
      
      //cout << "START OF SOLVER " << endl << endl;
      Mode_->SetModeDOF(v);
      Force2 = Mode_->ModeForce();
      
      /*
      //obtains outer product and performs predictor update
      //A = A + h^-1 * (H(V) - H(u))*t(A)^T
      for (i=0; i < count_minus_one; i++)
      {
      for (j=0; j< count; j++)
      {
      Temporary[i][j] = (Force2[i] - Force1[i]) * Tangent1_[j];
      }
      for (j=0; j< count; j++)
      {
      Stiffness[i][j] = Stiffness[i][j] + (1/CurrentDS_)*Temporary[i][j];
      }
      }
      */
      Stiffness = Mode_->ModeStiffness();
      
      QR(Stiffness, Q, R, 1);
      
      for(i=0;i< count;i++)
      {
	 Tangent2_[i] = Q[i][count_minus_one]*omega;
      }
      
      MoorePenrose(Q,R, Force2, Corrector);
      
      Magnitude1 = Corrector.Norm();
      
      //CORRECTOR LOOP STARTS HERE
      Converge_Test = 0;
      do
      {
	 for (i=0;i<count;i++)
	 {
	    w[i] = v[i] - Corrector[i];
	    difference[i] = w[i] - v[i];
	 }
	 
	 Mode_->SetModeDOF(w);
	 
	 Force2 = Mode_->ModeForce();
	 
	 MoorePenrose(Q,R, Force2,Corrector);
	 
	 Magnitude2 = Corrector.Norm();
	 
	 temp = 0;
	 for (i=0; i<count ;i++)
	 {
	    temp = temp + (Tangent1_[i] * Tangent2_[i]);
	 }
	 
	 //checks parameters for steplength adaptation
	 Kappa = sqrt((Magnitude2/ Magnitude1)/cont_rate_nom_);
	 //Alpha = sqrt(acos(Tangent1_ * Tangent2_)/alpha_nom_);
	 Alpha = sqrt(acos(temp)/alpha_nom_);
	 Delta = sqrt(Magnitude1/delta_nom_);
	 //cout << setprecision(15)<< endl;
	 //cout << "Kappa = " << Kappa << endl;
	 //cout << "Alpha = " << Alpha << endl;
	 //cout << "Delta = " << Delta << endl;
	 
	 temp = max(Kappa, Alpha);
	 f = max(temp, Delta);
	 
	 temp = min(f,2.0);
	 f = max(temp, 0.5);
	 
	 //cout << " f = " << f << endl;
	 //cout << setprecision(10) << endl;
	 
	 if(f >= 2.0)
	 {
	    //cout << "STEPLENGTH TOO LARGE " << endl << endl << endl;
	    //cout << "Previous_Solution_ = " << endl << setw(15) << Previous_Solution_
	    //<< endl << endl;
	    //cout << "v = " <<endl << setw(15) << v << endl << endl;
	    //cout << "w = " << endl << setw(15) << w << endl << endl;
	    //CurrentDS_ = CurrentDS_/f;
	    CurrentDS_ = CurrentDS_/2.0;
	    if(CurrentDS_/MaxDS_ < MinDSRatio_)
	    {
	       cout << "Minimum StepSize ratio violated. Exit Solver. "<< endl << endl;
	       exit(-53);
	    }
	    break;
	 }
	 else
	 {
	    for(i=0;i<count;i++)
	    {
	       v[i] = w[i];
	    }
	    
	    //cout << "STEPLENGTH OKAY " << endl << endl << endl;
	    
	    if (Force2.Norm() <= Converge_ && Corrector.Norm() <= Converge_)
	    {
	       Converge_Test = 1;
	       
	       //CurrentDS_ = CurrentDS_/f;
	       CurrentDS_ = CurrentDS_*2.0;
	       if(CurrentDS_ > MaxDS_)
	       {
		  CurrentDS_ = MaxDS_;
	       }
	    }
	    else
	    {
	       //cout << "HAS NOT CONVERGED " << endl << endl << endl;
	       
	       temp = 0;
	       for(i=0;i<count;i++)
	       {
		  temp = temp + (difference[i] * difference[i]);
	       }
	       
	       for (i=0; i < count_minus_one; i++)
	       {
		  for (j=0; j< count; j++)
		  {
		     Temporary[i][j] = Force2[i] * difference[j];
		  }
		  for (j=0; j< count; j++)
		  {
		     Stiffness[i][j] = Stiffness[i][j] + (1/temp)*Temporary[i][j];
		  }
	       }
	       
	       QR(Stiffness, Q, R, 1);
	       
	       MoorePenrose(Q,R, Force2,Corrector);
	    }
	 }
      }
      while (Converge_Test != 1);
   }
   while (f >= 2.0);
   
   //cout << "HAS CONVERGED " << endl << endl << endl;
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Mode_->ModeDOF() - FirstSolution_).Norm() < MaxDS_))
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
}

int NewtonUpdatePCSolution::BisectAlert(int LHN,double LHEV,int RHN,double RHEV,Lattice *Lat,
					char *datafile,const char *prefix,int Width,
					fstream &out)
{
   ArcLengthSolution S1(Mode_, datafile, "^", Previous_Solution_,Mode_->ModeDOF(), 1);
   //S1.SetCurrentDS((Mode_->ModeDOF()-Previous_Solution_).Norm());
   int sz=Previous_Solution_.Dim();
   Vector tmp_diff(sz),tmp_DOF(Mode_->ModeDOF());
   double tmp_ds=0.0;
   for (int i=0;i<sz-1;++i)
   {
      tmp_ds += (Previous_Solution_[i]-tmp_DOF[i])*(Previous_Solution_[i]-tmp_DOF[i]);
   }
   tmp_ds += (Previous_Solution_[sz-1]-tmp_DOF[sz-1])*(Previous_Solution_[sz-1]-tmp_DOF[sz-1])
      /(S1.GetAspect()*S1.GetAspect());
   S1.SetCurrentDS(sqrt(tmp_ds));
   S1.BisectAlert(LHN,LHEV,RHN,RHEV,Lat,datafile,"^",Width,out);
   
   return 1;
}

void NewtonUpdatePCSolution::MoorePenrose(const Matrix& Q,const Matrix& R,const Vector& Force,
					  Vector& Corrector)
{
   double sum;
   int i,j;
   int k=0;
   static int Size = Force.Dim()+1;
   static int count_minus_one = Size -1;
   static Vector y(Size);
   
   for ( i = 0; i < count_minus_one; i++)
   {
      sum = 0;
      for ( j = 0; j < k;j++)
      {
	 sum +=R[j][i] * y[j];
      }
      
      y[i] = (Force[i] - sum)/R[i][i];
      
      k++;
   }
   
   for (i = 0; i < Size; i++)
   {
      sum = 0;
      
      for(j = 0; j < Size; j++)
      {
	 sum = sum + Q[i][j] * y[j];
      }
      
      Corrector[i] = sum;
   }
}
