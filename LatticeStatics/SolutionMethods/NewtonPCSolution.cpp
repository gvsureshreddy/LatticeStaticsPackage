#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "ArcLengthSolution.h"

using namespace std;

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,const Vector &one,
                                   int CurrentSolution,int NumSolutions,double MaxDS,
                                   double CurrentDS,double cont_rate_nom,double delta_nom,
                                   double alpha_nom,double Converge,double MinDSRatio,
                                   const Vector &FirstSolution,
                                   int Direction,int ClosedLoopStart,int StopAtCPNum,int Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(CurrentSolution),
     NumSolutions_(NumSolutions),
     MaxDS_(MaxDS),
     CurrentDS_(CurrentDS),
     cont_rate_nom_(cont_rate_nom),
     delta_nom_(delta_nom),
     alpha_nom_(alpha_nom),
     Converge_(Converge),
     MinDSRatio_(MinDSRatio),
     ClosedLoopStart_(ClosedLoopStart),
     StopAtCPNum_(StopAtCPNum),
     Direction_(Direction),
     FirstSolution_(FirstSolution)
{
   int count = (Mode_->ModeDOF()).Dim();
   int count_minus_one = count -1;
   
   //QR Decomposition of Stiffness Matrix
   Matrix Q(count, count);
   Matrix R(count, count_minus_one);
   
   //Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
   QR(Mode->ModeStiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   for(int i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,PerlInput &Input,
                                   const Vector &one,int Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getPosInt(Hash,"ClosedLoopStart");
   }
   else
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }

   if (Input.ParameterOK(Hash,"StopAtCPNum"))
   {
      StopAtCPNum_ = Input.getInt(Hash,"StopAtCPNum");
   }
   else
   {
      // Set default value
      StopAtCPNum_ = -1;
   }

   if (Input.ParameterOK(Hash,"Direction"))
   {
      Direction_ = Input.getInt(Hash,"Direction");
      if ((Direction_ < -1) || (Direction_ > 1))
      {
         cerr << "Unknown value for " << Hash.Name << "{Direction}" << "\n";
         exit(-6);
      }
   }
   else
   {
      // Default to positive;
      Direction_ = 1;
   }

   CurrentDS_ = MaxDS_;
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   Mode_->SetModeDOF(one);
   
   Previous_Solution_.Resize(one.Dim());
   
   int count = (Mode_->ModeDOF()).Dim();
   int count_minus_one = count -1;
   
   //QR Decomposition of Stiffness Matrix
   Matrix Q(count, count);
   Matrix R(count, count_minus_one);
   
   //Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
   QR(Mode->ModeStiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   for(int i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,PerlInput &Input,int Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getPosInt(Hash,"ClosedLoopStart");
   }
   else
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }

   if (Input.ParameterOK(Hash,"StopAtCPNum"))
   {
      StopAtCPNum_ = Input.getInt(Hash,"StopAtCPNum");
   }
   else
   {
      // Set default value
      StopAtCPNum_ = -1;
   }

   if (Input.ParameterOK(Hash,"Direction"))
   {
      Direction_ = Input.getInt(Hash,"Direction");
      if ((Direction_ < -1) || (Direction_ > 1))
      {
         cerr << "Unknown value for " << Hash.Name << "{Direction}" << "\n";
         exit(-6);
      }
   }
   else
   {
      // Default to positive;
      Direction_ = 1;
   }

   CurrentDS_ = MaxDS_;

   const char *starttype = Input.getString("StartType","Type");
   if (!strcmp("Bifurcation",starttype))
   {
      // Bifurcation
      // Get solution1
      int count = (Mode_->ModeDOF()).Dim();
      int i;
      Vector one(count);
      Previous_Solution_.Resize(count);
      Tangent1_.Resize(count);
      Tangent2_.Resize(count);

      Input.getVector(Tangent1_,"StartType","Tangent");
      Input.getVector(one,"StartType","BifurcationPoint");
      // override direction with start file value
      if (Input.ParameterOK("StartType","Direction"))
      {
         Direction_ = Input.getInt("StartType","Direction");
      }
      
      FirstSolution_.Resize(one.Dim());
      FirstSolution_ = one;
      
      Mode_->SetModeDOF(one);
      
      for(i=0; i<count; i++)
      {
         Tangent1_[i] = Direction_ * Tangent1_[i];
         Tangent2_[i] = Tangent1_[i];
      }
   }
   else if (!strcmp("Continuation",starttype))
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

      Input.getVector(one,"StartType","Solution1");
      // override direction with start file value
      if (Input.ParameterOK("StartType","Direction"))
      {
         Direction_ = Input.getInt("StartType","Direction");
      }
      
      FirstSolution_.Resize(one.Dim());
      FirstSolution_ = one;
      Mode_->SetModeDOF(one);
      
      Matrix Q(count, count);
      Matrix R(count, count_minus_one);
      
      QR(Mode_->ModeStiffness(),Q,R,1);
      for(i=0;i<count;i++)
      {
         Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
      }
   }
   else if (!strcmp("ConsistenceCheck",starttype))
   {
      // ConsistencyCheck
      
      // do nothing for now
   }
   else
   {
      cerr << "Unknown StartType!" << "\n";
      exit(-1);
   }
}

int NewtonPCSolution::AllSolutionsFound()
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

int NewtonPCSolution::FindNextSolution()
{
   //Finds the next solution
   static int count = FirstSolution_.Dim();
   static int count_minus_one = count-1;
   static Vector v(count);
   static Vector w(count);
   static Vector Force(count_minus_one);
   static Vector Corrector(count);
   static Matrix Q(count,count);
   static Matrix R(count,count_minus_one);
   int good=0;
   int omega=1;
   int i, Converge_Test;
   double Kappa, Alpha, Delta, Magnitude1, Magnitude2, temp, f;
   
   Previous_Solution_ = Mode_->ModeDOF();
   
   temp=0.0;
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
   
   //Starts solver
   do
   {
      for (i=0;i< count;i++)
      {
         v[i] = Previous_Solution_[i] + CurrentDS_ * Tangent1_[i];
      }
      
      //Sets state to predicted point
      Mode_->SetModeDOF(v);
      Force = Mode_->ModeForce();
      
      QR(Mode_->ModeStiffness(), Q, R, 1);
      
      for(i=0;i<count;i++)
      {
         Tangent2_[i] = Direction_ * Q[i][count_minus_one]*omega;
      }
      
      MoorePenrose(Q,R, Force, Corrector);
      
      Magnitude1 = Corrector.Norm();
      
      //CORRECTOR LOOP STARTS HERE
      Converge_Test = 0;
      do
      {
         
         for (i=0;i<count;i++)
         {
            w[i] = v[i] - Corrector[i];
         }
         
         Mode_->SetModeDOF(w);
         
         Force = Mode_->ModeForce();
         
         MoorePenrose(Q,R, Force,Corrector);
         Magnitude2 = Corrector.Norm();
         
         temp = 0.0;
         for (i=0; i<count; i++)
         {
            temp = temp + (Tangent1_[i]*Tangent2_[i]);
         }
         
         //checks parameters for steplength adaptation
         Kappa = sqrt((Magnitude2/ Magnitude1)/cont_rate_nom_);
         Alpha = sqrt(acos(temp)/alpha_nom_);
         Delta = sqrt(Magnitude1/delta_nom_);
         //cout << "Kappa = " << Kappa << "\n";
         //cout << "Alpha = " << Alpha << "\n";
         //cout << "Delta = " << Delta << "\n";
         
         temp = max(Kappa, Alpha);
         f = max(temp, Delta);
         
         temp = min(f,2.0);
         f = max(temp, 0.5);
         
         if(f >= 2.0)
         {
            //cout << " STEPLENGTH TOO LARGE " << "\n" << "\n";
            CurrentDS_ = CurrentDS_/2.0;
            if(CurrentDS_/MaxDS_ < MinDSRatio_)
            {
               cout << "Minimum StepSize ratio violated. Exit Solver. "<< "\n" << "\n";
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
            
            // cout << "STEPLENGTH OKAY " << "\n" << "\n" << "\n";
            if (Force.Norm() <= Converge_ && Corrector.Norm() <= Converge_)
            {
               Converge_Test = 1;
               CurrentDS_ = CurrentDS_ * 2.0;
               if(CurrentDS_ > MaxDS_)
               {
                  CurrentDS_ = MaxDS_;
               }
            }
            else
            {
               //cout << "HAS NOT CONVERGED " << "\n" << "\n" << "\n";
               QR(Mode_->ModeStiffness(), Q, R, 1);
               MoorePenrose(Q,R, Force,Corrector);
            }
         }
      }
      while (Converge_Test != 1);
   }
   while (f >= 2.0);
   //cout << "HAS CONVERGED " << "\n" << "\n" << "\n";
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Mode_->ModeDOF() - FirstSolution_).Norm() < MaxDS_))
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

   // always have current solution point printed
   good = 1;

   return good;
}

int NewtonPCSolution::FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out)
{
   static int TotalNumCPs=0;
   int NumCPs;
   int sz=Previous_Solution_.Dim();
   Vector tmp_diff(sz),tmp_DOF(Mode_->ModeDOF());
   double tmp_ds=0.0;
   for (int i=0;i<sz;++i)
   {
      tmp_ds += (Previous_Solution_[i]-tmp_DOF[i])*(Previous_Solution_[i]-tmp_DOF[i]);
   }
   tmp_ds = sqrt(tmp_ds);
   
   //ArcLengthSolution S1(Mode_,Input,Previous_Solution_,Mode_->ModeDOF(),1);
   int MaxIter = 50;
   ArcLengthSolution S1(Mode_,Mode_->ModeDOF(),MaxIter,Converge_,0.1*Converge_,tmp_ds,
                        tmp_ds,tmp_ds,1.0,0.5,1.0,1,0,Previous_Solution_,
                        Mode_->ModeDOF()-Previous_Solution_,10,Echo_);
   NumCPs=S1.FindCriticalPoint(Lat,Input,Width,out);
   TotalNumCPs += NumCPs;

   // Check to see if we should stop
   if ((StopAtCPNum_ > -1) && (TotalNumCPs >= StopAtCPNum_))
      CurrentSolution_ = NumSolutions_;
   
   return NumCPs;
}

void NewtonPCSolution::MoorePenrose(const Matrix& Q,const Matrix& R,const Vector& Force,
                                    Vector& Corrector)
{
   double sum;
   int i,j;
   int k=0;
   static int Size = Force.Dim()+1;
   static int count_minus_one = Size-1;
   static Vector y(Size);
   
   for(i=0;i < Size; i++)
   {
      y[i] = 0.0;
   }
   
   for (i=0; i<count_minus_one; i++)
   {
      sum = 0;
      for (j=0; j<k; j++)
      {
         sum += R[j][i]*y[j];
      }
      y[i] = (Force[i] - sum)/R[i][i];
      k++;
   }
   
   for (i=0; i<Size; i++)
   {
      sum = 0.0;
      for(j=0; j<Size; j++)
      {
         sum = sum + Q[i][j]*y[j];
      }
      Corrector[i] = sum;
   }
}
