#include <cmath>
#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "ArcLengthSolution.h"

using namespace std;

NewtonPCSolution::NewtonPCSolution(LatticeMode* const Mode,
                                   Vector const& one,int const& CurrentSolution,
                                   int const& UpdateType,int const& NumSolutions,
                                   double const& MaxDS,double const& CurrentDS,
                                   double const& cont_rate_nom,double const& delta_nom,
                                   double const& alpha_nom,double const& Converge,
                                   double const& MinDSRatio,Vector const& FirstSolution,
                                   int const& Direction,int const& ClosedLoopStart,
                                   int const& StopAtCPNum,int const& Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(CurrentSolution),
     UpdateType_(UpdateType),
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
   //initialize "static" variables
   v_static.Resize(count);
   w_static.Resize(count);
   Force_static.Resize(count_minus_one);
   Corrector_static.Resize(count);
   difference_static.Resize(count);
   Q_static.Resize(count,count);
   R_static.Resize(count,count_minus_one);
   Stiff_static.Resize(count_minus_one,count);
   
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

NewtonPCSolution::NewtonPCSolution(LatticeMode* const Mode,PerlInput const& Input,
                                   Vector const& one,int const& Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   if (Input.ParameterOK(Hash,"UpdateType"))
   {
      const char* const updatetype=Input.getString(Hash,"UpdateType");
      if (!strcmp("None",updatetype))
      {
         UpdateType_ = 2;
      }
      else if (!strcmp("Stiffness",updatetype))
      {
         UpdateType_ = 1;
      }
      else if (!strcmp("QR",updatetype))
      {
         UpdateType_ = 0;
      }
      else
      {
         cerr << "Unknown UpdateType: " << updatetype << "\nExiting!\n";
         exit(-21);
      }
   }
   else
   {
      // default to QR
      UpdateType_ = 0;
   }
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   CurrentDS_ = Input.getDouble(Hash,"CurrentDS");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
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
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   Mode_->SetModeDOF(one);

   Previous_Solution_.Resize(one.Dim());
   
   int count = (Mode_->ModeDOF()).Dim();
   int count_minus_one = count - 1;
   //initialize "static" variables
   v_static.Resize(count);
   w_static.Resize(count);
   Force_static.Resize(count_minus_one);
   Corrector_static.Resize(count);
   difference_static.Resize(count);
   Q_static.Resize(count,count);
   R_static.Resize(count,count_minus_one);
   Stiff_static.Resize(count_minus_one,count);
   
   //QR Decomposition of Stiffness Matrix
   Matrix Q(count, count);
   Matrix R(count, count_minus_one);
   
   //Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
   QR(Mode->ModeStiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   
   for(int i=0;i< count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(LatticeMode* const Mode,PerlInput const& Input,
                                   int const& Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   if (Input.ParameterOK(Hash,"UpdateType"))
   {
      const char* const updatetype=Input.getString(Hash,"UpdateType");
      if (!strcmp("None",updatetype))
      {
         UpdateType_ = 2;
      }
      else if (!strcmp("Stiffness",updatetype))
      {
         UpdateType_ = 1;
      }
      else if (!strcmp("QR",updatetype))
      {
         UpdateType_ = 0;
      }
      else
      {
         cerr << "Unknown UpdateType: " << updatetype << "\nExiting!\n";
         exit(-21);
      }
   }
   else
   {
      // default to QR
      UpdateType_ = 0;
   }
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   CurrentDS_ = Input.getDouble(Hash,"CurrentDS");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
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
   
   int count = (Mode_->ModeDOF()).Dim();
   int count_minus_one = count - 1;
   //initialize "static" variables
   v_static.Resize(count);
   w_static.Resize(count);
   Force_static.Resize(count_minus_one);
   Corrector_static.Resize(count);
   difference_static.Resize(count);
   Q_static.Resize(count,count);
   R_static.Resize(count,count_minus_one);
   Stiff_static.Resize(count_minus_one,count);
   
   const char *starttype = Input.getString("StartType","Type");
   if (!strcmp("Bifurcation",starttype))
   {
      // Bifurcation
      // Get solution1
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

int NewtonPCSolution::AllSolutionsFound() const
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
   //Stiffness: NxN+1
   //DOF: N+1
   //Force: N
   int count = FirstSolution_.Dim();
   int count_minus_one = count -1;
   int oldprecision = cout.precision();
   int good=0;
   int omega=1;
   int i, Converge_Test;
   int iterations=0;
   double Kappa, Alpha, Delta, Magnitude1, Magnitude2, temp, f;
   double forcenorm;
   
   Previous_Solution_ = Mode_->ModeDOF();
   
   temp=0.0;
   for(i=0;i<count;i++)
   {
      temp = temp + Tangent1_[i] * Tangent2_[i];
   }
   
   if (temp < 0.0)
   {
      omega = -omega;
   }
   for(i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i] * omega;
   }
   
   do
   {
      if(CurrentDS_/MaxDS_ < MinDSRatio_)
      {
         cout << "Minimum StepSize ratio violated. Exit Solver.\n";
         exit(-53);
      }
      
      for (i=0;i<count;i++)
      {
         v_static[i] = Previous_Solution_[i] + CurrentDS_ * Tangent1_[i];
      }
      cout << "Taking Predictor Step. CurrentDS = " << CurrentDS_ << "\n";
      Mode_->SetModeDOF(v_static);
      Force_static = Mode_->ModeForce();
      
      Stiff_static = Mode_->ModeStiffness();
      QR(Stiff_static, Q_static, R_static, 1);
      
      for(i=0;i<count;i++)
      {
         Tangent2_[i] = Direction_* Q_static[i][count_minus_one]*omega;
      }
      
      MoorePenrose(Q_static,R_static, Force_static, Corrector_static);
      
      Magnitude1 = Corrector_static.Norm();
      
      //CORRECTOR LOOP STARTS HERE
      Converge_Test = 0;
      iterations = 0;
      do
      {
         for (i=0;i<count;i++)
         {
            w_static[i] = v_static[i] - Corrector_static[i];
            difference_static[i] = w_static[i] - v_static[i];
         }
         
         Mode_->SetModeDOF(w_static);
         
         Force_static = Mode_->ModeForce();
         forcenorm = Force_static.Norm();
         
         MoorePenrose(Q_static,R_static, Force_static,Corrector_static);
         
         Magnitude2 = Corrector_static.Norm();
         
         temp = 0.0;
         for (i=0;i<count;i++)
         {
            temp = temp + (Tangent1_[i] * Tangent2_[i]);
         }
         
         //checks parameters for steplength adaptation
         Kappa = sqrt((Magnitude2/ Magnitude1)/cont_rate_nom_);
         Alpha = sqrt(acos(temp)/alpha_nom_);
         Delta = sqrt(Magnitude1/delta_nom_);
         cout << setprecision(15);
         cout << "\tForceNorm = " << forcenorm
              << " \tCorrectorNorm = " << Magnitude2
              << " \tKappa = " << Kappa
              << " \tAlpha = " << Alpha
              << " \tDelta = " << Delta << "\n";
         cout << setprecision(oldprecision);
         
         temp = max(Kappa, Alpha);
         f = max(temp, Delta);
         
         temp = min(f,2.0);
         f = max(temp, 0.5);
         
         if(f >= 2.0)
         {
            CurrentDS_ = CurrentDS_/2.0;
            break;
         }
         else
         {
            for(i=0;i<count;i++)
            {
               v_static[i] = w_static[i];
            }
            
            if ((forcenorm <= Converge_) && (Magnitude2 <= Converge_))
            {
               Converge_Test = 1;
               
               CurrentDS_ = CurrentDS_/f;
               if(CurrentDS_ > MaxDS_)
               {
                  CurrentDS_ = MaxDS_;
               }
            }
            else
            {
               GetQR(Force_static,difference_static,Q_static,R_static);
               MoorePenrose(Q_static,R_static, Force_static,Corrector_static);
            }
         }
         
         ++iterations;
      }
      while (Converge_Test != 1);
      cout << "Corrector Iterations: " << iterations << "\n";
   }
   while (f >= 2.0);
   
   cout << "Converged with ForceNorm = " << forcenorm
        << " and CorrectorNorm = " << Magnitude2 << "\n";
   
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

int NewtonPCSolution::FindCriticalPoint(Lattice* const Lat,PerlInput const& Input,
                                        int const& Width,fstream& out)
{
   int TotalNumCPs=0;
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
   ArcLengthSolution S1(Mode_,Mode_->ModeDOF(),MaxIter,Converge_,Converge_,tmp_ds,
                        tmp_ds,tmp_ds,1.0,0.5,1.0,1,0,Previous_Solution_,
                        Mode_->ModeDOF()-Previous_Solution_,10,Echo_);
   NumCPs=S1.FindCriticalPoint(Lat,Input,Width,out);
   TotalNumCPs += NumCPs;
   
   // Check to see if we should stop
   if ((StopAtCPNum_ > -1) && (TotalNumCPs >= StopAtCPNum_))
      CurrentSolution_ = NumSolutions_;
   
   return NumCPs;
}

void NewtonPCSolution::MoorePenrose(Matrix const& Q,Matrix const& R,Vector const& Force,
                                    Vector& Corrector) const
{
   double sum;
   int i,j;
   int k=0;
   int Size = Force.Dim()+1;
   int count_minus_one = Size -1;
   
   y_static.Resize(Size,0.0);
   
   for ( i = 0; i < count_minus_one; i++)
   {
      sum = 0;
      for ( j = 0; j < k;j++)
      {
         sum +=R[j][i] * y_static[j];
      }
      
      y_static[i] = (Force[i] - sum)/R[i][i];
      k++;
   }
   
   for (i = 0; i < Size; i++)
   {
      sum = 0;
      for(j = 0; j < Size; j++)
      {
         sum = sum + Q[i][j] * y_static[j];
      }
      
      Corrector[i] = sum;
   }
}

void NewtonPCSolution::GetQR(Vector const& Force,Vector const& diff,Matrix& Q,Matrix& R) const
{
   int count = Mode_->ModeDOF().Dim();
   int count_minus_one = count - 1;
   int i,j;
   double temp = 0.0;
   
   for(i=0;i<count;i++)   switch (UpdateType_)
   {
      case 0:
         QRUpdate(Force,diff,Q,R);
         break;
      case 1:
         for(i=0;i<count;i++)
         {
            temp = temp + (diff[i] * diff[i]);
         }
         
         for (i=0; i < count_minus_one; i++)
         {
            for (j=0; j< count; j++)
            {
               Stiff_static[i][j] = Stiff_static[i][j] + (1/temp)*(Force[i] * diff[j]);
            }
         }
         
         QR(Stiff_static, Q, R, 1);
         break;
      case 2:
         QR(Mode_->ModeStiffness(),Q,R,1);
         break;
      default:
         cerr << "Unknown Update Type in NewtonPCSolution\n";
         exit(-20);
         break;
   }
}

void NewtonPCSolution::QRUpdate(Vector const& Force,Vector const& difference,Matrix& QBar,
                                Matrix& RBar) const
{
   int count(Force.Dim());
   int count_plus = count + 1;
   double norm;
   double C,S,r,A1,A2;
   double sum;
   
   u_static.Resize(count_plus);
   a_static.Resize(count);
   e_static.Resize(count_plus);
   
   norm =0.0;
   for(int i =0;i<count_plus;i++)
   {
      norm = norm + difference[i]*difference[i];
   }
   norm = sqrt(norm);
   
   //set a = Force/norm;
   for(int i=0;i<count;i++)
   {
      a_static[i] = Force[i]/norm;
      e_static[i] = difference[i]/norm;
   }
   e_static[count] = difference[count]/norm;
   
   //set Qbar =Q , Rbar = R and u = Q^T *e
   for(int i=0;i<count_plus;i++)
   {
      sum = 0.0;
      for(int w=0;w<count_plus;w++)
      {
         sum = sum + QBar[w][i]*e_static[w];
      }
      u_static[i] = sum;
   }
   
   //Algorithm 16.3.3 in Intro To Numerical Continuation Methods-- Algower, Georg
   for (int i = count-1;i>=0;i--)
   {
      //Calculate Rotation Cosine and Sine
      C = u_static[i];
      S = u_static[i+1];
      
      if(S!=0.0)
      {
         r = sqrt(C*C+S*S);
         C = C/r;
         S = S/r;
         
         //Givens Rotation on Rows of RBar
         for(int k=0;k<count;k++)
         {
            A1 = C * RBar[i][k] + S * RBar[i+1][k];
            A2 = C * RBar[i+1][k] - S * RBar[i][k];
            RBar[i][k] = A1;
            RBar[i+1][k] = A2;
            A1 = C * QBar[i][k] + S * QBar[i+1][k];
            A2 = C * QBar[i+1][k] - S * QBar[i][k];
            QBar[i][k] = A1;
            QBar[i+1][k] = A2;
            
         }
         A1 = C * QBar[i][count] + S * QBar[i+1][count];
         A2 = C * QBar[i+1][count] - S * QBar[i][count];
         QBar[i][count] = A1;
         QBar[i+1][count] = A2;
         
         //Givens Rotation on u
         A1 = C * u_static[i] + S * u_static[i+1];
         A2 = C * u_static[i+1] - S * u_static[i];
         u_static[i] = A1;
         u_static[i+1] = A2;
      }
   }
   
   for(int i=0;i<count;i++)
   {
      RBar[0][i] = RBar[0][i] + u_static[0]*a_static[i];
   }
   
   for(int i = 0;i<count;i++)
   {
      C = RBar[i][i];
      S = RBar[i+1][i];
      
      if(S!=0.0)
      {
         r = sqrt(C*C+S*S);
         C = C/r;
         S = S/r;
         
         //Givens Rotation on Rows of RBar
         for(int k=0;k<count;k++)
         {
            A1 = C * RBar[i][k] + S * RBar[i+1][k];
            A2 = C * RBar[i+1][k] - S * RBar[i][k];
            RBar[i][k] = A1;
            RBar[i+1][k] = A2;
            A1 = C * QBar[i][k] + S * QBar[i+1][k];
            A2 = C * QBar[i+1][k] - S * QBar[i][k];
            QBar[i][k] = A1;
            QBar[i+1][k] = A2;
         }
         A1 = C * QBar[i][count] + S * QBar[i+1][count];
         A2 = C * QBar[i+1][count] - S * QBar[i][count];
         QBar[i][count] = A1;
         QBar[i+1][count] = A2;
      }
   }
}
