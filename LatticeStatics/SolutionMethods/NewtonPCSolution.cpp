#include <cmath>
#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "ArcLengthSolution.h"

using namespace std;

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,
                                   Vector const& one,int const& CurrentSolution,
                                   UpdateType const& Type,int const& NumSolutions,
                                   double const& MaxDS,double const& CurrentDS,
                                   double const& MinDS,double const& cont_rate_max,
                                   double const& delta_max,double const& alpha_max,
                                   double const& Converge,Vector const& FirstSolution,
                                   int const& Direction,double const& accel_max,
                                   int const& ClosedLoopStart,int const& StopAtCPCrossingNum,
                                   int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(CurrentSolution),
     UpdateType_(Type),
     NumSolutions_(NumSolutions),
     MaxDS_(MaxDS),
     CurrentDS_(CurrentDS),
     MinDS_(MinDS),
     cont_rate_max_(cont_rate_max),
     delta_max_(delta_max),
     alpha_max_(alpha_max),
     Converge_(Converge),
     ClosedLoopStart_(ClosedLoopStart),
     StopAtCPCrossingNum_(StopAtCPCrossingNum),
     Direction_(Direction),
     Omega_(1.0),
     accel_max_(accel_max),
     FirstSolution_(FirstSolution),
     PreviousSolution_(FirstSolution_.Dim())
{
   if (cos(alpha_max_) <= 0.0)
   {
      cerr << "error: NewtonPCSolution::Angle too large!\n";
      exit(-22);
   }

   Restrict_->SetDOF(FirstSolution_);
   int count = (Restrict_->DOF()).Dim();
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
   QR(Restrict_->Stiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   for(int i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,PerlInput const& Input,
                                   Vector const& one,int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(0),
     Omega_(1.0),
     PreviousSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   if (Input.ParameterOK(Hash,"UpdateType"))
   {
      char const* const updatetype=Input.getString(Hash,"UpdateType");
      if (!strcmp("NoUpdate",updatetype))
      {
         UpdateType_ = NoUpdate;
      }
      else if (!strcmp("QRUpdate",updatetype))
      {
         UpdateType_ = QRUpdate;
      }
      else if (!strcmp("StiffnessUpdate",updatetype))
      {
         UpdateType_ = StiffnessUpdate;
      }
      else if (!strcmp("Exact",updatetype))
      {
         UpdateType_ = Exact;
      }
      else
      {
         cerr << "Unknown UpdateType: " << updatetype << "\nExiting!\n";
         exit(-21);
      }
   }
   else
   {
      Input.useString("QRUpdate",Hash,"UpdateType"); // Default Value
      UpdateType_ = QRUpdate;
   }
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   CurrentDS_ = Input.getDouble(Hash,"CurrentDS");
   MinDS_ = Input.getDouble(Hash,"MinDS");
   cont_rate_max_ = Input.getDouble(Hash,"Contraction");
   delta_max_ = Input.getDouble(Hash,"Distance");
   alpha_max_ = Input.getDouble(Hash,"Angle");
   if (cos(alpha_max_) <= 0.0)
   {
      cerr << "error: NewtonPCSolution::Angle too large!\n";
      exit(-22);
   }
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT,Hash,"ClosedLoopStart"); // Default Value
   }
   
   if (Input.ParameterOK(Hash,"StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash,"StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1,Hash,"StopAtCPCrossingNum"); // Default Value
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
      Direction_ = Input.useInt(1,Hash,"Direction");
   }

   if (Input.ParameterOK(Hash,"Acceleration"))
   {
      accel_max_ = Input.getDouble(Hash,"Acceleration");
      if (accel_max_ <= 0.0)
      {
         cerr << "Invalid value for " << Hash.Name << "{Acceleration}" << "\n";
         exit(-7);
      }
   }
   else
   {
      accel_max_ = Input.useDouble(2.0,Hash,"Acceleration"); // Default Value
   }
   Input.EndofInputSection();
   
   FirstSolution_.Resize(one.Dim());
   FirstSolution_ = one;
   Restrict_->SetDOF(one);

   PreviousSolution_.Resize(one.Dim());
   
   int count = (Restrict_->DOF()).Dim();
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
   QR(Restrict_->Stiffness(),Q,R,1);
   
   Tangent1_.Resize(count);
   Tangent2_.Resize(count);
   
   for(int i=0;i< count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,PerlInput const& Input,
                                   int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(0),
     Omega_(1.0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
   if (Input.ParameterOK(Hash,"UpdateType"))
   {
      const char* const updatetype=Input.getString(Hash,"UpdateType");
      if (!strcmp("NoUpdate",updatetype))
      {
         UpdateType_ = NoUpdate;
      }
      else if (!strcmp("QRUpdate",updatetype))
      {
         UpdateType_ = QRUpdate;
      }
      else if (!strcmp("StiffnessUpdate",updatetype))
      {
         UpdateType_ = StiffnessUpdate;
      }
      else if (!strcmp("Exact",updatetype))
      {
         UpdateType_ = Exact;
      }
      else
      {
         cerr << "Unknown UpdateType: " << updatetype << "\nExiting!\n";
         exit(-21);
      }
   }
   else
   {
      Input.useString("QRUpdate",Hash,"UpdateType"); // Default Value
      UpdateType_ = QRUpdate;
   }
   
   NumSolutions_ = Input.getPosInt(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"MaxDS");
   CurrentDS_ = Input.getDouble(Hash,"CurrentDS");
   MinDS_ = Input.getDouble(Hash,"MinDS");
   cont_rate_max_ = Input.getDouble(Hash,"Contraction");
   delta_max_ = Input.getDouble(Hash,"Distance");
   alpha_max_ = Input.getDouble(Hash,"Angle");
   if (cos(alpha_max_) <= 0.0)
   {
      cerr << "error: NewtonPCSolution::Angle too large!\n";
      exit(-22);
   }   
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getInt(Hash,"ClosedLoopStart");
   }
   else
   {
      ClosedLoopStart_ = Input.useInt(CLOSEDDEFAULT,Hash,"ClosedLoopStart"); // Default Value
   }
   
   if (Input.ParameterOK(Hash,"StopAtCPCrossingNum"))
   {
      StopAtCPCrossingNum_ = Input.getInt(Hash,"StopAtCPCrossingNum");
   }
   else
   {
      StopAtCPCrossingNum_ = Input.useInt(-1,Hash,"StopAtCPCrossingNum"); // Default Value
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
      Direction_ = Input.useInt(1,Hash,"Direction"); // Default Value
   }

   if (Input.ParameterOK(Hash,"Acceleration"))
   {
      accel_max_ = Input.getDouble(Hash,"Acceleration");
      if (accel_max_ <= 0.0)
      {
         cerr << "Invalid value for " << Hash.Name << "{Acceleration}" << "\n";
         exit(-7);
      }
   }
   else
   {
      accel_max_ = Input.useDouble(2.0,Hash,"Acceleration"); // Default Value
   }
   Input.EndofInputSection();
   
   int count = (Restrict_->DOF()).Dim();
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
      PreviousSolution_.Resize(count);
      Tangent1_.Resize(count);
      Tangent2_.Resize(count);

      Vector tan1tmp(Input.getArrayLength("StartType","Tangent"));
      Input.getVector(tan1tmp,"StartType","Tangent");
      Tangent1_ = Restrict_->TransformVector(tan1tmp);
      Tangent1_ = Tangent1_/Tangent1_.Norm();
      Vector biftmp(Input.getArrayLength("StartType","BifurcationPoint"));
      Input.getVector(biftmp,"StartType","BifurcationPoint");
      one = Restrict_->RestrictDOF(biftmp);
      
      FirstSolution_.Resize(one.Dim());
      FirstSolution_ = one;
      
      Restrict_->SetDOF(one);
      
      for(i=0; i<count; i++)
      {
         Tangent1_[i] = Tangent1_[i];
         Tangent2_[i] = Tangent1_[i];
      }
   }
   else if (!strcmp("Continuation",starttype))
   {
      // Continuation
      
      // Get solution1
      int i;
      Vector one(count);
      
      PreviousSolution_.Resize(count);
      Tangent1_.Resize(count);
      Tangent2_.Resize(count);

      Vector onetmp(Input.getArrayLength("StartType","Solution"));
      Input.getVector(onetmp,"StartType","Solution");
      one = Restrict_->RestrictDOF(onetmp);
      
      FirstSolution_.Resize(one.Dim());
      FirstSolution_ = one;
      Restrict_->SetDOF(one);
      
      Matrix Q(count, count);
      Matrix R(count, count_minus_one);
      
      QR(Restrict_->Stiffness(),Q,R,1);
      for(i=0;i<count;i++)
      {
         Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
      }
   }
   else if (!strcmp("ConsistencyCheck",starttype))
   {
      double ConsistencyEpsilon;
      int Width;
      Vector Solution(count);

      Vector onetmp(Input.getArrayLength("StartType","Solution"));
      Input.getVector(onetmp,"StartType","Solution");
      Solution = Restrict_->RestrictDOF(onetmp);
      // Get Epsilon and Width
      ConsistencyEpsilon = Input.getDouble("StartType","Epsilon");
      Width = Input.getPosInt("Main","FieldWidth");

      fstream::fmtflags oldflags=cout.flags();
      cout << scientific;
      Restrict_->ConsistencyCheck(Solution,ConsistencyEpsilon,Width,cout);
      cout.flags(oldflags);
      // We're done
      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      cerr << "Unknown StartType!" << "\n";
      exit(-1);
   }
   Input.EndofInputSection();
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
   int good=0;
   int i, Converge_Test;
   int predictions=1;
   int corrections=0;
   double Kappa=0.0;
   double Magnitude1=0.0;
   double Magnitude2=0.0;
   double Dot=0.0;
   double f;
   double forcenorm;
   double const eta = 0.1;
   
   PreviousSolution_ = Restrict_->DOF();
   
   for(i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i];
   }
   
   do
   {
      corrections = 0;
      
      if(CurrentDS_ < MinDS_)
      {
         cout << "Minimum StepSize (MinDS) violated. Exit Solver.\n";
         exit(-53);
      }

      // setup acceleration factor to as small as possible
      f = 1.0/accel_max_;
      
      for (i=0;i<count;i++)
      {
         v_static[i] = PreviousSolution_[i] + CurrentDS_*Omega_*Tangent1_[i];
      }
      cout << "Taking Predictor Step. CurrentDS = " << CurrentDS_;

      Restrict_->SetDOF(v_static);
      // get stiffness first for efficiency
      Stiff_static = Restrict_->Stiffness();
      Force_static = Restrict_->Force();
      forcenorm = Force_static.Norm();
      QR(Stiff_static, Q_static, R_static, 1);
      
      for(i=0;i<count;i++)
      {
         Tangent2_[i] = Direction_*Q_static[i][count_minus_one];
      }
      
      Dot=0.0;
      for (i=0;i<count;++i)
      {
         Dot += Tangent1_[i]*Tangent2_[i];
      }
      
      cout << " \tDot = " << Dot << " \tOmega = " << Omega_ << " \talpha = "
           << acos(fabs(Dot))*180.0/3.1415926 << "\n";
      if (acos(fabs(Dot)) > alpha_max_)
      {
         f = accel_max_;
         CurrentDS_ /= f;
         cout << "Prediction " << predictions << " Corrector Iterations: "
              << corrections << "\n";
         ++predictions;
         continue;
      }
      f = max(f,(acos(fabs(Dot))/alpha_max_)*accel_max_);

      MoorePenrose(Q_static,R_static, Force_static, Corrector_static);
      
      Magnitude1 = Corrector_static.Norm();
      Magnitude2 = Magnitude1;

      corrections++;
      cout << " \tForceNorm = " << forcenorm << " \tDeltaNorm = " << Magnitude1;
      if (Magnitude1 > delta_max_)
      {
         f = accel_max_;
         CurrentDS_ /= f;
         cout << "\nPrediction " << predictions << " Corrector Iterations: "
              << corrections << "\n";
         ++predictions;
         continue;
      }
      f = max(f,sqrt(Magnitude1/delta_max_)*accel_max_);
      cout << " \tAcceleration = " << 1.0/f << "\n";
      
      //CORRECTOR LOOP STARTS HERE
      Converge_Test = 0;
      do
      {
         for (i=0;i<count;i++)
         {
            w_static[i] = v_static[i] - Corrector_static[i];
            difference_static[i] = -Corrector_static[i];
         }
         
         Restrict_->SetDOF(w_static);
         
         Force_static = Restrict_->Force();
         forcenorm = Force_static.Norm();
         cout << "\tForceNorm = " << forcenorm;
         
         MoorePenrose(Q_static,R_static, Force_static,Corrector_static);

         Magnitude2 = Corrector_static.Norm();
         cout << " \tDeltaNorm = " << Magnitude2;
         if (Magnitude2 > delta_max_)
         {
            f = accel_max_;
            CurrentDS_ /= f;
            cout << "\n";
            break;
         }
         f = max(f,sqrt(Magnitude2/delta_max_)*accel_max_);

         Kappa = Magnitude2/(Magnitude1+Converge_*eta);
         cout << " \tContraction = " << Kappa;
         if (Kappa > cont_rate_max_)
         {
            f = accel_max_;
            CurrentDS_ /= f;
            cout << "\n";
            break;
         }
         f = max(f,sqrt(Kappa/cont_rate_max_)*accel_max_);
         cout << " \tAcceleration = " << 1.0/f << "\n";
         Magnitude1=Magnitude2;
         
         // accpet iteration
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
         
         ++corrections;
      }
      while (Converge_Test != 1);
      cout << "Prediction " << predictions << " Corrector Iterations: " << corrections << "\n";
      ++predictions;
   }
   while (f >= accel_max_);
   
   cout << "Converged with ForceNorm = " << forcenorm
        << " and CorrectorNorm = " << Magnitude2 << "\n";

   if (Dot <= 0.0)
   {
      Omega_ = -Omega_;
   }
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Restrict_->DOF() - FirstSolution_).Norm() < CurrentDS_))
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

void NewtonPCSolution::FindCriticalPoint(Lattice* const Lat,int& TotalNumCPCrossings,
                                         PerlInput const& Input,int const& Width,fstream& out)
{
   int sz=PreviousSolution_.Dim();
   Vector tmp_diff(sz),tmp_DOF(Restrict_->DOF());
   double tmp_ds=0.0;
   for (int i=0;i<sz;++i)
   {
      tmp_ds += (PreviousSolution_[i]-tmp_DOF[i])*(PreviousSolution_[i]-tmp_DOF[i]);
   }
   tmp_ds = sqrt(tmp_ds);
   
   //ArcLengthSolution S1(Restrict_,Input,PreviousSolution_,Restrict_->DOF(),1);
   int MaxIter = 50;
   ArcLengthSolution S1(Restrict_,Restrict_->DOF(),MaxIter,Converge_,Converge_,tmp_ds,
                        tmp_ds,tmp_ds,1.0,0.5,1.0,1,0,PreviousSolution_,
                        Restrict_->DOF()-PreviousSolution_,10,-1,Echo_);
   S1.FindCriticalPoint(Lat,TotalNumCPCrossings,Input,Width,out);
      
   // Check to see if we should stop
   if ((StopAtCPCrossingNum_ > -1) && (TotalNumCPCrossings >= StopAtCPCrossingNum_))
      CurrentSolution_ = NumSolutions_;
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
   int count = Restrict_->DOF().Dim();
   int count_minus_one = count - 1;
   int i,j;
   double temp = 0.0;

   switch (UpdateType_)
   {
      case NoUpdate:
         break;
      case QRUpdate:
         UpdateQR(Force,diff,Q,R);
         break;
      case StiffnessUpdate:
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
      case Exact:
         QR(Restrict_->Stiffness(),Q,R,1);
         break;
      default:
         cerr << "Unknown Update Type in NewtonPCSolution\n";
         exit(-20);
         break;
   }
}

void NewtonPCSolution::UpdateQR(Vector const& Force,Vector const& difference,Matrix& QBar,
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
