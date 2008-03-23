#include <cmath>
#include "NewtonQRUpdatePCSolution.h"
#include "Matrix.h"
#include "ArcLengthSolution.h"

using namespace std;

NewtonQRUpdatePCSolution::NewtonQRUpdatePCSolution(LatticeMode *Mode,
                                                   const Vector &one,unsigned CurrentSolution,
                                                   unsigned NumSolutions,double MaxDS,
                                                   double CurrentDS,double cont_rate_nom,
                                                   double delta_nom,double alpha_nom,
                                                   double Converge,double MinDSRatio,
                                                   const Vector &FirstSolution,int Direction,
                                                   unsigned ClosedLoopStart,int Echo)
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

NewtonQRUpdatePCSolution::NewtonQRUpdatePCSolution(LatticeMode *Mode,PerlInput &Input,
                                                   const Vector &one,int Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonQRUpdatePCSolution");
   NumSolutions_ = Input.getUnsigned(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"StepLength");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getUnsigned(Hash,"ClosedLoopStart");
   }
   else
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }

   if (Input.ParameterOK(Hash,"Direction"))
   {
      Direction_ = Input.getInt(Hash,"Direction");
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
      Tangent1_[i] = Tangent2_[i] = Direction_ * Q[i][count_minus_one];
   }
}

NewtonQRUpdatePCSolution::NewtonQRUpdatePCSolution(LatticeMode *Mode,PerlInput &Input,int Echo)
   : Mode_(Mode),
     Echo_(Echo),
     CurrentSolution_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonQRUpdatePCSolution");
   NumSolutions_ = Input.getUnsigned(Hash,"NumSolutions");
   MaxDS_ = Input.getDouble(Hash,"StepLength");
   cont_rate_nom_ = Input.getDouble(Hash,"Contraction");
   delta_nom_ = Input.getDouble(Hash,"Distance");
   alpha_nom_ = Input.getDouble(Hash,"Angle");
   Converge_ = Input.getDouble(Hash,"ConvergeCriteria");
   MinDSRatio_ = Input.getDouble(Hash,"MinDSRatio");
   if (Input.ParameterOK(Hash,"ClosedLoopStart"))
   {
      ClosedLoopStart_ = Input.getUnsigned(Hash,"ClosedLoopStart");
   }
   else
   {
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
   }

   if (Input.ParameterOK(Hash,"Direction"))
   {
      Direction_ = Input.getInt(Hash,"Direction");
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

      Input.getVector(one,"StartType","Solution1");
      Input.getVector(Tangent1_,"StartType","Tangent");
      // override direction with start file value
      if (Input.ParameterOK("StartType","Direction"))
      {
         Direction_ = Input.getUnsigned("StartType","Direction");
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
         Direction_ = Input.getUnsigned("StartType","Direction");
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
   if (!strcmp("ConsistenceCheck",starttype))
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

int NewtonQRUpdatePCSolution::AllSolutionsFound()
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

int NewtonQRUpdatePCSolution::FindNextSolution()
{
   //Finds the next solution
   //Stiffness: NxN+1
   //DOF: N+1
   //Force: N
   static int count = FirstSolution_.Dim();
   static int count_minus_one = count -1;
   static Vector v(count);
   static Vector w(count);
   static Vector Force(count_minus_one);
   static Vector Corrector(count);
   static Vector difference(count);
   static Matrix Q(count, count);
   static Matrix R(count, count_minus_one);
   int good=0;
   int omega=1;
   int i, Converge_Test ;
   double Kappa, Alpha, Delta, Magnitude1, Magnitude2, temp, f;
   
   Previous_Solution_ = Mode_->ModeDOF();
   
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
      //cout << "START OF SOLVER " << "\n" << "\n";
      Mode_->SetModeDOF(v);
      Force = Mode_->ModeForce();
      
      QR(Mode_->ModeStiffness(), Q, R, 1);
      
      for(i=0;i< count;i++)
      {
         Tangent2_[i] =Direction_* Q[i][count_minus_one]*omega;
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
            difference[i] = w[i] - v[i];
         }
         
         Mode_->SetModeDOF(w);
         
         Force = Mode_->ModeForce();
         
         MoorePenrose(Q,R, Force,Corrector);
         
         Magnitude2 = Corrector.Norm();
         
         temp = 0;
         for (i=0; i<count ;i++)
         {
            temp = temp + (Tangent1_[i] * Tangent2_[i]);
         }
         
         //checks parameters for steplength adaptation
         Kappa = sqrt((Magnitude2/ Magnitude1)/cont_rate_nom_);
         Alpha = sqrt(acos(temp)/alpha_nom_);
         Delta = sqrt(Magnitude1/delta_nom_);
         //cout << setprecision(15)<< "\n";
         //cout << "Kappa = " << Kappa << "\n";
         //cout << "Alpha = " << Alpha << "\n";
         //cout << "Delta = " << Delta << "\n";
         
         temp = max(Kappa, Alpha);
         f = max(temp, Delta);
         
         temp = min(f,2.0);
         f = max(temp, 0.5);
         
         if(f >= 2.0)
         {
            //cout << "STEPLENGTH TOO LARGE " << "\n" << "\n" << "\n";
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
            
            //cout << "STEPLENGTH OKAY " << "\n" << "\n" << "\n";
            if (Force.Norm() <= Converge_ && Corrector.Norm() <= Converge_)
            {
               Converge_Test = 1;
               
               CurrentDS_ = CurrentDS_*2.0;
               if(CurrentDS_ > MaxDS_)
               {
                  CurrentDS_ = MaxDS_;
               }
            }
            else
            {
               //cout << "HAS NOT CONVERGED " << "\n" << "\n" << "\n";
               QRUpdate(Force, difference, Q,R);
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

int NewtonQRUpdatePCSolution::FindCriticalPoint(Lattice *Lat,PerlInput &Input,
                                                int Width,fstream &out)
{
   int sz=Previous_Solution_.Dim();
   Vector tmp_diff(sz),tmp_DOF(Mode_->ModeDOF());
   double tmp_ds=0.0;
   for (int i=0;i<sz;++i)
   {
      tmp_ds += (Previous_Solution_[i]-tmp_DOF[i])*(Previous_Solution_[i]-tmp_DOF[i]);
   }
   
   //ArcLengthSolution S1(Mode_,Input,Previous_Solution_,Mode_->ModeDOF(),1);
   int MaxIter = 50;
   ArcLengthSolution S1(Mode_,Mode_->ModeDOF(),MaxIter,Converge_,Converge_,tmp_ds,
                        tmp_ds,tmp_ds,1.0,0.5,1.0,1,0,Previous_Solution_,
                        Mode_->ModeDOF()-Previous_Solution_,10,Echo_);
   S1.FindCriticalPoint(Lat,Input,Width,out);
   return 1;
}

void NewtonQRUpdatePCSolution::MoorePenrose(const Matrix& Q,const Matrix& R,const Vector& Force,
                                            Vector& Corrector)
{
   double sum;
   int i,j;
   int k=0;
   static int Size = Force.Dim()+1;
   static int count_minus_one = Size -1;
   static Vector y(Size);
   
   for(i=0;i < Size; i++)
   {
      y[i] = 0.0;
   }
   
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


void NewtonQRUpdatePCSolution::QRUpdate(const Vector& Force,  const Vector& difference,
                                        Matrix& QBar, Matrix& RBar)
{
   static int count(Force.Dim());
   static int count_plus = count + 1;
   static Vector u(count_plus);
   static Vector a(count);
   static Vector e(count_plus);
   double norm;
   double C,S,r,A1,A2;
   double sum;
   
   norm =0.0;
   for(int i =0;i<count_plus;i++)
   {
      norm = norm + difference[i]*difference[i];
   }
   norm = sqrt(norm);
   
   //set a = Force/norm;
   for(int i=0;i<count;i++)
   {
      a[i] = Force[i]/norm;
      e[i] = difference[i]/norm;
   }
   e[count] = difference[count]/norm;
   
   //set Qbar =Q , Rbar = R and u = Q^T *e
   for(int i =0;i< count_plus;i++)
   {
      sum = 0.0;
      for(int w=0;w<count_plus;w++)
      {
         sum = sum + QBar[w][i]*e[w];
      }
      u[i] = sum;
   }
   
   //Algorithm 16.3.3 in Intro To Numerical Continuation Methods-- Algower, Georg
   for (int i = count-1;i>=0;i--)
   {
      //Calculate Rotation Cosine and Sine
      C = u[i];
      S = u[i+1];
      
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
         A1 = C * u[i] + S * u[i+1];
         A2 = C * u[i+1] - S * u[i];
         u[i] = A1;
         u[i+1] = A2;
      }
   }
   
   for(int i=0;i<count;i++)
   {
      RBar[0][i] = RBar[0][i] + u[0]*a[i];
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
