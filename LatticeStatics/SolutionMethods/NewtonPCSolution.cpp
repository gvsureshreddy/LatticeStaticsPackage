#include <cmath>
#include <fstream>
#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "ArcLengthSolution.h"

using namespace std;

NewtonPCSolution::~NewtonPCSolution()
{
   cout.width(0);
   cout << "NewtonPCSolution Function Calls:\n"
        << "\tGetQR - " << counter_[0] << "\n"
        << "\tMoorePenrose - " << counter_[1] << "\n"
        << "\tUpdateQR - " << counter_[2] << "\n"
        << "\tAllSolutionsFound - " << counter_[3] << "\n"
        << "\tFindNextSolution - " << counter_[4] << "\n"
        << "\tFindCriticalPoint - " << counter_[5] << "\n";
}

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,
                                   Vector const& one,int const& CurrentSolution,
                                   UpdateType const& Type,int const& NumSolutions,
                                   double const& MaxDS,double const& CurrentDS,
                                   double const& MinDS,double const& cont_rate_max,
                                   double const& delta_max,double const& alpha_max,
                                   double const& Converge,ConvergeType CnvrgTyp,
                                   Vector const& FirstSolution,int const& Direction,
                                   double const& accel_max,int const& BifStartFlag,
                                   Vector const& BifTangent,int const& ClosedLoopStart,
                                   int const& ClosedLoopUseAsFirst,
                                   int const& StopAtCPCrossingNum,int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(CurrentSolution),
     UpdateType_(Type),
     NumSolutions_(NumSolutions),
     MaxDS_(MaxDS),
     PreviousDS_(CurrentDS),
     CurrentDS_(CurrentDS),
     MinDS_(MinDS),
     cont_rate_max_(cont_rate_max),
     delta_max_(delta_max),
     alpha_max_(alpha_max),
     Converge_(Converge),
     ConvergeType_(CnvrgTyp),
     BifStartFlag_(BifStartFlag),
     BifTangent_(BifTangent),
     ClosedLoopStart_(ClosedLoopStart),
     ClosedLoopUseAsFirst_(ClosedLoopUseAsFirst),
     StopAtCPCrossingNum_(StopAtCPCrossingNum),
     Direction_(Direction),
     Omega_(1.0),
     accel_max_(accel_max),
     FirstSolution_(FirstSolution),
     PreviousSolution_(FirstSolution_.Dim())
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
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
   double tansign = 1.0;
   for (int i=0;i<count_minus_one;++i)
   {
      tansign *= R[i][i]/fabs(R[i][i]);
   }
   for (int i=0;i<count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * tansign * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,PerlInput const& Input,
                                   Vector const& one,int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(0),
     BifStartFlag_(0),
     BifTangent_(),
     Omega_(1.0),
     PreviousSolution_(0)
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
   
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
   PreviousDS_ = CurrentDS_;
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
   if (Input.ParameterOK(Hash,"ConvergeType"))
   {
      char const* const cnvrgtyp=Input.getString(Hash,"ConvergeType");
      if (!strcmp("Both",cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force",cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement",cnvrgtyp))
      {
         ConvergeType_ = Displacement;
      }
      else
      {
         cerr << "Unknown ConvergeType: " << cnvrgtyp << "\nExiting!\n";
         exit(-22);
      }
   }
   else
   {
      Input.useString("Both",Hash,"ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
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
   
   FirstSolution_.Resize(one.Dim());
   if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
   {
      Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
      Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
      FirstSolution_ = Restrict_->RestrictDOF(clfs);
   }
   else
   {
      if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
      {
         ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
         if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
         {
            cerr << "Error: NewtonPCSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                 << endl;
            exit(-33);
         }
      }
      else
      {
         // Default Value
         ClosedLoopUseAsFirst_ = Input.useInt(0,"StartType","ClosedLoopUseAsFirst");
         FirstSolution_ = one;
      }
   }
   Input.EndofInputSection();

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
   double tansign = 1.0;
   for (int i=0;i<count_minus_one;++i)
   {
      tansign *= R[i][i]/fabs(R[i][i]);
   }   
   for(int i=0;i< count;i++)
   {
      Tangent1_[i] = Tangent2_[i] = Direction_ * tansign * Q[i][count_minus_one];
   }
}

NewtonPCSolution::NewtonPCSolution(Restriction* const Restrict,PerlInput const& Input,
                                   int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     CurrentSolution_(0),
     BifStartFlag_(0),
     BifTangent_(),
     Omega_(1.0)
{
   for (int i=0;i<nocounters_;++i) counter_[i]=0;
   
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
   PreviousDS_ = CurrentDS_;
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
   if (Input.ParameterOK(Hash,"ConvergeType"))
   {
      char const* const cnvrgtyp=Input.getString(Hash,"ConvergeType");
      if (!strcmp("Both",cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force",cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement",cnvrgtyp))
      {
         ConvergeType_ = Displacement;
      }
      else
      {
         cerr << "Unknown ConvergeType: " << cnvrgtyp << "\nExiting!\n";
         exit(-22);
      }
   }
   else
   {
      Input.useString("Both",Hash,"ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }
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
      BifStartFlag_ = 1;
      BifTangent_.Resize(count);
      // Get solution1
      int i;
      Vector one(count);
      PreviousSolution_.Resize(count);
      Tangent1_.Resize(count);
      Tangent2_.Resize(count);
      
      Vector tan1tmp(Input.getArrayLength("StartType","Tangent"));
      Input.getVector(tan1tmp,"StartType","Tangent");
      Tangent1_ = Restrict_->TransformVector(tan1tmp);
      Tangent1_ = Tangent1_/(double(Direction_)*Tangent1_.Norm());
      // set up projection vector BifTangent_
      BifTangent_ = Tangent1_;
      BifTangent_[count-1] = 0.0;
      Vector biftmp(Input.getArrayLength("StartType","BifurcationPoint"));
      Input.getVector(biftmp,"StartType","BifurcationPoint");
      one = Restrict_->RestrictDOF(biftmp);
      
      FirstSolution_.Resize(one.Dim());
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
         Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
            if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
            {
               cerr << "Error: NewtonPCSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                    << endl;
               exit(-33);
            }
         }
         else
         {
            // Default Value
            ClosedLoopUseAsFirst_ = Input.useInt(0,"StartType","ClosedLoopUseAsFirst");
            FirstSolution_ = one;
         }
      }
      
      cout << "Projection on BifTangent of BifurcationPoint = " << FirstSolution_*BifTangent_
           << "\n";
      
      Restrict_->SetDOF(one);
      
      for(i=0; i<count; i++)
      {
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
      if (Input.ParameterOK("StartType","ClosedLoopFirstSolution"))
      {
         Vector clfs(Input.getArrayLength("StartType","ClosedLoopFirstSolution"));
         Input.getVector(clfs,"StartType","ClosedLoopFirstSolution");
         FirstSolution_ = Restrict_->RestrictDOF(clfs);
      }
      else
      {
         if (Input.ParameterOK("StartType","ClosedLoopUseAsFirst"))
         {
            ClosedLoopUseAsFirst_ = Input.getPosInt("StartType","ClosedLoopUseAsFirst");
            if (ClosedLoopUseAsFirst_ >= ClosedLoopStart_)
            {
               cerr << "Error: NewtonPCSolution -- ClosedLoopUseAsFirst must be < ClosedLoopStart."
                    << endl;
               exit(-33);
            }
         }
         else
         {
            // Default Value
            ClosedLoopUseAsFirst_ = Input.useInt(0,"StartType","ClosedLoopUseAsFirst");
            FirstSolution_ = one;
         }
      }

      Restrict_->SetDOF(one);
      
      Matrix Q(count, count);
      Matrix R(count, count_minus_one);
      
      QR(Restrict_->Stiffness(),Q,R,1);
      double tansign = 1.0;
      for (i=0;i<count_minus_one;++i)
      {
         tansign *= R[i][i]/fabs(R[i][i]);
      }
      for(i=0;i<count;i++)
      {
         Tangent1_[i] = Tangent2_[i] = Direction_ * tansign * Q[i][count_minus_one];
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
   ++counter_[3];
   
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
   ++counter_[4];
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
      
      if (CurrentDS_ < MinDS_)
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
      cout << " \tForceNorm = " << forcenorm;
      
      QR(Stiff_static, Q_static, R_static, 1);
      
      double tansign = 1.0;
      for (i=0;i<count_minus_one;++i)
      {
         tansign *= R_static[i][i]/fabs(R_static[i][i]);
      }
      for(i=0;i<count;i++)
      {
         Tangent2_[i] = Direction_*tansign*Q_static[i][count_minus_one];
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

      // solve and update
      MoorePenrose(Q_static,R_static, Force_static, Corrector_static);
      Magnitude1 = Corrector_static.Norm();
      Magnitude2 = Magnitude1;
      for (i=0;i<count;i++)
      {
         w_static[i] = v_static[i] - Corrector_static[i];
         difference_static[i] = -Corrector_static[i];
      }
      Restrict_->SetDOF(w_static);
      Force_static = Restrict_->Force();
      forcenorm = Force_static.Norm();

      corrections++;
      cout << " \tCorrectorNorm = " << Magnitude2 << " \tForceNorm = " << forcenorm;
      if (Magnitude2 > delta_max_)
      {
         f = accel_max_;
         CurrentDS_ /= f;
         cout << "\nPrediction " << predictions << " Corrector Iterations: "
              << corrections << "\n";
         ++predictions;
         continue;
      }
      f = max(f,sqrt(Magnitude2/delta_max_)*accel_max_);
      cout << " \tAcceleration = " << 1.0/f << "\n";

      // accpet iteration
      for(i=0;i<count;i++)
      {
         v_static[i] = w_static[i];
      }

      //test for iteration 1 convergence
      Converge_Test = 0;
      switch (ConvergeType_)
      {
         case Both:
            if ((forcenorm <= Converge_) && (Magnitude2 <= Converge_))
            {
               Converge_Test = 1;
            }
            break;
         case Force:
            if (forcenorm <= Converge_)
            {
               Converge_Test = 1;
            }
            break;
         case Displacement:
            if (Magnitude2 <= Converge_)
            {
               Converge_Test = 1;
            }
            break;
      }
      if (Converge_Test == 1)
      {
         // keep track of DS used to find the current point
         PreviousDS_ = CurrentDS_;
         
         // adaptively change DS for next point
         CurrentDS_ = CurrentDS_/f;
         if(CurrentDS_ > MaxDS_)
         {
            CurrentDS_ = MaxDS_;
         }
      }
      //CORRECTOR LOOP STARTS HERE (iteration 2)
      do
      {
         GetQR(Force_static,difference_static,Q_static,R_static);
         MoorePenrose(Q_static,R_static, Force_static,Corrector_static);
         Magnitude2 = Corrector_static.Norm();
         cout << " \tCorrectorNorm = " << Magnitude2;
         // update
         for (i=0;i<count;i++)
         {
            w_static[i] = v_static[i] - Corrector_static[i];
            difference_static[i] = -Corrector_static[i];
         }
         Restrict_->SetDOF(w_static);
         Force_static = Restrict_->Force();
         forcenorm = Force_static.Norm();
         cout << " \tForceNorm = " << forcenorm;
         
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

         switch (ConvergeType_)
         {
            case Both:
               if ((forcenorm <= Converge_) && (Magnitude2 <= Converge_))
               {
                  Converge_Test = 1;
               }
               break;
            case Force:
               if (forcenorm <= Converge_)
               {
                  Converge_Test = 1;
               }
               break;
            case Displacement:
               if (Magnitude2 <= Converge_)
               {
                  Converge_Test = 1;
               }
               break;
         }
         if (Converge_Test == 1)
         {
            // keep track of DS used to find the current point
            PreviousDS_ = CurrentDS_;
            
            // adaptively change DS for next point
            CurrentDS_ = CurrentDS_/f;
            if(CurrentDS_ > MaxDS_)
            {
               CurrentDS_ = MaxDS_;
            }
         }
         
         ++corrections;
      }
      while (Converge_Test != 1);
      cout << "Prediction " << predictions << " Corrector Iterations: " << corrections << "\n";
      ++predictions;
   }
   while (f >= accel_max_);
   
   if (Dot <= 0.0)
   {
      Omega_ = -Omega_;
   }
   
   if (BifStartFlag_)
   {
      v_static = Omega_*Tangent2_;
      v_static[count-1] = 0.0;
      v_static /= v_static.Norm();
      cout << "Projection on BifTangent = " << Restrict_->DOF()*BifTangent_
           << ",     Angle (deg.) with BifTangent = "
           << acos(v_static*BifTangent_)*(57.2957795130823) << "\n";
   }
   
   cout << "Converged with CorrectorNorm = " << Magnitude2 
        << ",     ForceNorm = " << forcenorm << "\n";
   
   if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Restrict_->DOF() - FirstSolution_).Norm() < CurrentDS_))
   {
      // We are done -- set currentsolution to numsolutions
      cout << "NOTE: Closed Loop detected at Solution # " << CurrentSolution_
           << " --- Terminating!" << "\n";
      
      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      CurrentSolution_++;
   }
   if (CurrentSolution_ == ClosedLoopUseAsFirst_)
   {
      // set First Solution for use with Closed Loop check.
      FirstSolution_ = Restrict_->DOF();
   }
   
   // always have current solution point printed
   good = 1;
   
   return good;
}

void NewtonPCSolution::FindCriticalPoint(Lattice* const Lat,int& TotalNumCPCrossings,
                                         PerlInput const& Input,int const& Width,ostream& out)
{
   ++counter_[5];
   
   int sz=PreviousSolution_.Dim();
   Vector tmp_diff(sz),tmp_DOF(Restrict_->DOF());
   double tmp_ds=0.0;
   for (int i=0;i<sz;++i)
   {
      tmp_ds += (PreviousSolution_[i]-tmp_DOF[i])*(PreviousSolution_[i]-tmp_DOF[i]);
   }
   tmp_ds = sqrt(tmp_ds);
   
   //ArcLengthSolution S1(Restrict_,Input,PreviousSolution_,Restrict_->DOF(),1);
   int MaxIter = 20;
   ArcLengthSolution::ConvergeType CT = ArcLengthSolution::Both; // initialize to avoid compiler complaint
   switch (ConvergeType_)
   {
      case Both:
         CT = ArcLengthSolution::Both;
         break;
      case Force:
         CT = ArcLengthSolution::Force;
         break;
      case Displacement:
         CT = ArcLengthSolution::Displacement;
         break;
   }
   ArcLengthSolution S1(Restrict_,Restrict_->DOF(),MaxIter,Converge_,CT,tmp_ds,tmp_ds,
                        tmp_ds,1.0,0.5,1.0,1,0,PreviousSolution_,
                        Restrict_->DOF()-PreviousSolution_,0,Vector(),10,0,-1,Echo_);
   S1.FindCriticalPoint(Lat,TotalNumCPCrossings,Input,Width,out);
   
   // Check to see if we should stop
   if ((StopAtCPCrossingNum_ > -1) && (TotalNumCPCrossings >= StopAtCPCrossingNum_))
      CurrentSolution_ = NumSolutions_;
}

void NewtonPCSolution::MoorePenrose(Matrix const& Q,Matrix const& R,Vector const& Force,
                                    Vector& Corrector) const
{
   ++counter_[1];
   
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
   ++counter_[0];
   
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
               // special form of Broyden update; simplified for Newton case
               Stiff_static[i][j] = Stiff_static[i][j] + (1/temp)*(Force[i] * diff[j]);
            }
         }
         
         QR(Stiff_static, Q, R, 1);  // Stiff_static^T = Q*R
         break;
      case Exact:
         // Stiffness^T = Q*R
         QR(Restrict_->Stiffness(),Q,R,1);
         break;
      default:
         cerr << "Unknown Update Type in NewtonPCSolution\n";
         exit(-20);
         break;
   }
}

// Remember A^T = Q*R
void NewtonPCSolution::UpdateQR(Vector const& Force,Vector const& difference,Matrix& QBar,
                                Matrix& RBar) const
{
   ++counter_[2];
   
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
   // modified appropratiely for A^T = Q*R
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
            A1 = C * QBar[k][i] + S * QBar[k][i+1];
            A2 = C * QBar[k][i+1] - S * QBar[k][i];
            QBar[k][i] = A1;
            QBar[k][i+1] = A2;
         }
         A1 = C * QBar[count][i] + S * QBar[count][i+1];
         A2 = C * QBar[count][i+1] - S * QBar[count][i];
         QBar[count][i] = A1;
         QBar[count][i+1] = A2;
         
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
            A1 = C * QBar[k][i] + S * QBar[k][i+1];
            A2 = C * QBar[k][i+1] - S * QBar[k][i];
            QBar[k][i] = A1;
            QBar[k][i+1] = A2;
         }
         A1 = C * QBar[count][i] + S * QBar[count][i+1];
         A2 = C * QBar[count][i+1] - S * QBar[count][i];
         QBar[count][i] = A1;
         QBar[count][i+1] = A2;
      }
   }
}
