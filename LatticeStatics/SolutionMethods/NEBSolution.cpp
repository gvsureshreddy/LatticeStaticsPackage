#include <fstream>
#include "NEBSolution.h"
#include "Matrix.h"

using namespace std;

NEBSolution::NEBSolution(Restriction* const Restrict, PerlInput const& Input,
                         int const& Echo) :
   Restrict_(Restrict),
   Echo_(Echo),
   SolutionFound_(0),
   NumSolutions_(1)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "NEBSolution");
   Converge_ = Input.getDouble(Hash, "ConvergeCriteria");

   if (Input.ParameterOK(Hash, "ConvergeType"))
   {
      char const* const cnvrgtyp = Input.getString(Hash, "ConvergeType");
      if (!strcmp("Both", cnvrgtyp))
      {
         ConvergeType_ = Both;
      }
      else if (!strcmp("Force", cnvrgtyp))
      {
         ConvergeType_ = Force;
      }
      else if (!strcmp("Displacement", cnvrgtyp))
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
      Input.useString("Both", Hash, "ConvergeType");  // Default Value
      ConvergeType_ = Both;
   }

   Width_ = Input.getInt("Main", "FieldWidth");

   SpringK_ = Input.getDouble(Hash, "SpringK");
   QDMass_ = Input.getDouble(Hash, "QDMass");
   Deltat_ = Input.getDouble(Hash, "dt");
   tFinal_ = Input.getDouble(Hash, "tFinal");

   DimDOFS_ = (Restrict_->DOF()).Dim();

   if (Input.ParameterOK(Hash, "IntermediateStateMatrix"))
   {
      NumInterStates_ = Input.getArrayLength(Hash, "IntermediateStateMatrix");
      InterStates_.Resize(NumInterStates_, DimDOFS_);
      Input.getMatrix(InterStates_, Hash, "IntermediateStateMatrix");
   }
   else
   {
      NumInterStates_ = 0;
   }

   NumReplicas_ = Input.getPosInt(Hash, "NumReplicas");

   if (Input.ParameterOK(Hash, "VelocityScaling"))
   {
      VScaling_ = Input.getDouble(Hash, "VelocityScaling");
   }
   else
   {
      cout << "Warning, Velocity Scaling not in Input file, using default value of zero initial velocity" << "\n";
      VScaling_ = 0.0;
   }

   if (VScaling_ != 0.0)
   {
      if (Input.ParameterOK(Hash, "RandomSeed"))
      {
         VSeed_ = Input.getPosInt(Hash, "RandomSeed");
      }
      else
      {
         cerr << "RandomSeed not in input file. Option:(1) Set VelocityScaling to 0.0 or, (2) input RandomSeed in input file. Exiting. \n";
         exit(-22);
      }
   }

   InitState_.Resize(DimDOFS_);
   FinalState_.Resize(DimDOFS_);

   Input.getVector(InitState_, Hash, "InitialState");
   Input.getVector(FinalState_, Hash, "FinalState");

   if ((InitState_.Dim()) != DimDOFS_)
   {
      cerr << "Dimension error in initial state in Input file. Exiting \n";
      exit(-22);
   }
   if ((FinalState_.Dim()) != DimDOFS_)
   {
      cerr << "Dimension error in initial state in Input file. Exiting \n";
      exit(-22);
   }

   TempState_.Resize(DimDOFS_);
   DOF_.Resize(DimDOFS_);
   InitDOF_.Resize(DimDOFS_);
   FinalDOF_.Resize(DimDOFS_);

   DimReplicas_ = DimDOFS_ - 1;

   iTangent_.Resize(DimReplicas_);
   FiTot_.Resize(DimReplicas_);

   QuenchedState_.Resize(NumReplicas_, DimReplicas_);

   rReplica1_.Resize(NumReplicas_, DimReplicas_);
   rReplica2_.Resize(NumReplicas_, DimReplicas_);
   vReplica1_.Resize(NumReplicas_, DimReplicas_);
   vReplica2_.Resize(NumReplicas_, DimReplicas_);


   if (VScaling_ == 0.0)
   {
      for (int i = 0; i < NumReplicas_; ++i)
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            vReplica1_[i][j] = 0.0;
         }
      }
   }
   else
   {
      srand(VSeed_);

      for (int i = 0; i < NumReplicas_; ++i)
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            vReplica1_[i][j] = (2 * (double) rand() / (double) RAND_MAX - 1) / VScaling_;
         }
      }
   }

   if (NumInterStates_ == 0)
   {
      EnergyBarrier_.Resize(2 + NumReplicas_, 0.0);
      StateMatrix_.Resize(2 + NumReplicas_, DimDOFS_, 0.0);
   }
   else
   {
      EnergyBarrier_.Resize(2 + NumInterStates_ + (NumReplicas_ * (NumInterStates_ + 1)), 0.0);
      StateMatrix_.Resize(2 + NumInterStates_ + (NumReplicas_ * (NumInterStates_ + 1)), DimDOFS_, 0.0);
   }
}

int NEBSolution::FindNextSolution(PerlInput const& Input, int const& Width, ostream& out)
{
   // set dofs to next guess
   if (SolutionFound_ < NumSolutions_)
   {
      Restrict_->SetDOF(Restrict_->RestrictDOF(InitState_));
   }
   else
   {
      cerr << "NEBSolution::FindNextSolution() called too many times.\n";
      return 0;
   }

   // Obtains Energies of InitStates and FinalStates
   // and stores it in EnergyBarrier_[NumReplicas_] and  EnergyBarrier_[NumReplicas_ + 1] respectively
   InitDOF_ = RefineState(InitState_);
   Restrict_->SetDOF(InitDOF_);
   EnergyBarrier_[0] = Restrict_->Energy();

   FinalDOF_ = RefineState(FinalState_);
   Restrict_->SetDOF(FinalDOF_);
   EnergyBarrier_[1] = Restrict_->Energy();
   cout << "InitDOF Energy = " << EnergyBarrier_[0] << "\n";
   cout << "FinalDOF Energy = " << EnergyBarrier_[1] << "\n";

   if ((InitDOF_ - InitState_).Norm() > 1000 * Converge_)
   {
      cout << "\n WARNING INITIAL DOF FAR AWAY FROM INPUT POINT \n" << "\n";
   }
   if ((FinalDOF_ - FinalState_).Norm() > 1000 * Converge_)
   {
      cout << "\n WARNING FINAL DOF FAR AWAY FROM INPUT POINT \n" << "\n";
   }

   for (int cols = 0; cols < DimDOFS_; ++cols)
   {
      StateMatrix_[0][cols] = InitDOF_[cols];
      StateMatrix_[1][cols] = FinalDOF_[cols];
   }

   // Takes care of replicas
   Replicas_.Resize(NumReplicas_, DimDOFS_);

   Matrix ReplicaTemp(NumReplicas_, DimReplicas_);

   if (NumInterStates_ == 0)     // If there are no intermediate states
   {
      // sets state of replicas
      for (int rows = 0; rows < NumReplicas_; ++rows)
      {
         for (int cols = 0; cols < DimDOFS_; ++cols)
         {
            Replicas_[rows][cols] = InitDOF_[cols] + ((rows + 1.0) / (NumReplicas_ + 1.0)) * (FinalDOF_[cols] - InitDOF_[cols]);
         }
      }

      ReplicaTemp = NEBQuenchedDynamics(InitDOF_, FinalDOF_, Replicas_, Deltat_, tFinal_, QDMass_);
      for (int rows = 0; rows < NumReplicas_; ++rows)
      {
         for (int cols = 0; cols < DimReplicas_; ++cols)
         {
            TempState_[cols] = ReplicaTemp[rows][cols];
            StateMatrix_[rows + 2][cols] = TempState_[cols];
         }
         TempState_[DimReplicas_] = InitDOF_[DimReplicas_];
         StateMatrix_[rows + 2][DimReplicas_] = InitDOF_[DimReplicas_];

         Restrict_->SetDOF(TempState_);
         EnergyBarrier_[rows + 2] = Restrict_->Energy();
      }
   }
   else      // If there are Intermediate States
   {
      // Sets Rows [2, NumInterStates_ - 1] for intermediate states
      for (int rows = 0; rows < NumInterStates_; ++rows)
      {
         for (int cols = 0; cols < DimDOFS_; ++cols)
         {
            TempState_[cols] = InterStates_[rows][cols];
            StateMatrix_[rows + 2][cols] = TempState_[cols];
         }
         Restrict_->SetDOF(TempState_);
         EnergyBarrier_[rows + 2] = Restrict_->Energy();
      }

      if (NumInterStates_ == 1)           // If there is only one intermediate state
      {
         int count2 = 0;
         Vector InitDOFTemp = InitDOF_;
         Vector FinalDOFTemp = FinalDOF_;
         for (int count = 0; count < NumInterStates_ + 1; ++count)
         {
            if (count == 0)
            {
               InitDOF_ = InitDOFTemp;
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  FinalDOF_[cols] = InterStates_[0][cols];
               }
            }
            else
            {
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  InitDOF_[cols] = InterStates_[0][cols];
               }
               FinalDOF_ = FinalDOFTemp;
            }
            // sets state of replicas
            for (int rows = 0; rows < NumReplicas_; ++rows)
            {
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  Replicas_[rows][cols] = InitDOF_[cols] + ((rows + 1.0) / (NumReplicas_ + 1.0)) * (FinalDOF_[cols] - InitDOF_[cols]);
               }
            }
            ReplicaTemp = NEBQuenchedDynamics(InitDOF_, FinalDOF_, Replicas_, Deltat_, tFinal_, QDMass_);

            for (int rows = 0; rows < NumReplicas_; ++rows)
            {
               for (int cols = 0; cols < DimReplicas_; ++cols)
               {
                  TempState_[cols] = ReplicaTemp[rows][cols];
                  StateMatrix_[count2 + NumInterStates_ + 2][cols] = TempState_[cols];
               }
               TempState_[DimReplicas_] = InitDOF_[DimReplicas_];
               StateMatrix_[count2 + NumInterStates_ + 2][DimReplicas_] = TempState_[DimReplicas_];

               Restrict_->SetDOF(TempState_);
               EnergyBarrier_[count2 + NumInterStates_ + 2] = Restrict_->Energy();
               count2++;
            }
         }
      }
      else           // If there are more than one intermediate state
      {
         int count2 = 0;
         Vector InitDOFTemp = InitDOF_;
         Vector FinalDOFTemp = FinalDOF_;
         for (int count = 0; count < NumInterStates_ + 1; ++count)
         {
            if (count == 0)
            {
               InitDOF_ = InitDOFTemp;
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  FinalDOF_[cols] = InterStates_[count][cols];
               }
            }
            else if (count == NumInterStates_)
            {
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  InitDOF_[cols] = InterStates_[count - 1][cols];
               }
               FinalDOF_ = FinalDOFTemp;
            }
            else
            {
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  InitDOF_[cols] = InterStates_[count - 1][cols];
                  FinalDOF_[cols] = InterStates_[count][cols];
               }
            }

            // sets state of replicas
            for (int rows = 0; rows < NumReplicas_; ++rows)
            {
               for (int cols = 0; cols < DimDOFS_; ++cols)
               {
                  Replicas_[rows][cols] = InitDOF_[cols] + ((rows + 1.0) / (NumReplicas_ + 1.0)) * (FinalDOF_[cols] - InitDOF_[cols]);
               }
            }
            ReplicaTemp = NEBQuenchedDynamics(InitDOF_, FinalDOF_, Replicas_, Deltat_, tFinal_, QDMass_);

            for (int rows = 0; rows < NumReplicas_; ++rows)
            {
               for (int cols = 0; cols < DimReplicas_; ++cols)
               {
                  TempState_[cols] = ReplicaTemp[rows][cols];
                  StateMatrix_[count2 + NumInterStates_ + 2][cols] = TempState_[cols];
               }
               TempState_[DimReplicas_] = InitDOF_[DimReplicas_];
               StateMatrix_[count2 + NumInterStates_ + 2][DimReplicas_] = TempState_[DimReplicas_];

               Restrict_->SetDOF(TempState_);
               EnergyBarrier_[count2 + NumInterStates_ + 2] = Restrict_->Energy();
               count2++;
            }
         }
      }
   }


   ++SolutionFound_;

   cout << " \n Solutions Found = " << SolutionFound_ << "\n";

   out << "\n \n ==========================NEB OUTPUT========================= \n";
   out << "State Matrix" << setw(Width) << StateMatrix_ << "\n";
   
   Vector Temp(2);
   Temp[0] = EnergyBarrier_[0];
   Temp[1] = EnergyBarrier_[1];		   
   long double MinValue = Temp.MinElement();
   
   if (Echo_)
   {
      cout << "StateMatrix" << setw(Width) << StateMatrix_ << "\n";
      cout << "Energy Vector: \n" << setw(Width) << EnergyBarrier_ << "\n";
   }
   
   out << "Energy Vector: \n" << setw(Width) << EnergyBarrier_ << "\n";

   Vector EBrel = EnergyBarrier_;
   for (int i = 0; i < EBrel.Dim(); ++i)
   {
      EBrel[i] -= MinValue;
   }
   out << "Normalized Energy Vector: \n" << setw(Width) << EBrel << "\n";
   out << "Energy Values offset by " << MinValue << "\n";
   
   if (Echo_)
   {
      cout << "Normalized Energy Vector: \n" << setw(Width) << EBrel << "\n";
      cout << "Energy Values offset by " << MinValue << "\n";			      
   }

   return 1;
}

Vector const& NEBSolution::RefineState(Vector const& CurrentState)
{
   Restrict_->SetDOF(CurrentState);

   Vector DOF = Restrict_->DOF();
   Vector dx(Restrict_->Force().Dim(), 0.0);
   Vector Stress = Restrict_->Force();
   Matrix Stiff(dx.Dim(), dx.Dim(), 0.0);
   Matrix tmpStiff(dx.Dim(), dx.Dim() + 1, 0.0);
   int itr = 0;
   int Converged = 0;
   // dxnorm initial value: should indicate a problem if this value is ever printed out...
   double dxnorm = -1.0;
   double forcenorm = Stress.Norm();

   const int MaxItr = 20;

   cout << "ForceNorm = " << forcenorm << "\n";

   while ((itr < MaxItr) && (!Converged))
   {
      ++itr;

      tmpStiff = Restrict_->Stiffness();
      for (int i = 0; i < dx.Dim(); ++i)
      {
         for (int j = 0; j < dx.Dim(); ++j)
         {
            Stiff[i][j] = tmpStiff[i][j];
         }
      }
#ifdef SOLVE_SVD
      dx = SolveSVD(Stiff, Stress, MAXCONDITION, Echo_);
#else
      dx = SolvePLU(Stiff, Stress);
#endif

      for (int i = 0; i < dx.Dim(); ++i)
      {
         DOF[i] -= dx[i];
      }
      Restrict_->SetDOF(DOF);
      Stress = Restrict_->Force();

      dxnorm = dx.Norm();
      forcenorm = Stress.Norm();
      cout << "\tCorrectorNorm = " << dxnorm
           << " \tForceNorm = " << forcenorm << "\n";

      switch (ConvergeType_)
      {
         case Both:
            if ((forcenorm <= Converge_) && (dxnorm <= Converge_))
            {
               Converged = 1;
            }
            break;
         case Force:
            if (forcenorm <= Converge_)
            {
               Converged = 1;
            }
            break;
         case Displacement:
            if (dxnorm <= Converge_)
            {
               Converged = 1;
            }
            break;
      }
   }

   DOF_ = Restrict_->DOF();

   return DOF_;
}

// Computes force on i-th replica given its nearest neighbors iMinusReplica and iPlusReplica
Vector const& NEBSolution::NEBForce(Vector const& iMinusReplica, Vector const& iPlusReplica, Vector const& iReplica)
{
   Vector FiPot(DimReplicas_), FiSpring(DimReplicas_);
   Vector iReplicaTemp = InitDOF_;

   Vector iTangent = NEBTangent(iMinusReplica, iPlusReplica, iReplica);

   for (int i = 0; i < DimReplicas_; ++i)
   {
      iReplicaTemp[i] = iReplica[i];
   }

   Restrict_->SetDOF(iReplicaTemp);

   Vector StateForce = -1.0 * (Restrict_->Force());
   FiPot = StateForce - (StateForce * iTangent) * iTangent;
   FiSpring = (SpringK_ * (((iPlusReplica - iReplica).Norm()) - ((iReplica - iMinusReplica).Norm()))) * iTangent;

   FiTot_ = FiPot + FiSpring;
   return FiTot_;
}

// Obtains the tangent vector iTangent
Vector const& NEBSolution::NEBTangent(Vector const& iMinusReplica, Vector const& iPlusReplica, Vector const& iReplica)
{
   double EnergyiPlus, EnergyiMinus, Energyi, DeltaEMax, DeltaEMin;
   int casetype = 0;
   Vector Temp(2), iTangentTemp(DimReplicas_), TangentiPlus(DimReplicas_), TangentiMinus(DimReplicas_);
   Vector iReplicaTemp = InitDOF_;

   TangentiPlus = iPlusReplica - iReplica;
   TangentiMinus = iReplica - iMinusReplica;

   for (int i = 0; i < DimReplicas_; ++i)
   {
      iReplicaTemp[i] = iReplica[i];
   }
   Restrict_->SetDOF(iReplicaTemp);
   Energyi = Restrict_->Energy();

   for (int i = 0; i < DimReplicas_; ++i)
   {
      iReplicaTemp[i] = iPlusReplica[i];
   }
   Restrict_->SetDOF(iReplicaTemp);
   EnergyiPlus = Restrict_->Energy();

   for (int i = 0; i < DimReplicas_; ++i)
   {
      iReplicaTemp[i] = iMinusReplica[i];
   }
   Restrict_->SetDOF(iReplicaTemp);
   EnergyiMinus = Restrict_->Energy();

   if (((EnergyiMinus >= Energyi) && (Energyi >= EnergyiPlus)) || ((EnergyiPlus >= Energyi) && (Energyi >= EnergyiMinus)))
   {
      casetype = 0;
   }
   else
   {
      casetype = 1;
   }

   if (casetype == 0)
   {
      if ((EnergyiMinus <= Energyi) && (Energyi <= EnergyiPlus))
      {
         iTangentTemp = TangentiPlus;
      }
      else
      {
         iTangentTemp = TangentiMinus;
      }
   }
   else
   {
      Temp[0] = abs(EnergyiPlus - Energyi);
      Temp[1] = abs(EnergyiMinus - Energyi);

      DeltaEMax = Temp.MaxElement();
      DeltaEMin = Temp.MinElement();

      if (EnergyiMinus <= EnergyiPlus)
      {
         iTangentTemp = TangentiPlus * DeltaEMax + TangentiMinus * DeltaEMin;
      }
      else
      {
         iTangentTemp = TangentiPlus * DeltaEMin + TangentiMinus * DeltaEMax;
      }
   }
   iTangent_ = iTangentTemp / iTangentTemp.Norm();
   return iTangent_;
}

Matrix const& NEBSolution::NEBQuenchedDynamics(Vector const& InitDOF, Vector const& FinalDOF, Matrix const& Replicas_, double const& Deltat, double const& tFinal, double const& M)
{
   double t = 0.0;

   for (int rows = 0; rows < NumReplicas_; ++rows)
   {
      for (int cols = 0; cols < DimReplicas_; ++cols)
      {
         rReplica1_[rows][cols] = Replicas_[rows][cols];
         rReplica2_[rows][cols] = Replicas_[rows][cols];
      }
   }

   while (t <= tFinal)
   {
      cout << "t = " << t << "\n";
      NEBQDStep(InitDOF, FinalDOF, Deltat, M);
      cout << "Normalized TotalReplicaForce_ = " << TotalReplicaForce_ / NumReplicas_ << "\n \n";
      t = t + Deltat;
   }

   for (int i = 0; i < NumReplicas_; ++i)
   {
      for (int j = 0; j < DimReplicas_; ++j)
      {
         QuenchedState_[i][j] = rReplica1_[i][j];
      }
   }
   return QuenchedState_;
}


Matrix const& NEBSolution::NEBQDStep(Vector const& InitDOF, Vector const& FinalDOF, double const& Deltat, double const& M)
{
   double vp;
   Vector Force(DimReplicas_), V(DimReplicas_), iPlus(DimReplicas_), iMinus(DimReplicas_), iReplica(DimReplicas_);

   for (int i = 0; i < NumReplicas_; ++i)
   {
      if (i == 0)
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = InitDOF[j];
            iReplica[j] = rReplica1_[i][j];
            iPlus[j] = rReplica1_[i + 1][j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         TotalReplicaForce_ = Force.Norm();
         for (int j = 0; j < DimReplicas_; ++j)
         {
            rReplica2_[i][j] = rReplica1_[i][j] + vReplica1_[i][j] * Deltat + 0.5 * Deltat * Deltat * Force[j] / M;
            vReplica2_[i][j] = vReplica1_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
      }
      else if (i == (NumReplicas_ - 1))
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = rReplica1_[i - 1][j];
            iReplica[j] = rReplica1_[i][j];
            iPlus[j] = FinalDOF[j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         TotalReplicaForce_ = TotalReplicaForce_ + Force.Norm();
         for (int j = 0; j < DimReplicas_; ++j)
         {
            rReplica2_[i][j] = rReplica1_[i][j] + vReplica1_[i][j] * Deltat + 0.5 * Deltat * Deltat * Force[j] / M;
            vReplica2_[i][j] = vReplica1_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
      }
      else if ((i > 0) && (i < (NumReplicas_ - 1)))
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = rReplica1_[i - 1][j];
            iReplica[j] = rReplica1_[i][j];
            iPlus[j] = rReplica1_[i + 1][j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         TotalReplicaForce_ = TotalReplicaForce_ + Force.Norm();
         for (int j = 0; j < DimReplicas_; ++j)
         {
            rReplica2_[i][j] = rReplica1_[i][j] + vReplica1_[i][j] * Deltat + 0.5 * Deltat * Deltat * Force[j] / M;
            vReplica2_[i][j] = vReplica1_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
      }
   }
   for (int i = 0; i < NumReplicas_; ++i)
   {
      if (i == 0)
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = InitDOF[j];
            iReplica[j] = rReplica2_[i][j];
            iPlus[j] = rReplica2_[i + 1][j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         for (int j = 0; j < DimReplicas_; ++j)
         {
            V[j] = vReplica2_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
         vp = V * Force;
         if (vp < 0)
         {
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = 0.0;
            }
         }
         else
         {
            V = vp * Force / Force.Norm();
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = V[j];
            }
         }
      }
      else if (i == (NumReplicas_ - 1))
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = rReplica2_[i - 1][j];
            iReplica[j] = rReplica2_[i][j];
            iPlus[j] = FinalDOF[j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         for (int j = 0; j < DimReplicas_; ++j)
         {
            V[j] = vReplica2_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
         vp = V * Force;
         if (vp < 0)
         {
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = 0.0;
            }
         }
         else
         {
            V = vp * Force / Force.Norm();
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = V[j];
            }
         }
      }
      else if ((i > 0) && (i < (NumReplicas_ - 1)))
      {
         for (int j = 0; j < DimReplicas_; ++j)
         {
            iMinus[j] = rReplica2_[i - 1][j];
            iReplica[j] = rReplica2_[i][j];
            iPlus[j] = rReplica2_[i + 1][j];
         }
         Force = NEBForce(iMinus, iPlus, iReplica);
         for (int j = 0; j < DimReplicas_; ++j)
         {
            V[j] = vReplica2_[i][j] + 0.5 * Deltat * Force[j] / M;
         }
         vp = V * Force;
         if (vp < 0)
         {
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = 0.0;
            }
         }
         else
         {
            V = vp * Force / Force.Norm();
            for (int j = 0; j < DimReplicas_; ++j)
            {
               vReplica1_[i][j] = V[j];
            }
         }
      }
   }

   for (int i = 0; i < NumReplicas_; ++i)
   {
      for (int j = 0; j < DimReplicas_; ++j)
      {
         rReplica1_[i][j] = rReplica2_[i][j];
      }
   }
   return rReplica1_;
}
