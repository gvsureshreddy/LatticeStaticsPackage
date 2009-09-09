#include <fstream>
#include "RefineEqbmSolution.h"
#include "Matrix.h"

using namespace std;

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict,Vector const& one,
                                       double const& Converge,ConvergeType CnvrgTyp,
                                       int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     SolutionFound_(0),
     Converge_(Converge),
     ConvergeType_(CnvrgTyp)
{
   Restrict_->SetDOF(one);
}

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict,PerlInput const& Input,
                                       Vector const& one,int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     SolutionFound_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","RefineEqbmSolution");
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

   Restrict_->SetDOF(one);
}

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict,PerlInput const& Input,
                                       int const& Echo)
   : Restrict_(Restrict),
     Echo_(Echo),
     SolutionFound_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","RefineEqbmSolution");
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

   if (Input.ParameterOK(Hash,"Solution"))
   {
      Vector onetmp(Input.getArrayLength(Hash,"Solution"));
      Input.getVector(onetmp,Hash,"Solution");
      Restrict_->SetDOF(Restrict_->RestrictDOF(onetmp));
   }
}

int RefineEqbmSolution::FindNextSolution()
{
   Vector DOF = Restrict_->DOF();
   Vector dx(Restrict_->Force().Dim(),0.0);
   Vector Stress=Restrict_->Force();
   Matrix Stiff(dx.Dim(),dx.Dim(),0.0);
   Matrix tmpStiff(dx.Dim(),dx.Dim()+1,0.0);
   int itr=0;
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
      for (int i=0;i<dx.Dim();++i)
         for (int j=0;j<dx.Dim();++j)
            Stiff[i][j] = tmpStiff[i][j];
#ifdef SOLVE_SVD
      dx = SolveSVD(Stiff,Stress,MAXCONDITION,Echo_);
#else
      dx = SolvePLU(Stiff,Stress);
#endif

      for (int i=0;i<dx.Dim();++i)
         DOF[i] -= dx[i];
      Restrict_->SetDOF(DOF);
      Stress=Restrict_->Force();

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
   cout << "Corrector Iterations: " << itr << "\n"
        << "Converged with CorrectorNorm = " << dxnorm << ",     ForceNorm = " << forcenorm
        << "\n";
   SolutionFound_ = 1;

   return 1;
}
