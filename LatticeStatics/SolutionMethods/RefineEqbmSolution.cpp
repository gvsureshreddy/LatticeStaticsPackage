#include <fstream>
#include "RefineEqbmSolution.h"
#include "Matrix.h"

using namespace std;

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict, Vector const& one,
                                       double const& Converge, ConvergeType CnvrgTyp,
                                       int const& Echo) :
   Restrict_(Restrict),
   Echo_(Echo),
   SolutionFound_(0),
   NumSolutions_(1),
   Converge_(Converge),
   ConvergeType_(CnvrgTyp)
{
   Guesses_ = new Vector(one.Dim());
   Guesses_[0] = one;
}

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict, PerlInput const& Input,
                                       Vector const& one, int const& Echo) :
   Restrict_(Restrict),
   Echo_(Echo),
   SolutionFound_(0),
   NumSolutions_(1)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "RefineEqbmSolution");
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


   Guesses_ = new Vector(one.Dim());
   Guesses_[0] = one;
}

RefineEqbmSolution::RefineEqbmSolution(Restriction* const Restrict, PerlInput const& Input,
                                       int const& Echo) :
   Restrict_(Restrict),
   Echo_(Echo),
   SolutionFound_(0)
{
   // get needed parameters
   PerlInput::HashStruct Hash = Input.getHash("SolutionMethod", "RefineEqbmSolution");
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

   if (Input.ParameterOK(Hash, "Solution"))
   {
      int r, c;
      r = Input.getArrayLength(Hash, "Solution");
      c = Input.getArrayLength(Hash, "Solution", 0);
      Guesses_ = new Vector[r];
      NumSolutions_ = r;
      for (int i = 0; i < r; ++i)
      {
         Guesses_[i].Resize(c);
         Input.getVector(Guesses_[i], Hash, "Solution", i);
      }
   }
   else
   {
      NumSolutions_ = 1;
      Guesses_ = new Vector[1];
      Guesses_[0].Resize(Restrict_->DOF().Dim());
      Guesses_[0] = Restrict_->DOF();
   }
}

int RefineEqbmSolution::FindNextSolution(PerlInput const& Input, int const& Width, ostream& out)
{
   // set dofs to next guess
   if (SolutionFound_ < NumSolutions_)
   {
      Restrict_->SetDOF(Restrict_->RestrictDOF(Guesses_[SolutionFound_]));
   }
   else
   {
      cerr << "RefineEqbmSolution::FindNextSolution() called too many times.\n";
      exit(-31);
   }
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
   cout << "Corrector Iterations: " << itr << "\n"
        << "Converged with CorrectorNorm = " << dxnorm << ",     ForceNorm = " << forcenorm
        << "\n";

   ++SolutionFound_;

   return 1;
}
