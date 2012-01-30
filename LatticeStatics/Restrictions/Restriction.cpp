#include "Restriction.h"

Restriction::~Restriction()
{
   if (SymmetryCheckCount_ > 0)
   {
      delete[] SymmetryCheck_;
   }
}

Restriction::Restriction(PerlInput const& Input) :
   Lattice_(NULL)
{
   PerlInput::HashStruct Hash = Input.getHash("Restriction");
   if (Input.ParameterOK(Hash, "SymmetryCheckProjectionMatrices"))
   {
      SymmetryCheckCount_ = Input.getArrayLength(Hash, "SymmetryCheckProjectionMatrices");
      if (SymmetryCheckCount_ == 0)
      {
         cerr << "Error. " << Name()
              << " SymmetryCheckProjectionMatrices is empty\n";
         exit(-37);
      }

      SymmetryCheck_ = new SparseMatrix[SymmetryCheckCount_];

      for (int i = 0; i < SymmetryCheckCount_; ++i)
      {
         int Rows = Input.getArrayLength(Hash, "SymmetryCheckProjectionMatrices", i);
         int Cols = Input.getArrayLength(Hash, "SymmetryCheckProjectionMatrices", i, 0);
         Matrix SCPM(Rows, Cols);
         Input.getMatrix(SCPM, Hash, "SymmetryCheckProjectionMatrices", i);
         int nononzero = 0;
         for (int j = 0; j < Rows; ++j)
         {
            for (int k = 0; k < Cols; ++k)
            {
               if (fabs(SCPM[j][k]) > 1.0e-15)
               {
                  ++nononzero;
               }
            }
         }

         SymmetryCheck_[i].Resize(Rows, Cols + 1, nononzero); // Cols+1 to ignore load value

         int count = 0;
         for (int j = 0; j < Rows; ++j)
         {
            for (int k = 0; k < Cols; ++k)
            {
               if (fabs(SCPM[j][k]) > 1.0e-15)
               {
                  SymmetryCheck_[i].SetNonZeroEntry(count, j, k, SCPM[j][k]);
                  ++count;
               }
            }
         }
      }
   }
   else
   {
      SymmetryCheckCount_ = 0;
   }

   if (Input.ParameterOK(Hash, "SymmetryCheckTolerance"))
   {
      SymmetryCheckTol_ = Input.getDouble(Hash, "SymmetryCheckTolerance");
   }
   else
   {
      // Default Value
      SymmetryCheckTol_ = Input.useDouble(1.0e-14, Hash, "SymmetryCheckTolerance");
   }
}

int Restriction::SymmetryOK() const
{
   int retval = SymmetryCheckCount_;
   for (int i = 0; i < SymmetryCheckCount_; ++i)
   {
      if (!((SymmetryCheck_[i] * DOF()).Norm() < SymmetryCheckTol_))
      {
         --retval;
      }
   }

   return !retval;
}

void Restriction::ConsistencyCheck(Vector const& dof, double const& ConsistencyEpsilon,
                                   int const& Width, ostream& out)
{
   double potential;
   int Dim = DOF().Dim();
   Matrix
   Stiff(Dim - 1, Dim),
   PerturbedStiff(Dim - 1, Dim);
   Vector
   Frc(Dim - 1),
   PerturbedForce(Dim - 1),
   pert(Dim, 0.0),
   RHS(Dim - 1);

   // Set state to dof
   SetDOF(dof);

   // Do Consistency check
   for (int i = 0; i < 70; i++)
   {
      out << "=";
   }
   out << "\n";
   out << "Restriction Consistency Check." << "\n";
   out << "Epsilon = " << ConsistencyEpsilon << "\n";
   out << "Force(dof) * Epsilon" << "\n";
   Frc = ConsistencyEpsilon * Force();
   out << setw(Width) << Frc << "\n";
   out << "Stiffness(dof) * Epsilon" << "\n";
   Stiff = ConsistencyEpsilon * Stiffness();
   out << setw(Width) << Stiff << "\n";
   for (int i = 0; i < Dim; i++)
   {
      // Get RHS
      SetDOF(dof);
      potential = Energy();
      RHS = Force();

      // Perturb the lattice state
      pert.Resize(Dim, 0.0);
      pert[i] = 1.0;
      SetDOF(dof + ConsistencyEpsilon * pert);
      // Get Check
      potential = Energy() - potential;
      if (i < Dim - 1)
      {
         PerturbedForce[i] = potential;
      }
      RHS = Force() - RHS;
      for (int j = 0; j < Dim - 1; j++)
      {
         PerturbedStiff[j][i] = RHS[j];
      }
   }

   out << "Energy(dof) - Energy(dof + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedForce << "\n" << "\n";
   out << "Force(dof) - Force(dof + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedStiff;
   out << "Difference" << "\n";
   out << setw(Width) << Frc - PerturbedForce << "\n";
   out << setw(Width) << Stiff - PerturbedStiff;

   for (int i = 0; i < 70; i++)
   {
      out << "=";
   }
   out << "\n";
}
