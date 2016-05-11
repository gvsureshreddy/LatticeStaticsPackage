#include "Restriction.h"
#include <math.h>
#include <cstdlib>

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

void Restriction::ConsistencyCheckRidders(Vector const& dof,double const& ConsistencyEpsilon,
                                   int const& Width, ostream& out)
{
   double potential;
   int Dim = DOF().Dim();
   int Ntab = 5;
   double con = 1.4;
   double fac;
   double safe = 2.0;
   Matrix
   Stiff(Dim - 1, Dim),
   PerturbedStiff(Dim - 1, Dim),
   Aa(Ntab,Ntab),
   Bb(Ntab,Ntab);
   Vector
   Frc(Dim - 1),
   PerturbedForce(Dim - 1),
   pert(Dim, 0.0),
   RHS(Dim - 1),
   A(Dim - 1),
   B(Dim - 1);
   double errt;
   double initial = 1.0e30;
   double err = initial;
   double eps = ConsistencyEpsilon;
   double a,b;


   for (int i = 0; i < 70; i++)
   {
      out << "=";
   }


   // Set state to dof
   SetDOF(dof);

   out << "\n" << setw(Width) << *Lattice_ << "\n";

   out << "\n";
   out << "Restriction Consistency Check." << "\n";
   out << "Epsilon = " << ConsistencyEpsilon << "\n";
   out << "Force(dof)" << "\n";
   Frc = Force();
   out << setw(Width) << Frc << "\n";
   out << "Stiffness(dof)" << "\n";
   Stiff = Stiffness();
   out << setw(Width) << Stiff << "\n";



   for (int k = 0; k < Dim; k++)
   {



      err = initial;
      // Get RHS
      SetDOF(dof);
      potential = Energy();
      RHS = Force();

      // Perturb the lattice state
      pert.Resize(Dim, 0.0);
      pert[k] = 1.0;

    if (k < Dim - 1)
    {
      SetDOF(dof + eps * pert);
      a = Energy();
      SetDOF(dof - eps * pert);
      b = Energy();


      Aa[0][0] = (a-b)/(2.0* eps);

      for (int i=1; i < Ntab; ++i)
      {


         eps /= con;

         SetDOF(dof + eps * pert);
         a = Energy();
         SetDOF(dof - eps * pert);
         b = Energy();

         Aa[0][i] = (a-b)/(2.0*eps);

         fac = con * con;
         for (int j = 1; j <= i; ++j)
         {

             Aa[j][i] = (Aa[j-1][i] * fac - Aa[j-1][i-1])/(fac - 1.0);
             fac = con * con * fac;
             errt = fmax( fabs(Aa[j][i] - Aa[j-1][i]), fabs( Aa[j][i] - Aa[j-1][i-1]));
             if (errt <= err )
             {
                err = errt;
                PerturbedForce[k]=Aa[j][i];
             }
         }
         if (fabs(Aa[i][i]-Aa[i-1][i-1])  >= ( safe * err )) break;
      }
     }


      for (int l =  0; l < Dim - 1; ++l)
      {

         err = initial;
         eps = ConsistencyEpsilon;

         SetDOF(dof + eps * pert);
         A = Force();
         SetDOF(dof - eps * pert);
         B = Force();

         Bb[0][0] = (A[l]-B[l])/(2.0*eps);

         for (int i=1; i < Ntab; ++i)
         {

            eps /= con;

            SetDOF(dof + eps * pert);
            A = Force();
            SetDOF(dof - eps * pert);
            B = Force();

            Bb[0][i] = (A[l]-B[l])/(2.0*eps);

            fac = con * con;
            for (int j = 1; j <= i; ++j)
            {

               Bb[j][i] = (Bb[j-1][i] * fac - Bb[j-1][i-1])/(fac - 1.0);
               fac = con * con * fac;
               errt = fmax( fabs(Bb[j][i] - Bb[j-1][i]), fabs( Bb[j][i] - Bb[j-1][i-1]));
               if (errt <= err )
               {
                  err = errt;
                  PerturbedStiff[l][k]=Bb[j][i];
               }
            }
         if (fabs(Bb[i][i]-Bb[i-1][i-1]) >= safe * err ) break;
         }


       }
   }

   out << "Energy Derivative" << "\n";
   out << setw(Width) << PerturbedForce << "\n" << "\n";
   out << "Force Derivative" << "\n";
   out << setw(Width) << PerturbedStiff;
   out << "Difference" << "\n";
   out << setw(Width) << Frc - PerturbedForce << "\n";
   out << setw(Width) << Stiff - PerturbedStiff << "\n";
   out << "\n";
   out << "ForceNorm = " << Frc.Norm() << ", Diff Norm = " << (Frc - PerturbedForce).Norm() << "\n";
//   out << "StiffnessNorm = " << Stiff.Norm() << ", Diff Norm = " << (Stiff - PerturbedStiff).Norm() << "\n";


   for (int i = 0; i < 70; i++)
   {
      out << "=";
   }
   out << "\n";

}
