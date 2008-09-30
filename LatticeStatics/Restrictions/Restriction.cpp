#include "Restriction.h"

void Restriction::ConsistencyCheck(Vector const& dof,double const& ConsistencyEpsilon,
                                   int const& Width,ostream& out)
{
   double potential;
   int Dim=DOF().Dim();
   Matrix
      Stiff(Dim-1,Dim),
      PerturbedStiff(Dim-1,Dim);
   Vector
      Frc(Dim-1),
      PerturbedForce(Dim-1),
      pert(Dim,0.0),
      RHS(Dim-1);
   
   // Set state to dof
   SetDOF(dof);
   
   // Do Consistency check
   for (int i=0;i<70;i++) out << "="; out << "\n";
   out << "Restriction Consistency Check." << "\n";
   out << "Epsilon = " << ConsistencyEpsilon << "\n";
   out << "Force(dof) * Epsilon" << "\n";
   Frc = ConsistencyEpsilon*Force();
   out << setw(Width) << Frc << "\n";
   out << "Stiffness(dof) * Epsilon" << "\n";
   Stiff = ConsistencyEpsilon*Stiffness();
   out << setw(Width) << Stiff << "\n";
   for (int i=0;i<Dim;i++)
   {
      // Get RHS
      SetDOF(dof);
      potential = Energy();
      RHS = Force();
      
      // Perturb the lattice state
      pert.Resize(Dim,0.0);
      pert[i]=1.0;
      SetDOF(dof + ConsistencyEpsilon*pert);
      // Get Check
      potential = Energy() - potential;      
      if (i < Dim-1) PerturbedForce[i] = potential;
      RHS = Force() - RHS;
      for (int j=0;j<Dim-1;j++)
         PerturbedStiff[j][i] = RHS[j];
   }
   
   out << "Energy(dof) - Energy(dof + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedForce << "\n" << "\n";;
   out << "Force(dof) - Force(dof + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedStiff;
   out << "Difference" << "\n";
   out << setw(Width) << Frc - PerturbedForce << "\n";
   out << setw(Width) << Stiff - PerturbedStiff;
   
   for (int i=0;i<70;i++) out << "="; out << "\n";
}

