#include "CBKinematics.h"

void CBKinematics:: InfluenceRegion(double *InfluenceRegion)
{
   Matrix Eigvals(1,DIM3);
   double tmp;
   // find largest eigenvalue of the inverse transformation
   // (i.e. from current to ref) and use influence cube of
   // that size...
   //
   // Use the fact that eigs of Uinv = 1/ eigs of U.
   //
   // Use U_*RefLattice_ as def grad.  This takes an
   // orthonormal lattice to the current config.
   // Thus, allowing non-square unit cells....
   //
   // Use F*F^T and take sqrt of eigvecs.
   Eigvals = SymEigVal(F_*(*RefLattice_)*((F_*(*RefLattice_)).Transpose()));
   tmp = sqrt(Eigvals[0][0]);
   for (int i=0;i<3;i++)
      if (sqrt(Eigvals[0][i]) < tmp) tmp = sqrt(Eigvals[0][i]);
   
   // Set to inverse eigenvalue
   tmp = 1.0/tmp;
   for (int i=0;i<DIM3;i++)
   {
      InfluenceRegion[i]=tmp;
   }
}
