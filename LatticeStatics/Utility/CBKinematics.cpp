#include "CBKinematics.h"

void CBKinematics::Reset()
{
   int i,q,p;
   U_[0][0] = (*DOF_)[0];
   U_[1][1] = (*DOF_)[1];
   U_[2][2] = (*DOF_)[2];
   U_[0][1] = U_[1][0] = (*DOF_)[3];
   U_[0][2] = U_[2][0] = (*DOF_)[4];
   U_[1][2] = U_[2][1] = (*DOF_)[5];
   
   S_[0][0] = 0.0;
   S_[0][1] = 0.0;
   S_[0][2] = 0.0;
   i=6;
   for (q=1;q<InternalAtoms_;++q)
   {
      for (p=0;p<3;p++)
      {
	 S_[q][p] = (*DOF_)[i++];
      }
   }
}

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
   Eigvals = SymEigVal(U_*(*RefLattice_)*((U_*(*RefLattice_)).Transpose()));
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
