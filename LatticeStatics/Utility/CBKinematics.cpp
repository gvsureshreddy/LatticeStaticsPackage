#include "CBKinematics.h"

#define LINELENGTH 600

int const CBKinematics::DIM3 = 3;

CBKinematics::CBKinematics(int const& InternalAtoms,Matrix& RefLattice,
                           Vector* const AtomPositions)
   : InternalAtoms_(InternalAtoms),
     RefLattice_(RefLattice)
{
   // Set AtomPositions_
   InternalPOS_ = new Vector[InternalAtoms_];
   for (int i=0;i<InternalAtoms_;++i)
   {
      InternalPOS_[i].Resize(DIM3);
      InternalPOS_[i] = AtomPositions[i];
   }
}

CBKinematics::CBKinematics(PerlInput const& Input,PerlInput::HashStruct const* const ParentHash)
{
   PerlInput::HashStruct Hash;
   if (ParentHash != 0)
   {
      Hash = Input.getHash(*ParentHash,"CBKinematics");
   }
   else
   {
      Hash = Input.getHash("CBKinematics");
   }

   // Set number of atoms in unit cell
   InternalAtoms_ = Input.getPosInt(Hash,"InternalAtoms");
   
   // Set RefLattice_
   RefLattice_.Resize(DIM3,DIM3);
   Input.getMatrix(RefLattice_,Hash,"LatticeBasis");
   
   // Set AtomPositions_
   InternalPOS_ = new Vector[InternalAtoms_];
   for (int i=0;i<InternalAtoms_;++i)
   {
      InternalPOS_[i].Resize(DIM3);
      Input.getVector(InternalPOS_[i],Hash,"AtomPositions",i);
   }
}

void CBKinematics::InfluenceRegion(double* const InfluenceRegion)
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
   Eigvals = SymEigVal(F_*RefLattice_*((F_*RefLattice_).Transpose()));
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

void CBKinematics::SetReferenceDOFs()
{
   DOF_.Resize(DOFS(),0.0);
   for (int i=0;i<DIM3;++i)
   {
      DOF_[INDF(i,i)] = 1.0;
   }
   Reset();
}

void CBKinematics::SetReferenceToCurrent()
{
   Matrix CurrentLattice(DIM3,DIM3,0.0);
   for (int i=0;i<DIM3;++i)
   {
      for (int j=0;j<DIM3;++j)
         for (int k=0;k<DIM3;++k)
         {
            CurrentLattice[i][j] += F_[j][k]*RefLattice_[i][k];
         }
   }
   RefLattice_ = CurrentLattice;
   
   for (int i=0;i<InternalAtoms_;++i)
   {
      for (int j=0;j<DIM3;++j)
      {
         InternalPOS_[i][j] = InternalPOS_[i][j] + S_[i][j];
      }
   }
   
   SetReferenceDOFs();
}

Vector CBKinematics::CurrentLatticeVec(int const& p) const
{
   Vector tmp(DIM3,0.0);
   
   for (int i=0;i<DIM3;++i)
   {
      for (int j=0;j<DIM3;++j)
         tmp[i] += F_[i][j]*RefLattice_[p][j];
   }
   
   return tmp;
}
