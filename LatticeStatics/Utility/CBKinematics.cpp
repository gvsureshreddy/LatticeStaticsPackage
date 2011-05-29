#include "CBKinematics.h"

#define LINELENGTH 600

int const CBKinematics::DIM3 = 3;

CBKinematics::CBKinematics(int const& InternalAtoms, Matrix& RefLattice,
                           Vector* const AtomPositions) :
   InternalAtoms_(InternalAtoms),
   RefLattice_(RefLattice)
{
   // Set AtomPositions_
   InternalPOS_ = new Vector[InternalAtoms_];
   for (int i = 0; i < InternalAtoms_; ++i)
   {
      InternalPOS_[i].Resize(DIM3);
      InternalPOS_[i] = AtomPositions[i];
   }
}

CBKinematics::CBKinematics(PerlInput const& Input, PerlInput::HashStruct const* const ParentHash)
{
   PerlInput::HashStruct Hash;
   if (ParentHash != 0)
   {
      Hash = Input.getHash(*ParentHash, "CBKinematics");
   }
   else
   {
      Hash = Input.getHash("CBKinematics");
   }

   // Set RefLattice_
   RefLattice_.Resize(DIM3, DIM3);
   Input.getMatrix(RefLattice_, Hash, "LatticeBasis"); // g[i] = RefLattice_[i][j]*e[j]

   // Set number of atoms in unit cell
   InternalAtoms_ = Input.getPosInt(Hash, "InternalAtoms");
   if (InternalAtoms_ > CBK_MAX_ATOMS)
   {
      cerr << "Error: CBK, InternalAtoms > CBK_MAX_ATOMS=" << CBK_MAX_ATOMS << ", exiting...\n";
      exit(-23);
   }

   // Set AtomPositions_
   InternalPOS_ = new Vector[InternalAtoms_];
   for (int i = 0; i < InternalAtoms_; ++i)
   {
      InternalPOS_[i].Resize(DIM3);
      Input.getVector(InternalPOS_[i], Hash, "AtomPositions", i);
   }

   // Set AtomSpecies_
   Input.getPosIntVector(AtomSpecies_, InternalAtoms_, Hash, "AtomSpecies");
   NumberofSpecies_ = AtomSpecies_[0];
   for (int i = 1; i < InternalAtoms_; ++i)
   {
      if (NumberofSpecies_ < AtomSpecies_[i])
      {
         NumberofSpecies_ = AtomSpecies_[i];
      }
   }
   NumberofSpecies_++;

   // Check for supercell specification
   if (Input.ParameterOK(Hash, "SuperCell"))
   {
      int mu[DIM3][DIM3];
      int latrange[DIM3][2];
      int TmpIntAtoms;
      // g+[i] = SuperCell[i][j]*g[j]
      Input.getIntMatrix(&(mu[0][0]), DIM3, DIM3, Hash, "SuperCell");
      Matrix Mu(DIM3, DIM3);
      Matrix MuInvT(DIM3, DIM3);
      Matrix TmpRefLat(DIM3, DIM3, 0.0);
      for (int i = 0; i < DIM3; ++i)
      {
         latrange[i][0] = 0;
         latrange[i][1] = 0;
         for (int j = 0; j < DIM3; ++j)
         {
            Mu[i][j] = double(mu[i][j]);

            if (mu[i][j] < 0)
            {
               if (mu[j][i] < latrange[i][0])
               {
                  latrange[i][0] = mu[j][i];
               }
            }
            else
            {
               if (mu[j][i] > latrange[i][1])
               {
                  latrange[i][1] = mu[j][i];
               }
            }
         }
      }
      MuInvT = (Mu.Inverse()).Transpose();

      // Find lattice vectors in supercell
      int det = int(Mu.Det());
      int cnt = 0;
      Vector* CellVecs = new Vector[det];
      Vector L(DIM3);
      Vector l(DIM3);
      for (int i = latrange[0][0]; i <= latrange[0][1]; ++i)
      {
         for (int j = latrange[1][0]; j <= latrange[1][1]; ++j)
         {
            for (int k = latrange[2][0]; k <= latrange[2][1]; ++k)
            {
               L[0] = i;
               L[1] = j;
               L[2] = k;

               l = MuInvT * L; // l[i]*g+[i] = l[i]*Mu[i][j]*g[j] = L[j]g[j]

               if (((l[0] >= 0.0) && (l[0] < 1.0)) &&
                   ((l[1] >= 0.0) && (l[1] < 1.0)) &&
                   ((l[2] >= 0.0) && (l[2] < 1.0)))
               {
                  CellVecs[cnt].Resize(DIM3);
                  CellVecs[cnt][0] = double(i);
                  CellVecs[cnt][1] = double(j);
                  CellVecs[cnt][2] = double(k);
                  cnt++;
               }
            }
         }
      }

      // Overwrite with new valuse and add to input file
      TmpRefLat = Mu * RefLattice_;
      Input.useMatrix(TmpRefLat, Hash, "LatticeBasis"); // should change this so it doesn't print 'default value'
      RefLattice_ = TmpRefLat;

      TmpIntAtoms = cnt * InternalAtoms_;
      Input.usePosInt(TmpIntAtoms, Hash, "InternalAtoms"); // see above note
      if (InternalAtoms_ > CBK_MAX_ATOMS)
      {
         cerr << "Error: CBK, InternalAtoms > CBK_MAX_ATOMS=" << CBK_MAX_ATOMS
              << ", exiting...\n";
         exit(-23);
      }

      Vector* TmpIntPOS = new Vector[TmpIntAtoms];
      for (int i = 0; i < cnt; ++i)
      {
         for (int j = 0; j < InternalAtoms_; ++j)
         {
            TmpIntPOS[i * InternalAtoms_ + j].Resize(DIM3);
            TmpIntPOS[i * InternalAtoms_ + j] = MuInvT * (CellVecs[i] + InternalPOS_[j]);
            for (int k = 0; k < DIM3; ++k)
            {
               // for the case where atom sits on boundary of cell,
               // make sure it is on "lower left"
               if (TmpIntPOS[i * InternalAtoms_ + j][k] >= 1.0)
               {
                  TmpIntPOS[i * InternalAtoms_ + j][k]--;
               }
            }
            Input.useVector(TmpIntPOS[i * InternalAtoms_ + j], Hash, "AtomPositions",
                            i * InternalAtoms_ + j); // see above note
         }
      }
      delete[] CellVecs;

      for (int i = 1; i < cnt; ++i)
      {
         for (int j = 0; j < InternalAtoms_; ++j)
         {
            AtomSpecies_[i * InternalAtoms_ + j] = AtomSpecies_[j];
         }
      }

      delete[] InternalPOS_;
      InternalPOS_ = TmpIntPOS;

      InternalAtoms_ = TmpIntAtoms;
   }
}

void CBKinematics::InfluenceRegion(double* const InfluenceRegion)
{
   Matrix Eigvals(1, DIM3);
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
   Eigvals = SymEigVal(F_ * RefLattice_ * ((F_ * RefLattice_).Transpose()));
   tmp = sqrt(Eigvals[0][0]);
   for (int i = 0; i < 3; i++)
   {
      if (sqrt(Eigvals[0][i]) < tmp)
      {
         tmp = sqrt(Eigvals[0][i]);
      }
   }

   // Set to inverse eigenvalue
   tmp = 1.0 / tmp;
   for (int i = 0; i < DIM3; i++)
   {
      InfluenceRegion[i] = tmp;
   }
}

void CBKinematics::SetReferenceDOFs()
{
   DOF_.Resize(DOFS(), 0.0);
   for (int i = 0; i < DIM3; ++i)
   {
      DOF_[INDF(i, i)] = 1.0;
   }
   Reset();
}

void CBKinematics::SetReferenceToCurrent()
{
   Matrix CurrentLattice(DIM3, DIM3, 0.0);
   for (int i = 0; i < DIM3; ++i)
   {
      for (int j = 0; j < DIM3; ++j)
      {
         for (int k = 0; k < DIM3; ++k)
         {
            CurrentLattice[i][j] += F_[j][k] * RefLattice_[i][k];
         }
      }
   }
   RefLattice_ = CurrentLattice;

   for (int i = 0; i < InternalAtoms_; ++i)
   {
      for (int j = 0; j < DIM3; ++j)
      {
         InternalPOS_[i][j] = InternalPOS_[i][j] + S_[i][j];
      }
   }

   SetReferenceDOFs();
}

Vector CBKinematics::CurrentLatticeVec(int const& p) const
{
   Vector tmp(DIM3, 0.0);

   for (int i = 0; i < DIM3; ++i)
   {
      for (int j = 0; j < DIM3; ++j)
      {
         tmp[i] += F_[i][j] * RefLattice_[p][j];
      }
   }

   return tmp;
}
