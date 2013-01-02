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
      Matrix SuperCellRefLattice(DIM3, DIM3);
      Vector* TmpIntPOS = NULL;
      int TmpIntAtoms;
      int TmpAtomSpecies[CBK_MAX_ATOMS];

      Input.getIntMatrix(&(mu[0][0]), DIM3, DIM3, Hash, "SuperCell");
      SuperCellInfo(mu, SuperCellRefLattice, TmpIntAtoms, TmpIntPOS, TmpAtomSpecies);

      // Overwrite with new valuse and add to input file
      Input.useMatrix(SuperCellRefLattice, Hash, "LatticeBasis"); // should change this so it doesn't print 'default value'
      RefLattice_ = SuperCellRefLattice;

      Input.usePosInt(TmpIntAtoms, Hash, "InternalAtoms"); // see above note
      InternalAtoms_ = TmpIntAtoms;
      for (int i = 0; i < InternalAtoms_; ++i)
      {
         Input.useVector(TmpIntPOS[i], Hash, "AtomPositions", i); // see above note
      }
      delete[] InternalPOS_;
      InternalPOS_ = TmpIntPOS;

      Input.useIntVector(TmpAtomSpecies, InternalAtoms_, Hash, "AtomSpecies");
      for (int i = 0; i < InternalAtoms_; ++i)
      {
         AtomSpecies_[i] = TmpAtomSpecies[i];
      }
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

Vector CBKinematics::CBKtoCoords() const
{
	//THIS IS USED TO OBTAIN ABSOLUTE COORDINATES (FOR USE IN MultiLatticeKIM). RIGHT NOW IT ONLY WORKS WITH SymLagrangeWTransCB
	Vector Coords(3);
	Vector coords(3*InternalAtoms_);

	Vector AtomShift(3);
	for(int i = 0;i < InternalAtoms_; i++)
	{
		for(int j = 0;j < DIM3; j++)
		{
			AtomShift[j]=DOF_[INDS(i,j)];
		}
		Coords = F_ * (RefLattice_.Transpose()) * (InternalPOS_[i]+AtomShift);
		for (int j = 0; j< DIM3; j++)
		{
			coords[i*3+j] = Coords[j];
		}
	}
   /*
    cout << "coords = " << endl;
    for (int i = 0;i < InternalAtoms_; i++)
    {
        for (int j = 0; j< DIM3; j++)
		{
                cout << coords[i*3+j] << ", ";
		}
    }
    cout << endl;
    */
	return coords;
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

void CBKinematics::SuperCellInfo(int const SuperCell[3][3], Matrix& SuperCellRefLattice,
                                 int& SuperCellInternalAtoms, Vector*& SuperCellInternalPOS,
                                 int* const SuperCellAtomSpecies) const
{
   // g+[i] = SuperCell[i][j]*g[j]
   int latrange[DIM3][2];
   Matrix Mu(DIM3, DIM3);
   Matrix MuInvT(DIM3, DIM3);
   Matrix TmpRefLat(DIM3, DIM3, 0.0);
   for (int i = 0; i < DIM3; ++i)
   {
      latrange[i][0] = 0;
      latrange[i][1] = 0;
      for (int j = 0; j < DIM3; ++j)
      {
         Mu[i][j] = double(SuperCell[i][j]);

         if (SuperCell[j][i] < 0)
         {
            if (SuperCell[j][i] < latrange[i][0])
            {
               latrange[i][0] = SuperCell[j][i];
            }
         }
         else
         {
            if (SuperCell[j][i] > latrange[i][1])
            {
               latrange[i][1] = SuperCell[j][i];
            }
         }
      }
   }
   MuInvT = (Mu.Inverse()).Transpose();

   // Find lattice vectors in supercell
   int det = int(Mu.Det());
   int count = 0;
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
               CellVecs[count].Resize(DIM3);
               CellVecs[count][0] = double(i);
               CellVecs[count][1] = double(j);
               CellVecs[count][2] = double(k);
               count++;
            }
         }
      }
   }

   SuperCellRefLattice = Mu * RefLattice_;
   SuperCellInternalAtoms = count * InternalAtoms_;
   if (SuperCellInternalAtoms > CBK_MAX_ATOMS)
   {
      cerr << "Error: CBK::SuperCellInfo, SuperCellInternalAtoms > CBK_MAX_ATOMS="
           << CBK_MAX_ATOMS << ", exiting...\n";
      exit(-23);
   }

   SuperCellInternalPOS = new Vector[SuperCellInternalAtoms];
   for (int i = 0; i < count; ++i)
   {
      for (int j = 0; j < InternalAtoms_; ++j)
      {
         SuperCellInternalPOS[i * InternalAtoms_ + j].Resize(DIM3);
         SuperCellInternalPOS[i * InternalAtoms_ + j] = MuInvT * (CellVecs[i] + InternalPOS_[j]);
         for (int k = 0; k < DIM3; ++k)
         {
            // for the case where atom sits on boundary of cell,
            // make sure it is on "lower left"
            if (SuperCellInternalPOS[i * InternalAtoms_ + j][k] >= 1.0)
            {
               SuperCellInternalPOS[i * InternalAtoms_ + j][k]--;
            }
         }
      }
   }
   delete[] CellVecs;

   for (int i = 0; i < count; ++i)
   {
      for (int j = 0; j < InternalAtoms_; ++j)
      {
         SuperCellAtomSpecies[i * InternalAtoms_ + j] = AtomSpecies_[j];
      }
   }
}
