#include "EulerCB.h"

EulerCB::EulerCB(int const& InternalAtoms, Matrix& RefLattice, Vector* const AtomPositions) :
   CBKinematics(InternalAtoms, RefLattice, AtomPositions)
{
   F_.Resize(DIM3, DIM3);
   S_.Resize(InternalAtoms_, DIM3);

   SetReferenceDOFs();
   Reset();
}

EulerCB::EulerCB(PerlInput const& Input, PerlInput::HashStruct const* const ParentHash) :
   CBKinematics(Input, ParentHash)
{
   F_.Resize(DIM3, DIM3);
   S_.Resize(InternalAtoms_, DIM3);

   SetReferenceDOFs();
   Reset();
}

void EulerCB::Reset()
{
   int i, j, q, p;
   for (i = 0; i < DIM3; ++i)
   {
      for (j = 0; j < DIM3; ++j)
      {
         F_[i][j] = DOF_[INDF(i, j)];
      }
   }

   i = Fsize();
   for (q = 0; q < InternalAtoms_; ++q)
   {
      for (p = 0; p < DIM3; p++)
      {
         S_[q][p] = DOF_[i++];
      }
   }
}

Vector EulerCB::FractionalPosVec(int const& p) const
{
   Vector pos(DIM3, 0.0), fracpos(DIM3, 0.0), tmp(DIM3, 0.0);
   Matrix CurrentLattice(DIM3, DIM3, 0.0), InverseLattice(DIM3, DIM3);

   for (int i = 0; i < DIM3; ++i)
   {
      tmp = CurrentLatticeVec(i);
      for (int j = 0; j < DIM3; ++j)
      {
         CurrentLattice[i][j] = tmp[j];
      }

      pos[i] = InternalPOS_[p][i] + S_[p][i];
   }
   InverseLattice = CurrentLattice.Inverse();


   for (int i = 0; i < DIM3; ++i)
   {
      for (int j = 0; j < DIM3; ++j)
      {
         fracpos[i] += pos[j] * InverseLattice[j][i];
      }
   }

   return fracpos;
}

double EulerCB::DX(double const* const X, int const& p, int const& q, int const& i) const
{
   double tmp = 0.0;

   for (int j = 0; j < DIM3; ++j)
   {
      tmp += X[j] * RefLattice_[j][i];
   }

   return tmp;
}

double EulerCB::Dx(double const* const X, int const& p, int const& q, int const& i) const
{
   double tmp = 0.0;

   for (int j = 0; j < DIM3; ++j)
   {
      tmp += F_[i][j] * DX(X, p, q, j);
   }
   tmp += InternalPOS_[q][i] - InternalPOS_[p][i] + S_[q][i] - S_[p][i];

   return tmp;
}

double EulerCB::DyDF(double const* const Dx, double const* const DX, int const& r, int const& s)
const
{
   return 2.0 * Dx[r] * DX[s];
}

double EulerCB::D2yDFF(double const* const DX, int const& r, int const& s, int const& t,
                       int const& u) const
{
   return 2.0 * Del(r, t) * DX[s] * DX[u];
}

double EulerCB::DyDS(double const* const Dx, int const& p, int const& q, int const& i,
                     int const& j) const
{
   return 2.0 * DELTA(i, p, q) * Dx[j];
}

double EulerCB::D2yDSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                       int const& l) const
{
   return 2.0 * DELTA(i, p, q) * DELTA(k, p, q) * Del(j, l);
}

double EulerCB::D2yDFS(double const* const Dx, double const* const DX, int const& p, int const& q,
                       int const& i, int const& j, int const& k, int const& l) const
{
   return 2.0 * DELTA(k, p, q) * Del(i, l) * DX[j];
}

double EulerCB::D3yDFFS(double const* const DX, int const& p, int const& q, int const& i,
                        int const& j, int const& k, int const& l, int const& m, int const& n) const
{
   return 0.0;
}

double EulerCB::D3yDSSF(int const& p, int const& q, int const& i, int const& j, int const& k,
                        int const& l, int const& m, int const& n) const
{
   return 0.0;
}

double EulerCB::D4yDFFSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                         int const& l, int const& m, int const& n, int const& a, int const& b) const
{
   return 0.0;
}

