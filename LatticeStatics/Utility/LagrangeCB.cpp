#include "LagrangeCB.h"

LagrangeCB::LagrangeCB(int const& InternalAtoms, Matrix& RefLattice, Vector* const AtomPositions) :
   CBKinematics(InternalAtoms, RefLattice, AtomPositions)
{
   F_.Resize(DIM3, DIM3);
   S_.Resize(InternalAtoms_, DIM3);

   SetReferenceDOFs();
   Reset();
}

LagrangeCB::LagrangeCB(PerlInput const& Input, PerlInput::HashStruct const* const ParentHash) :
   CBKinematics(Input, ParentHash)
{
   F_.Resize(DIM3, DIM3);
   S_.Resize(InternalAtoms_, DIM3);

   SetReferenceDOFs();
   Reset();
}

void LagrangeCB::Reset()
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

Vector LagrangeCB::FractionalPosVec(int const& p) const
{
   Vector pos(DIM3, 0.0);

   for (int i = 0; i < DIM3; ++i)
   {
      pos[i] = InternalPOS_[p][i] + S_[p][i];
   }

   return pos;
}

double LagrangeCB::DX(double const* const X, int const& p, int const& q, int const& i) const
{
   double tmp = 0.0;

   for (int j = 0; j < DIM3; ++j)
   {
      tmp += (X[j] + ((InternalPOS_[q][j] + S_[q][j])
                      - (InternalPOS_[p][j] + S_[p][j])))
             * RefLattice_[j][i];
   }

   return tmp;
}

double LagrangeCB::Dx(double const* const X, int const& p, int const& q, int const& i) const
{
   double tmp = 0.0;

   for (int k = 0; k < DIM3; ++k)
   {
      tmp += F_[i][k] * DX(X, p, q, k);
   }

   return tmp;
}

double LagrangeCB::DyDF(double const* const Dx, double const* const DX, int const& r, int const& s)
const
{
   return 2.0 * Dx[r] * DX[s];
}

double LagrangeCB::D2yDFF(double const* const DX, int const& r, int const& s, int const& t,
                          int const& u) const
{
   return 2.0 * Del(r, t) * DX[s] * DX[u];
}

double LagrangeCB::DyDS(double const* const Dx, int const& p, int const& q, int const& i,
                        int const& j) const
{
   double ret = 0;

   ret = 0;
   if (DELTA(i, p, q))
   {
      for (int r = 0; r < DIM3; r++)
      {
         for (int k = 0; k < DIM3; k++)
         {
            ret += F_[k][r] * RefLattice_[j][r] * Dx[k];
         }
      }
      ret *= 2.0 * DELTA(i, p, q);
   }

   return ret;
}

double LagrangeCB::D2yDSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                          int const& l) const
{
   double tmp = 0;
   if (DELTA(i, p, q) * DELTA(k, p, q))
   {
      for (int t = 0; t < DIM3; t++)
      {
         for (int r = 0; r < DIM3; r++)
         {
            for (int s = 0; s < DIM3; s++)
            {
               tmp += F_[t][r] * RefLattice_[j][r] * F_[t][s] * RefLattice_[l][s];
            }
         }
      }
      tmp *= 2.0 * DELTA(i, p, q) * DELTA(k, p, q);
   }

   return tmp;
}

double LagrangeCB::D2yDFS(double const* const Dx, double const* const DX, int const& p,
                          int const& q, int const& i, int const& j, int const& k, int const& l)
const
{
   double tmp = 0;

   if (DELTA(k, p, q))
   {
      for (int r = 0; r < DIM3; r++)
      {
         tmp += F_[i][r] * RefLattice_[l][r];
      }
      tmp = 2.0 * DELTA(k, p, q) * (Dx[i] * RefLattice_[l][j] + tmp * DX[j]);
   }

   return tmp;
}

double LagrangeCB::D3yDFFS(double const* const DX, int const& p, int const& q, int const& i,
                           int const& j, int const& k, int const& l, int const& m, int const& n)
const
{
   return 2.0 * DELTA(m, p, q) * Del(i, k) * (RefLattice_[n][j] * DX[l] + DX[j] * RefLattice_[n][l]);
}

double LagrangeCB::D3yDSSF(int const& p, int const& q, int const& i, int const& j, int const& k,
                           int const& l, int const& m, int const& n) const
{
   double tmp = 0;

   if (DELTA(i, p, q) * DELTA(k, p, q))
   {
      for (int s = 0; s < DIM3; s++)
      {
         tmp += F_[m][s] * (RefLattice_[j][n] * RefLattice_[l][s]
                            + RefLattice_[j][s] * RefLattice_[l][n]);
      }
      tmp *= 2.0 * DELTA(i, p, q) * DELTA(k, p, q);
   }

   return tmp;
}

double LagrangeCB::D4yDFFSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                            int const& l, int const& m, int const& n, int const& a, int const& b)
const
{
   return (2.0 * DELTA(m, p, q) * DELTA(a, p, q) * Del(i, k) *
           (RefLattice_[n][j] * RefLattice_[b][l]
            + RefLattice_[n][l] * RefLattice_[b][j]));
}
