#include "SCLDChainSum.h"

int SCLDCHAINSUMcomp(void const* const a, void const* const b);
int SCLDCHAINSUMind(double const& i, double const& j);
#include <cstdlib>

using namespace std;

SCLDChainSum::SCLDChainSum(Vector const* const DOF, int const& LagrangeCB, int const& Translations,
                           Matrix const* const RefLat, int const& InternalAtoms,
                           Vector const* const InternalPOS,
                           PairPotentials const* const* const* const PairPot,
                           double const* const InfluDist, double const* const Ntemp) :
   Recalc_(0),
   InfluanceDist_(InfluDist),
   DOF_(DOF),
   LagrangeCB_(LagrangeCB),
   RefLattice_(RefLat),
   InternalAtoms_(InternalAtoms),
   InternalPOS_(InternalPOS),
   Potential_(PairPot),
   Ntemp_(Ntemp),
   CurrentPOS_(0), Pairs_(0),
   Translations_(Translations),
   V_(InternalAtoms),
   RelPosDATA_(int (2 * (*InfluDist) * InternalAtoms * InternalAtoms), SCLDCHAINSUMdatalen)
{
   Initialize();
}

void SCLDChainSum::operator()(Vector const* const DOF, int const& LagrangeCB,
                              int const& Translations, Matrix const* const RefLat,
                              int const& InternalAtoms, Vector const* const InternalPOS,
                              PairPotentials const* const* const* const PairPot,
                              double const* const InfluDist, double const* const Ntemp)
{
   DOF_ = DOF;
   LagrangeCB_ = LagrangeCB;
   Translations_ = Translations;
   RefLattice_ = RefLat;
   InternalAtoms_ = InternalAtoms;
   InternalPOS_ = InternalPOS;
   InfluanceDist_ = InfluDist;
   V_.Resize(InternalAtoms);
   Recalc_ = 0;
   CurrentPOS_ = 0;
   Pairs_ = 0;
   Potential_ = PairPot;
   Ntemp_ = Ntemp;
   RelPosDATA_.Resize(int(2 * (*InfluDist) * InternalAtoms * InternalAtoms), SCLDCHAINSUMdatalen);

   Initialize();
}

void SCLDChainSum::Reset()
{
   if (Recalc_)
   {
      Initialize();
   }
   else
   {
      CurrentPOS_ = 0;
   }
}

void SCLDChainSum::Initialize()
{
   double X;
   double Influancedist, tmp;
   int p, q;
   int Top, Bottom, CurrentInfluanceDist;

   F_ = (*DOF_)[0];

   if (Translations_)
   {
      for (q = 0; q < InternalAtoms_; ++q)
      {
         V_[q] = (*DOF_)[q + 1];
      }
   }
   else
   {
      V_[0] = 0.0;
      for (q = 1; q < InternalAtoms_; ++q)
      {
         V_[q] = (*DOF_)[q];
      }
   }

   // Set to inverse eigenvalue
   tmp = 1.0 / (((*RefLattice_)[0][0]) * F_);
   Influancedist = tmp * (*InfluanceDist_);

   tmp = 1;
   // set influance distance based on cube size
   CurrentInfluanceDist = int(ceil(Influancedist));
   tmp *= 2.0 * CurrentInfluanceDist;

   Top = CurrentInfluanceDist;
   Bottom = -CurrentInfluanceDist;

   // set tmp to the number of pairs in the cell to be scanned
   tmp *= InternalAtoms_ * InternalAtoms_;
   // Vol of sphere of R=0.5 is 0.52
   // make sure there is enough memory to store a sphere (which fits inside
   // the cube) of pairs.
   if (RelPosDATA_.Rows() < tmp)
   {
      RelPosDATA_.Resize(int(tmp), SCLDCHAINSUMdatalen);
      cerr << "Resizing RELPOSDATA matrix in SCLDChainSum object to " << tmp << "\n";
   }

   Pairs_ = 0;
   for (p = 0; p < InternalAtoms_; p++)
   {
      for (q = 0; q < InternalAtoms_; q++)
      {
         for (X = Bottom; X <= Top; X++)
         {
            RelPosDATA_[Pairs_][SCLDCHAINSUMatomstart] = p;
            RelPosDATA_[Pairs_][SCLDCHAINSUMatomstart + 1] = q;

            // reference position
            RelPosDATA_[Pairs_][SCLDCHAINSUMdXrefstart] =
               (X) *(*RefLattice_)[0][0];

            if (LagrangeCB_)
            {
               RelPosDATA_[Pairs_][SCLDCHAINSUMdXstart] =
                  (X + ((InternalPOS_[q][0] + V_[q])
                        - (InternalPOS_[p][0] + V_[p])))
                  * (*RefLattice_)[0][0];
            }
            else
            {
               RelPosDATA_[Pairs_][SCLDCHAINSUMdXstart] =
                  (X + InternalPOS_[q][0] - InternalPOS_[p][0]) * (*RefLattice_)[0][0];
            }

            if (LagrangeCB_)
            {
               RelPosDATA_[Pairs_][SCLDCHAINSUMdxstart]
                  = F_ * RelPosDATA_[Pairs_][SCLDCHAINSUMdXstart];
            }
            else
            {
               RelPosDATA_[Pairs_][SCLDCHAINSUMdxstart]
                  = F_ * RelPosDATA_[Pairs_][SCLDCHAINSUMdXstart] + V_[q] - V_[p];
            }
            RelPosDATA_[Pairs_][SCLDCHAINSUMr2start] =
               RelPosDATA_[Pairs_][SCLDCHAINSUMdxstart] * RelPosDATA_[Pairs_][SCLDCHAINSUMdxstart];

            // Only use Sphere of Influance (current)
            if ((RelPosDATA_[Pairs_][SCLDCHAINSUMr2start] != 0)
                && (RelPosDATA_[Pairs_][SCLDCHAINSUMr2start]
                    <= (*InfluanceDist_) * (*InfluanceDist_)))
            {
               // calculate phi1, phi2, phi3, phi4, phi5, and phi6
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi1start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::DY);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi2start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D2Y);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi3start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D3Y);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi4start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D4Y);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi5start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D5Y);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi6start] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D6Y);

               // calculate phi1T, phi2T, phi3T, phi4T, and phi5T
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi1Tstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::DY, PairPotentials::DT);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi2Tstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D2Y, PairPotentials::DT);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi3Tstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D3Y, PairPotentials::DT);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi4Tstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D4Y, PairPotentials::DT);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi5Tstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D5Y, PairPotentials::DT);

               // calculate phi1TT, phi2TT, phi3TT, and phi4TT
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi1TTstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::DY, PairPotentials::D2T);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi2TTstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D2Y, PairPotentials::D2T);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi3TTstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D3Y, PairPotentials::D2T);
               RelPosDATA_[Pairs_][SCLDCHAINSUMphi4TTstart] = Potential_[p][q]->PairPotential(
                  *Ntemp_, RelPosDATA_[Pairs_][SCLDCHAINSUMr2start], PairPotentials::D4Y, PairPotentials::D2T);

               ++Pairs_;
            }
         }
      }
   }

   Recalc_ = 0;
   CurrentPOS_ = 0;
}

Matrix SCLDChainSum::NeighborDistances(int const& cutoff, double const& eps)
{
   Reset();
   Matrix NeighborInfo(Pairs_, 3);

   for (int i = 0; i < Pairs_; ++i)
   {
      NeighborInfo[i][0] = RelPosDATA_[i][SCLDCHAINSUMr2start];
      NeighborInfo[i][1] = RelPosDATA_[i][SCLDCHAINSUMatomstart];
      NeighborInfo[i][2] = RelPosDATA_[i][SCLDCHAINSUMatomstart + 1];
   }
   qsort(&(NeighborInfo[0][0]), Pairs_, 3 * sizeof(double), &SCLDCHAINSUMcomp);

   Matrix NeighborDist(cutoff, (InternalAtoms_ * (InternalAtoms_ + 1)) / 2 + 1, 0.0);
   int i = 0, j = 0;
   double CurrentDist;
   while (i < cutoff)
   {
      CurrentDist = NeighborInfo[j][0];
      NeighborDist[i][0] = sqrt(CurrentDist);
      while (fabs(CurrentDist - NeighborInfo[j][0]) < eps)
      {
         ++NeighborDist[i][SCLDCHAINSUMind(NeighborInfo[j][1], NeighborInfo[j][2])];
         ++j;
      }
      ++i;
   }

   return NeighborDist;
}

int SCLDCHAINSUMcomp(void const* const a, void const* const b)
{
   double t;
   if (*((double*) a) == *((double*) b))
   {
      return 0;
   }
   else
   {
      t = *((double*) a) - *((double*) b);
      t /= fabs(t);
      return int(t);
   }
}

int SCLDCHAINSUMind(double const& i, double const& j)
{
   int I = int(i + 1), J = int(j + 1);

   // Make sure I<=J
   if (I > J)
   {
      int s;
      s = I;
      I = J;
      J = s;
   }

   return ((J - 1) * J) / 2 + I;
}
