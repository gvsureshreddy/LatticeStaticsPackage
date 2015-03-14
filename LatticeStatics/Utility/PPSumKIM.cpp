#include "PPSumKIM.h"

#include <cstdlib>

using namespace std;

void PPSumKIM::operator()(CBKinematics* const CBK, int const& InternalAtoms,
                          double* const InfluDist)
{
   CBK_ = CBK;
   InternalAtoms_ = InternalAtoms;
   InfluenceDist_ = InfluDist;
   Recalc_ = 0;
   CurrentPOS_ = 0;

   numNeigh_ = NULL;
   nListAtom_ = NULL;
   nListRVec_ = NULL;
   Initialize();
}

void PPSumKIM::Reset()
{
   if (Recalc_)
   {
      delete[] numNeigh_;
      delete[] nListAtom_;
      delete[] nListRVec_;

      Initialize();
   }
   else
   {
      CurrentPOS_ = 0;
   }
}

void PPSumKIM::Initialize()
{
   double X[3];
   double Influencedist[3], tmp;
   int p, q, i;
   int Top[3], Bottom[3], CurrentInfluenceDist;
   int ListMax, ListVecMax;
   double AtomicDensity;
   double SphereVol;

   CBK_->InfluenceRegion(Influencedist);
   for (i = 0; i < 3; i++)
   {
      Influencedist[i] *= (*InfluenceDist_);
   }

   tmp = 1;
   for (p = 0; p < 3; p++)
   {
      // set influence distance based on cube size
      //
      // also setup to be large enough to encompass Eulerian sphere
      CurrentInfluenceDist = int(ceil(Influencedist[p]));
      tmp *= 2.0 * CurrentInfluenceDist;

      Top[p] = CurrentInfluenceDist;
      Bottom[p] = -CurrentInfluenceDist;
   }

   // set tmp to the number of pairs in the sphere to be scanned
   AtomicDensity = InternalAtoms_ / ((CBK_->DeltaVolume()) * (CBK_->RefVolume()));
   SphereVol = (4.0 * 3.15 / 3.0) * (*InfluenceDist_) * (*InfluenceDist_) * (*InfluenceDist_);
   tmp = ceil(1.15 * InternalAtoms_ * AtomicDensity * SphereVol);

   // make sure there is enough memory to store the pairs.
   ListMax = InternalAtoms_ * int(tmp);
   ListVecMax = InternalAtoms_ * 3 * int(tmp);
   numNeigh_ = new int[InternalAtoms_];
   nListAtom_ = new int[ListMax];
   nListRVec_ = new double[ListVecMax];

   int numTemp0, numTemp1, numTemp2;
   double r2;
   numTemp0 = 0;
   numTemp2 = 0;
   for (p = 0; p < InternalAtoms_; p++)
   {
      numTemp1 = 0;
      for (q = 0; q < InternalAtoms_; q++)
      {
         for (X[0] = Bottom[0]; X[0] <= Top[0]; X[0]++)
         {
            for (X[1] = Bottom[1]; X[1] <= Top[1]; X[1]++)
            {
               for (X[2] = Bottom[2]; X[2] <= Top[2]; X[2]++)
               {
                  r2 = 0.0;
                  for (i = 0; i < 3; i++)
                  {
                     r2 += (CBK_->Dx(X, p, q, i)) * (CBK_->Dx(X, p, q, i));
                  }
                  if ((r2 != 0) && (r2 <= (*InfluenceDist_) * (*InfluenceDist_)))
                  {
                     for (int component = 0; component < 3; component++)
                     {
                        nListRVec_[numTemp2] = CBK_->Dx(X, p, q, component);
                        ++numTemp2;
                     }

                     nListAtom_[numTemp0] = q;
                     ++numTemp0;
                     ++numTemp1;
                  }
               }
            }
         }
      }
      numNeigh_[p] = numTemp1;
   }
   if ((numTemp0 > ListMax) ||
       (numTemp1 > ListMax) ||
       (numTemp2 > ListVecMax))
   {
     cerr << "Error: PPSumKIM::Initialize() - Memory overrun and corrupted." << endl;
     exit(-2);
   }

   Recalc_ = 0;
   CurrentPOS_ = 0;
}
