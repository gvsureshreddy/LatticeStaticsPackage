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

   Initialize();
}

void PPSumKIM::Reset()
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

void PPSumKIM::Initialize()
{
   double X[3];
   double Influencedist[3], tmp;
   int p, q, i;
   int Top[3], Bottom[3], CurrentInfluenceDist;

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

   // clear the containers to be sure we start with fresh lists
   numNeigh_.clear();
   nListAtom_.clear();
   nListRVec_.clear();

   int numTemp;
   double r2;
   for (p = 0; p < InternalAtoms_; p++)
   {
      numTemp = 0;
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
                        nListRVec_.push_back(CBK_->Dx(X, p, q, component));
                     }

                     nListAtom_.push_back(q);
                     ++numTemp;
                  }
               }
            }
         }
      }
      numNeigh_.push_back(numTemp);
   }

   Recalc_ = 0;
   CurrentPOS_ = 0;
}
