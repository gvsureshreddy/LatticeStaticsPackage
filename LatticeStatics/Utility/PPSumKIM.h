#ifndef RSE__PPSumKIM
#define RSE__PPSumKIM

#include <Matrix.h>
#include <Vector.h>
#include <vector>
#include "CBKinematics.h"

#define PPSUMKIMdatalen 11
#define PPSUMKIMatomstart 0
#define PPSUMKIMdXstart 2
#define PPSUMKIMdxstart 5
#define PPSUMKIMr2start 8
#define PPSUMKIMphi1start 9
#define PPSUMKIMphi2start 10


class PPSumKIM
{
private:
   int Recalc_;
   double* InfluenceDist_;
   int InternalAtoms_;
   CBKinematics* CBK_;

   int CurrentPOS_;

   std::vector<int> numNeigh_;
   std::vector<int>  nListAtom_;
   std::vector<double> nListRVec_;

   void Initialize();

public:
   PPSumKIM()
   {
   }

   PPSumKIM(CBKinematics* const CBK, int const& InternalAtoms,
            double* const InfluDist);
   ~PPSumKIM()
   {
   }

   void operator()(CBKinematics* const CBK, int const& InternalAtoms,
                   double* const InfluDist);

   void Reset();
   void Recalc()
   {
      Recalc_ = 1;
   }

   int Done() const
   {
      return CurrentPOS_ >= InternalAtoms_;
   }
   int CurrentPOS() const
   {
      return CurrentPOS_;
   }

   void operator++()
   {
      ++CurrentPOS_;
   }

   int& numNeigh()
   {
      return numNeigh_[CurrentPOS_];
   }

   int* nListAtom()
   {
      int temp = 0;
      if (CurrentPOS_ != 0)
      {
         for (int i = 1; i <= (CurrentPOS_); i++)
         {
            temp += numNeigh_[i - 1];
         }
      }
      return &(nListAtom_[temp]);
   }

   double* nListRVec()
   {
      int temp = 0;
      if (CurrentPOS_ != 0)
      {
         for (int i = 1; i <= (CurrentPOS_); i++)
         {
            temp += numNeigh_[i - 1];
         }
         temp = temp * 3;
      }
      return &(nListRVec_[temp]);
   }
};

#endif
