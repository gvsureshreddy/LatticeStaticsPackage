#ifndef RSE__SCLDChainSum
#define RSE__SCLDChainSum

#include <Matrix.h>
#include <Vector.h>
#include "KnownPairPotentials.h"

#define SCLDCHAINSUMdatalen 21
#define SCLDCHAINSUMatomstart 0
#define SCLDCHAINSUMdXrefstart 2
#define SCLDCHAINSUMdXstart 3
#define SCLDCHAINSUMdxstart 4
#define SCLDCHAINSUMr2start 5
#define SCLDCHAINSUMphi1start 6
#define SCLDCHAINSUMphi2start 7
#define SCLDCHAINSUMphi3start 8
#define SCLDCHAINSUMphi4start 9
#define SCLDCHAINSUMphi5start 10
#define SCLDCHAINSUMphi6start 11
#define SCLDCHAINSUMphi1Tstart 12
#define SCLDCHAINSUMphi2Tstart 13
#define SCLDCHAINSUMphi3Tstart 14
#define SCLDCHAINSUMphi4Tstart 15
#define SCLDCHAINSUMphi5Tstart 16
#define SCLDCHAINSUMphi1TTstart 17
#define SCLDCHAINSUMphi2TTstart 18
#define SCLDCHAINSUMphi3TTstart 19
#define SCLDCHAINSUMphi4TTstart 20

class SCLDChainSum
{
private:
   int Recalc_;
   double const* InfluanceDist_;
   Vector const* DOF_;
   int LagrangeCB_;
   Matrix const* RefLattice_;
   int InternalAtoms_;
   Vector const* InternalPOS_;
   PairPotentials const* const* const* Potential_;
   double const* Ntemp_;

   int CurrentPOS_;
   int Pairs_;

   double F_;
   int Translations_;
   Vector V_;
   Matrix RelPosDATA_;

   void Initialize();

public:
   SCLDChainSum()
   {
   }
   SCLDChainSum(Vector const* const DOF, int const& LagrangeCB, int const& Translations,
                Matrix const* const RefLat, int const& InternalAtoms,
                Vector const* const InternalPOS,
                PairPotentials const* const* const* const PairPot,
                double const* const InfluDist, double const* const Ntemp);
   ~SCLDChainSum()
   {
   }

   void operator()(Vector const* const DOF, int const& LagrangeCB, int const& Translations,
                   Matrix const* const RefLat, int const& InternalAtoms,
                   Vector const* const InternalPOS,
                   PairPotentials const* const* const* const PairPot,
                   double const* const InfluDist, double const* const Ntemp);

   void Reset();
   void Recalc()
   {
      Recalc_ = 1;
   }
   int Done() const
   {
      return CurrentPOS_ >= Pairs_;
   }
   void operator++()
   {
      ++CurrentPOS_;
   }

   double DXref(int const& i) const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMdXrefstart + i];
   }
   double DX(int const& i) const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMdXstart + i];
   }
   double const* const pDX() const
   {
      return &(RelPosDATA_[CurrentPOS_][SCLDCHAINSUMdXstart]);
   }
   double Dx(int const& i) const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMdxstart + i];
   }
   double const* const pDx() const
   {
      return &(RelPosDATA_[CurrentPOS_][SCLDCHAINSUMdxstart]);
   }
   double r2() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMr2start];
   }
   int Atom(int const& i) const
   {
      return int(RelPosDATA_[CurrentPOS_][SCLDCHAINSUMatomstart + i]);
   }
   double phi1() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi1start];
   }
   double phi2() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi2start];
   }
   double phi3() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi3start];
   }
   double phi4() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi4start];
   }
   double phi5() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi5start];
   }
   double phi6() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi6start];
   }
   double phi1T() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi1Tstart];
   }
   double phi2T() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi2Tstart];
   }
   double phi3T() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi3Tstart];
   }
   double phi4T() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi4Tstart];
   }
   double phi5T() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi5Tstart];
   }
   double phi1TT() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi1TTstart];
   }
   double phi2TT() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi2TTstart];
   }
   double phi3TT() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi3TTstart];
   }
   double phi4TT() const
   {
      return RelPosDATA_[CurrentPOS_][SCLDCHAINSUMphi4TTstart];
   }
   double J() const
   {
      return F_;
   }

   Matrix NeighborDistances(int const& cutoff, double const& eps);

   int Pairs() const
   {
      return Pairs_;
   }
   int Capacity() const
   {
      return RelPosDATA_.Rows();
   }
};

#endif

