#ifndef RSE__PPSum
#define RSE__PPSum

#include <Matrix.h>
#include <Vector.h>
#include "KnownPairPotentials.h"
#include "CBKinematics.h"

#define PPSUMdatalen 11
#define PPSUMatomstart 0
#define PPSUMdXstart 2
#define PPSUMdxstart 5
#define PPSUMr2start 8
#define PPSUMphi1start 9
#define PPSUMphi2start 10


class PPSum
{
private:
   int Recalc_;
   double *InfluenceDist_;
   int InternalAtoms_;
   PairPotentials ***Potential_;
   double *Ntemp_;
   CBKinematics *CBK_;
   
   int CurrentPOS_;
   int Pairs_;
   
   Matrix RelPosDATA_;
   
   void Initialize();
   
public:
   PPSum() {}
   PPSum(CBKinematics* const CBK,int const& InternalAtoms,PairPotentials*** const PairPot,
         double* const InfluDist,double* const Ntemp);
   ~PPSum() {}
   
   void operator()(CBKinematics* const CBK,int const& InternalAtoms,
                   PairPotentials*** const PairPot,double* const InfluDist,double* const Ntemp);
   
   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() const {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int const& i) const {return RelPosDATA_[CurrentPOS_][PPSUMdXstart+i];}
   double const* const pDX() const {return &(RelPosDATA_[CurrentPOS_][PPSUMdXstart]);}
   double Dx(int const i) const {return RelPosDATA_[CurrentPOS_][PPSUMdxstart+i];}
   double const* const pDx() const {return &(RelPosDATA_[CurrentPOS_][PPSUMdxstart]);}
   double r2() const {return RelPosDATA_[CurrentPOS_][PPSUMr2start];}
   int Atom(int const& i) const {return int(RelPosDATA_[CurrentPOS_][PPSUMatomstart+i]);}
   double phi1() const {return RelPosDATA_[CurrentPOS_][PPSUMphi1start];}
   double phi2() const {return RelPosDATA_[CurrentPOS_][PPSUMphi2start];}
   
   Matrix NeighborDistances(int const& cutoff,double const& eps);
   
   int Pairs() const {return Pairs_;}
   int Capacity() const {return RelPosDATA_.Rows();}
};

#endif
