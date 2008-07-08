#ifndef __ChainSum
#define __ChainSum

#include <Matrix.h>
#include <Vector.h>
#include "KnownPairPotentials.h"

#define CHAINSUMdatalen 7
#define CHAINSUMatomstart 0
#define CHAINSUMdXstart 2
#define CHAINSUMdxstart 3
#define CHAINSUMr2start 4
#define CHAINSUMphi1start 5
#define CHAINSUMphi2start 6


class ChainSum
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
   ChainSum() {}
   ChainSum(Vector const* const DOF,int const& LagrangeCB,int const& Translations,
            Matrix const* const RefLat,int const& InternalAtoms,Vector const* const InternalPOS,
            PairPotentials const* const* const* const PairPot,
            double const* const InfluDist,double const* const Ntemp);
   ~ChainSum() {}
   
   void operator()(Vector const* const DOF,int const& LagrangeCB,int const& Translations,
                   Matrix const* const RefLat,int const& InternalAtoms,
                   Vector const* const InternalPOS,
                   PairPotentials const* const* const* const PairPot,
                   double const* const InfluDist,double const* const Ntemp);
   
   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() const {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int const& i) const {return RelPosDATA_[CurrentPOS_][CHAINSUMdXstart+i];}
   double const* const pDX() const {return &(RelPosDATA_[CurrentPOS_][CHAINSUMdXstart]);}
   double Dx(int const& i) const {return RelPosDATA_[CurrentPOS_][CHAINSUMdxstart+i];}
   double const* const pDx() const {return &(RelPosDATA_[CurrentPOS_][CHAINSUMdxstart]);}
   double r2() const {return RelPosDATA_[CurrentPOS_][CHAINSUMr2start];}
   int Atom(int const& i) const {return int(RelPosDATA_[CurrentPOS_][CHAINSUMatomstart+i]);}
   double phi1() const {return RelPosDATA_[CurrentPOS_][CHAINSUMphi1start];}
   double phi2() const {return RelPosDATA_[CurrentPOS_][CHAINSUMphi2start];}
   double J() const {return F_;}
   
   Matrix NeighborDistances(int const& cutoff,double const& eps);
   
   int Pairs() const {return Pairs_;}
   int Capacity() const {return RelPosDATA_.Rows();}
};

#endif
