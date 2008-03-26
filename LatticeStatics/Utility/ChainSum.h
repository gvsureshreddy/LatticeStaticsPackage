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
   double *InfluanceDist_;
   Vector *DOF_;
   int LagrangeCB_;
   Matrix *RefLattice_;
   int InternalAtoms_;
   Vector *InternalPOS_;
   PairPotentials ***Potential_;
   double *Ntemp_;
   
   int CurrentPOS_;
   int Pairs_;
   
   double F_;
   int Translations_;
   Vector V_;
   Matrix RelPosDATA_;
   
   void Initialize();
   
public:
   ChainSum() {}
   ChainSum(Vector *DOF,int LagrangeCB,int Translations,Matrix *RefLat,
            int InternalAtoms,Vector *InternalPOS,PairPotentials ***PairPot,
            double *InfluDist,double *Ntemp);
   ~ChainSum() {}
   
   void operator()(Vector *DOF,int LagrangeCB,int Translations,Matrix *RefLat,
                   int InternalAtoms,Vector *InternalPOS,PairPotentials ***PairPot,
                   double *InfluDist,double *Ntemp);
   
   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int i) {return RelPosDATA_[CurrentPOS_][CHAINSUMdXstart+i];}
   double *pDX() {return &(RelPosDATA_[CurrentPOS_][CHAINSUMdXstart]);}
   double Dx(int i) {return RelPosDATA_[CurrentPOS_][CHAINSUMdxstart+i];}
   double *pDx() {return &(RelPosDATA_[CurrentPOS_][CHAINSUMdxstart]);}
   double r2() {return RelPosDATA_[CurrentPOS_][CHAINSUMr2start];}
   int Atom(int i) {return int(RelPosDATA_[CurrentPOS_][CHAINSUMatomstart+i]);}
   double phi1() {return RelPosDATA_[CurrentPOS_][CHAINSUMphi1start];}
   double phi2() {return RelPosDATA_[CurrentPOS_][CHAINSUMphi2start];}
   double J() {return F_;}
   
   Matrix NeighborDistances(int cutoff,double eps);
   
   int Pairs() {return Pairs_;}
   int Capacity() {return RelPosDATA_.Rows();}
};

#endif
