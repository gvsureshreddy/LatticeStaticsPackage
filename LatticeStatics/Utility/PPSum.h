#ifndef __PPSum
#define __PPSum

#include <Matrix.h>
#include <Vector.h>
#include "KnownPairPotentials.h"

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
   unsigned *InfluanceDist_;
   Vector *DOF_;
   Matrix *RefLattice_;
   int InternalAtoms_;
   Vector *InternalPOS_;
   PairPotentials ***Potential_;
   double *Ntemp_;

   unsigned CurrentPOS_;
   unsigned Pairs_;

   Matrix U_;
   Matrix V_;
   Matrix RelPosDATA_;

   void Initialize();
   
public:
   PPSum() {}
   PPSum(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS,
	 PairPotentials ***PairPot,unsigned *InfluDist,double *Ntemp);
   ~PPSum() {}

   void operator()(Vector *DOF,Matrix *RefLat,int InternalAtoms,
		   Vector *InternalPOS,PairPotentials ***PairPot,
		   unsigned *InfluDist,double *Ntemp);

   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int i) {return RelPosDATA_[CurrentPOS_][PPSUMdXstart+i];}
   double *pDX() {return &(RelPosDATA_[CurrentPOS_][PPSUMdXstart]);}
   double Dx(int i) {return RelPosDATA_[CurrentPOS_][PPSUMdxstart+i];}
   double *pDx() {return &(RelPosDATA_[CurrentPOS_][PPSUMdxstart]);}
   double r2() {return RelPosDATA_[CurrentPOS_][PPSUMr2start];}
   int Atom(int i) {return int(RelPosDATA_[CurrentPOS_][PPSUMatomstart+i]);}
   double phi1() {return RelPosDATA_[CurrentPOS_][PPSUMphi1start];}
   double phi2() {return RelPosDATA_[CurrentPOS_][PPSUMphi2start];}
   double J() {return U_.Det();}
   Matrix U() {return U_;}
   Matrix UInv() {return U_.Inverse();}

   Matrix NeighborDistances(int cutoff,double eps);
   
   int Pairs() {return Pairs_;}
   int Capacity() {return RelPosDATA_.Rows();}
};

#endif
