#ifndef __PPSum
#define __PPSum

#include <Matrix.h>
#include <Vector.h>
#include "KnownPairPotentials.h"


#define DATALEN 11
#define ATOMSTART 0
#define DXSTART 2
#define DxSTART 5
#define R2START 8
#define PHI1START 9
#define PHI2START 10

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
   
   double DX(int i) {return RelPosDATA_[CurrentPOS_][DXSTART+i];}
   double *pDX() {return &(RelPosDATA_[CurrentPOS_][DXSTART]);}
   double Dx(int i) {return RelPosDATA_[CurrentPOS_][DxSTART+i];}
   double *pDx() {return &(RelPosDATA_[CurrentPOS_][DxSTART]);}
   double r2() {return RelPosDATA_[CurrentPOS_][R2START];}
   int Atom(int i) {return int(RelPosDATA_[CurrentPOS_][ATOMSTART+i]);}
   double phi1() {return RelPosDATA_[CurrentPOS_][PHI1START];}
   double phi2() {return RelPosDATA_[CurrentPOS_][PHI2START];}
   double J() {return U_.Det();}
   Matrix U() {return U_;}
   Matrix UInv() {return U_.Inverse();}
   
   int Pairs() {return Pairs_;}
   int Capacity() {return RelPosDATA_.Rows();}
};

#endif
