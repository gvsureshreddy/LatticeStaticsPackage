#ifndef __PPSum
#define __PPSum

#include <Matrix.h>
#include <Vector.h>

#define DATALEN 9

class PPSum
{
private:
   int Recalc_;
   unsigned *InfluanceDist_;
   Vector *DOF_;
   Matrix *RefLattice_;
   int InternalAtoms_;
   Vector *InternalPOS_;

   unsigned CurrentPOS_;
   unsigned Pairs_;

   Matrix U_;
   Matrix V_;
   Matrix RelPosDATA_;

   void Initialize();
   
public:
   PPSum() {}
   PPSum(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS,
	 unsigned *InfluDist);
   ~PPSum() {}

   void operator()(Vector *DOF,Matrix *RefLat,int InternalAtoms,
		   Vector *InternalPOS,unsigned *InfluDist);

   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int i) {return RelPosDATA_[CurrentPOS_][i];}
   double *pDX() {return RelPosDATA_[CurrentPOS_];}
   double Dx(int i) {return RelPosDATA_[CurrentPOS_][3+i];}
   double *pDx() {return &(RelPosDATA_[CurrentPOS_][3]);}
   double r2() {return RelPosDATA_[CurrentPOS_][6];}
   int Atom(int i) {return int(RelPosDATA_[CurrentPOS_][7+i]);}
   double J() {return U_.Det();}
   Matrix U() {return U_;}
   Matrix UInv() {return U_.Inverse();}
   
   int Pairs() {return Pairs_;}
   int Capacity() {return RelPosDATA_.Rows();}
};

#endif
