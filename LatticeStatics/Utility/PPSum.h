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
   Matrix *InternalPOS_;

   unsigned CurrentPOS_;
   unsigned Pairs_;

   Matrix U;
   Matrix V;
   Matrix RelPosDATA_;

   void Initialize();
   
public:
   PPSum() {}
   PPSum(Vector *DOF,Matrix *RefLat,int InternalAtoms,Matrix *InternalPOS,
	 unsigned *InfluDist);
   ~PPSum() {}

   void Reset();
   void Recalc() {Recalc_ = 1;}
   int Done() {return CurrentPOS_ >= Pairs_;}
   void operator++() {++CurrentPOS_;}
   
   double DX(int i) {return RelPosDATA_[CurrentPOS_][i];}
   double Dx(int i) {return RelPosDATA_[CurrentPOS_][3+i];}
   double r() {return RelPosDATA_[CurrentPOS_][6];}
   int operator()(int i) {return int(RelPosDATA_[CurrentPOS_][7+i]);}

   int Pairs() {return Pairs_;}
   int Capacity() {return RelPosDATA_.Rows();}
};

#endif
