#ifndef __UnitCellIterator
#define __UnitCellIterator

#include <Matrix.h>
#include <Vector.h>

class UnitCellIterator
{
private:
   const Matrix *LatticeVec_;
   const double *RefLen_;

   // Cubic grid of size ((2^(GridSize_)) + 1)
   int GridSize_;
   int VectorsLen_;
   double **Vectors_;
   int CurrentPOS_;

public:

   UnitCellIterator() {}
   UnitCellIterator(int GridSize,Matrix *LatticeVec,double *RefLen);
   ~UnitCellIterator();

   void Reset() {CurrentPOS_ = 0;}

   int Done() {return VectorsLen_ <= CurrentPOS_;}
   const double operator[](int i) {return Vectors_[CurrentPOS_][i];}
   void operator++() {++CurrentPOS_;}
};

#endif

