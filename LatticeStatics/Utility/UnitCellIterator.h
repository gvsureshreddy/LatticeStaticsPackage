#ifndef __UnitCellIterator
#define __UnitCellIterator

#include <Matrix.h>
#include <Vector.h>

class UnitCellIterator
{
private:
   // Cubic grid of size (GridSize_)^3
   int GridSize_;
   int VectorsLen_;
   double **Vectors_;
   int CurrentPOS_;
   
   void Initialize(int GridSize,int DoHalfOnly,int SkipZero);

public:

   UnitCellIterator(): Vectors_(NULL) {}
   UnitCellIterator(int GridSize,int DoHalfOnly=1,int SkipZero=1)
      :Vectors_(NULL) {Initialize(GridSize,DoHalfOnly,SkipZero);}
   ~UnitCellIterator();

   void operator()(int GridSize,int DoHalfOnly=1,int SkipZero=1)
   {Initialize(GridSize,DoHalfOnly,SkipZero);}

   void Reset() {CurrentPOS_ = 0;}

   int Done() {return VectorsLen_ <= CurrentPOS_;}
   const double operator[](int i) {return Vectors_[CurrentPOS_][i];}
   void operator++() {++CurrentPOS_;}
};

#endif

