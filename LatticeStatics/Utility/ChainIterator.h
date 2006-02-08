#ifndef __ChainIterator
#define __ChainIterator

#include <Matrix.h>
#include <Vector.h>

class ChainIterator
{
private:
   // linear grid of size GridSize_
   int GridSize_;
   int VectorsLen_;
   double **Vectors_;
   int CurrentPOS_;
   
   void Initialize(int GridSize,int DoHalfOnly,int SkipZero);

public:

   ChainIterator(): Vectors_(NULL) {}
   ChainIterator(int GridSize,int DoHalfOnly=1,int SkipZero=1)
      :Vectors_(NULL) {Initialize(GridSize,DoHalfOnly,SkipZero);}
   ~ChainIterator();

   void operator()(int GridSize,int DoHalfOnly=1,int SkipZero=1)
   {Initialize(GridSize,DoHalfOnly,SkipZero);}

   void Reset() {CurrentPOS_ = 0;}

   int Done() {return VectorsLen_ <= CurrentPOS_;}
   double operator[](int i) {return Vectors_[CurrentPOS_][i];}
   void operator++() {++CurrentPOS_;}
};

#endif

