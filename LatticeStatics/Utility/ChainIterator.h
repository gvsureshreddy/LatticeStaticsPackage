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
   
   void Initialize(int const& GridSize,int const& DoHalfOnly,int const& SkipZero);
   
public:
   
   ChainIterator(): Vectors_(0) {}
   ChainIterator(int const& GridSize,int const& DoHalfOnly=1,int const& SkipZero=1)
      :Vectors_(0) {Initialize(GridSize,DoHalfOnly,SkipZero);}
   ~ChainIterator();
   
   void operator()(int const& GridSize,int const& DoHalfOnly=1,int const& SkipZero=1)
   {Initialize(GridSize,DoHalfOnly,SkipZero);}
   
   void Reset() {CurrentPOS_ = 0;}
   
   int Done() const {return VectorsLen_ <= CurrentPOS_;}
   double const& operator[](int const& i) const {return Vectors_[CurrentPOS_][i];}
   void operator++() {++CurrentPOS_;}
};

#endif

