#include "ChainIterator.h"
#include <cstdlib>

void ChainIterator::Initialize(int const& GridSize,int const& DoHalfOnly,int const& SkipZero)
{
   GridSize_ = GridSize;
   CurrentPOS_ = 0;
   
   if ((GridSize_ < 2))
   {
      cerr << "ChainIterator GridSize must be >= 2"
           << "\n";
      exit(-1);
   }
   VectorsLen_ = DoHalfOnly ?
      (GridSize/2+1) - SkipZero : GridSize - SkipZero;
   
   // Be extra sure there are no memory leaks
   if ( 0 != Vectors_)
   {
      delete [] Vectors_[0];
      delete [] Vectors_;
   }
   
   Vectors_ = new double*[VectorsLen_];
   Vectors_[0] = new double[VectorsLen_];
   for (int i=1;i<VectorsLen_;++i)
   {
      Vectors_[i] = Vectors_[i-1] + 1;
   }
   
   int even = 1-GridSize_%2;
   for (int k = DoHalfOnly ? 0 : -GridSize_/2-(!DoHalfOnly&&even); k <= GridSize_/2; ++k)
   {
      if (!(SkipZero && (k == 0)))
      {
         Vectors_[CurrentPOS_][0] = k/double(GridSize_);
         ++CurrentPOS_;
      }
   }
   
   CurrentPOS_ = 0;
}

ChainIterator::~ChainIterator()
{
   delete [] Vectors_[0];
   delete [] Vectors_;
}
