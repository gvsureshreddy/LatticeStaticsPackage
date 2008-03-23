#include "UnitCellIterator.h"

void UnitCellIterator::Initialize(unsigned GridSize,int DoHalfOnly,int SkipZero)
{
   GridSize_ = GridSize;
   CurrentPOS_ = 0;
   
   if ((GridSize_ < 2))
   {
      cerr << "UnitCellIterator GridSize must be >= 2"
           << "\n";
      exit(-1);
   }
   VectorsLen_ = DoHalfOnly ?
      GridSize*GridSize*(GridSize/2+1) - SkipZero
      : GridSize*GridSize*GridSize - SkipZero;
   
   // Be extra sure there are no memory leaks
   if ( NULL != Vectors_)
   {
      delete [] Vectors_[0];
      delete [] Vectors_;
   }
   
   Vectors_ = new double*[VectorsLen_];
   Vectors_[0] = new double[3*VectorsLen_];
   for (int i=1;i<VectorsLen_;++i)
   {
      Vectors_[i] = Vectors_[i-1] + 3;
   }
   
   int even = 1-GridSize_%2;
   for (unsigned k= DoHalfOnly ? 0 : -GridSize_/2-(!DoHalfOnly&&even); k <= GridSize_/2; ++k)
      for (unsigned j=-GridSize_/2;j<=GridSize_/2-even;++j)
         for (unsigned i=-GridSize_/2;i<=GridSize_/2-even;++i)
         {
            if (!(SkipZero && (i == 0) && (j == 0) && (k == 0)))
            {
               Vectors_[CurrentPOS_][0] = i/double(GridSize_);
               Vectors_[CurrentPOS_][1] = j/double(GridSize_);
               Vectors_[CurrentPOS_][2] = k/double(GridSize_);
               ++CurrentPOS_;
            }
         }
   
   CurrentPOS_ = 0;
}

UnitCellIterator::~UnitCellIterator()
{
   delete [] Vectors_[0];
   delete [] Vectors_;
}
