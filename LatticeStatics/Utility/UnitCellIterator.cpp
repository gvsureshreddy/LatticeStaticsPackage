#include "UnitCellIterator.h"
#include <cstdlib>

void UnitCellIterator::Initialize(int const& GridSize,int const& DoHalfOnly,int const& SkipZero)
{
   GridSize_ = GridSize;
   CurrentPOS_ = 0;
   if (SkipZero > 1 || SkipZero < 0)
   {
      cerr << "UnitCellIterator SkipZero must be zero or one." << "\n";
      exit(-1);
   }
   
   if ((GridSize_ < 3) || !(GridSize_%2))
   {
      cerr << "UnitCellIterator GridSize must be > 2 and odd"
           << "\n";
      exit(-1);
   }
   VectorsLen_ = DoHalfOnly ?
      (GridSize_*GridSize_*(GridSize_/2) // Everything above z=0
       + GridSize_*(GridSize_/2)         // z=0 Everything above y=0
       + GridSize_/2+1                   // z=0 y=0 Everything from x=0 to x=GridSize/2+1
       - SkipZero)                       // use (0,0,0) or not
      : (GridSize_*GridSize_*GridSize_ - SkipZero); // Everything

   // Be extra sure there are no memory leaks
   if ( 0 != Vectors_)
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

   double norm = 2.0*(GridSize_/2);
   if (DoHalfOnly)
   {
      for (int i=(SkipZero?1:0);i<=GridSize_/2;++i)
      {
         Vectors_[CurrentPOS_][0] = i/norm;
         Vectors_[CurrentPOS_][1] = Vectors_[CurrentPOS_][2] = 0.0;
         ++CurrentPOS_;
      }

      for (int j=1;j<=GridSize_/2;++j)
         for (int i=-GridSize_/2;i<=GridSize_/2;++i)
         {
            Vectors_[CurrentPOS_][0] = i/norm;
            Vectors_[CurrentPOS_][1] = j/norm;
            Vectors_[CurrentPOS_][2] = 0.0;
            ++CurrentPOS_;
         }

      for (int k=1;k<=GridSize_/2;++k)
         for (int j=-GridSize_/2;j<=GridSize_/2;++j)
            for (int i=-GridSize_/2;i<=GridSize_/2;++i)
            {
               Vectors_[CurrentPOS_][0] = i/norm;
               Vectors_[CurrentPOS_][1] = j/norm;
               Vectors_[CurrentPOS_][2] = k/norm;
               ++CurrentPOS_;
            }
   }
   else
   {
      for (int k= -GridSize_/2; k <= GridSize_/2; ++k)
         for (int j=-GridSize_/2;j<=GridSize_/2;++j)
            for (int i=-GridSize_/2;i<=GridSize_/2;++i)
            {
               if (!(SkipZero && (i == 0) && (j == 0) && (k == 0)))
               {
                  Vectors_[CurrentPOS_][0] = i/norm;
                  Vectors_[CurrentPOS_][1] = j/norm;
                  Vectors_[CurrentPOS_][2] = k/norm;
                  ++CurrentPOS_;
               }
            }
   }
   CurrentPOS_ = 0;
}

UnitCellIterator::~UnitCellIterator()
{
   delete [] Vectors_[0];
   delete [] Vectors_;
}
