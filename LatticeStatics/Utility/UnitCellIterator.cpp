#include "UnitCellIterator.h"

void UnitCellIterator::Initialize(int GridSize,int DoHalfOnly,int SkipZero)
{
   GridSize_ = GridSize;
   VectorsLen_ = DoHalfOnly ?
      int(pow(pow(2,GridSize)+1,2)*(pow(2,GridSize-1) + 1) - SkipZero)
      :
      int(pow(pow(2,GridSize)+1,3) - SkipZero);
   CurrentPOS_ = 0;

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
   
   double offset;
   double onehalf = 0.5;
   
   int LowerK = DoHalfOnly ? 1 : 0;
   // Place nodes at corners of cube
   for (int k=LowerK;k<2;++k)
      for (int j=0;j<2;++j)
	 for (int i=0;i<2;++i)
	 {
	    Vectors_[CurrentPOS_][0] = (i-1/2.0);
	    Vectors_[CurrentPOS_][1] = (j-1/2.0);
	    Vectors_[CurrentPOS_][2] = (k-1/2.0);
	    ++CurrentPOS_;
	 }

   for (int n=0;n<GridSize_;++n)
   {
      double twon = pow(2,n),
	 twon1 = pow(2,n+1);
      int twonm1 = int((n!=0) ? pow(2,n-1) : 0);
      int centroidstart = (n!=0) ? ( DoHalfOnly ? twonm1 : 0) : SkipZero;
      
      offset = (twon-1)/twon1;

      // Place nodes at centroid of n level squares.
      for (int k=centroidstart;k<twon;++k)
      {
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = ((i/twon)-offset);
	       Vectors_[CurrentPOS_][1] = ((j/twon)-offset);
	       Vectors_[CurrentPOS_][2] = ((k/twon)-offset);
	       ++CurrentPOS_;
	    }
	 }
      }

      // Place nodes at corners of n+1 level squares.
      LowerK = DoHalfOnly ? twonm1 : 0;
      for (int k=LowerK;k<twon;++k)
      {
	 for (int j=0;j<=twon;++j)
	 {
	    for (int i=0;i<=twon1;++i)
	    {
	       Vectors_[CurrentPOS_][0] = ((i/twon1) - onehalf);
	       Vectors_[CurrentPOS_][1] = ((j/twon) - onehalf);
	       Vectors_[CurrentPOS_][2] = ((k/twon) - offset);
	       ++CurrentPOS_;
	    }
	 }
	 
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<=twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = ((i/twon) - onehalf);
	       Vectors_[CurrentPOS_][1] = ((j/twon) - offset);
	       Vectors_[CurrentPOS_][2] = ((k/twon) - offset);
	       ++CurrentPOS_;
	    }
	 }
      }

      if (twonm1 == 0) ++twonm1;
      LowerK = DoHalfOnly ? twonm1 : 0;
      for (int k=LowerK;k<=twon;++k)
      {
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<=twon1;++i)
	    {
	       Vectors_[CurrentPOS_][0] = ((i/twon1) - onehalf);
	       Vectors_[CurrentPOS_][1] = ((j/twon) - offset);
	       Vectors_[CurrentPOS_][2] = ((k/twon) - onehalf);
	       ++CurrentPOS_;
	    }
	 }

	 for (int j=0;j<=twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = ((i/twon) - offset);
	       Vectors_[CurrentPOS_][1] = ((j/twon) - onehalf);
	       Vectors_[CurrentPOS_][2] = ((k/twon) - onehalf);
	       ++CurrentPOS_;
	    }
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
