#include "UnitCellIterator.h"

UnitCellIterator::UnitCellIterator(int GridSize,Matrix *LatticeVec,double *RefLen)
   : GridSize_(GridSize), LatticeVec_(LatticeVec), RefLen_(RefLen),
     VectorsLen_(int(pow(pow(2,GridSize)+1,2)*(pow(2,GridSize-1) + 1))),
     CurrentPOS_(0)
{
   Vectors_ = new double*[VectorsLen_];
   Vectors_[0] = new double[3*VectorsLen_];
   for (int i=1;i<VectorsLen_;++i)
   {
      Vectors_[i] = Vectors_[i-1] + 3;
   }
   
   double offset;
   double onehalf = 0.5;
   
   // I have changed the k values so as to only
   // cover 1/2 of the cube (positive Z)

   
   // Place nodes at corners of cube
   for (int k=1;k<2;++k)
      for (int j=0;j<2;++j)
	 for (int i=0;i<2;++i)
	 {
	    Vectors_[CurrentPOS_][0] = i-1/2.0;
	    Vectors_[CurrentPOS_][1] = j-1/2.0;
	    Vectors_[CurrentPOS_][2] = k-1/2.0;
	    ++CurrentPOS_;
	 }

   for (int n=0;n<GridSize_;++n)
   {
      double twon = pow(2,n),
	 twon1 = pow(2,n+1);
      int twonm1 = int((n!=0) ? pow(2,n-1) : 0);
      
      offset = (twon-1)/twon1;

      // Place nodes at centroid of n level squares.
      for (int k=twonm1;k<twon;++k)
      {
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = (i/twon)-offset;
	       Vectors_[CurrentPOS_][1] = (j/twon)-offset;
	       Vectors_[CurrentPOS_][2] = (k/twon)-offset;
	       ++CurrentPOS_;
	    }
	 }
      }

      // Place nodes at corners of n+1 level squares.
      for (int k=twonm1;k<twon;++k)
      {
	 for (int j=0;j<=twon;++j)
	 {
	    for (int i=0;i<=twon1;++i)
	    {
	       Vectors_[CurrentPOS_][0] = (i/twon1) - onehalf;
	       Vectors_[CurrentPOS_][1] = (j/twon) - onehalf;
	       Vectors_[CurrentPOS_][2] = (k/twon) - offset;
	       ++CurrentPOS_;
	    }
	 }
	 
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<=twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = (i/twon) - onehalf;
	       Vectors_[CurrentPOS_][1] = (j/twon) - offset;
	       Vectors_[CurrentPOS_][2] = (k/twon) - offset;
	       ++CurrentPOS_;
	    }
	 }
      }

      if (twonm1 == 0) ++twonm1;
      for (int k=twonm1;k<=twon;++k)
      {
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<=twon1;++i)
	    {
	       Vectors_[CurrentPOS_][0] = (i/twon1) - onehalf;
	       Vectors_[CurrentPOS_][1] = (j/twon) - offset;
	       Vectors_[CurrentPOS_][2] = (k/twon) - onehalf;
	       ++CurrentPOS_;
	    }
	 }

	 for (int j=0;j<=twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       Vectors_[CurrentPOS_][0] = (i/twon) - offset;
	       Vectors_[CurrentPOS_][1] = (j/twon) - onehalf;
	       Vectors_[CurrentPOS_][2] = (k/twon) - onehalf;
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
