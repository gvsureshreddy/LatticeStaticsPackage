#include "UnitCellIterator.h"

void UnitCellIterator::Initialize(int GridSize,Matrix *LatticeVec,int SkipZero)
{
   GridSize_ = GridSize;
   LatticeVec_ = LatticeVec;
   VectorsLen_ = int(pow(pow(2,GridSize)+1,2)*(pow(2,GridSize-1) + 1) - SkipZero);
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
   
   // I have changed the k values so as to only
   // cover 1/2 of the cube (positive Z)

   
   // Place nodes at corners of cube
   for (int k=1;k<2;++k)
      for (int j=0;j<2;++j)
	 for (int i=0;i<2;++i)
	 {
	    for (int p=0;p<3;++p)
	    {
	       Vectors_[CurrentPOS_][p] = (i-1/2.0)*(*LatticeVec_)[0][p];
	       Vectors_[CurrentPOS_][p] += (j-1/2.0)*(*LatticeVec_)[1][p];
	       Vectors_[CurrentPOS_][p] += (k-1/2.0)*(*LatticeVec_)[2][p];
	    }
	    ++CurrentPOS_;
	 }

   for (int n=0;n<GridSize_;++n)
   {
      double twon = pow(2,n),
	 twon1 = pow(2,n+1);
      int twonm1 = int((n!=0) ? pow(2,n-1) : 0);
      int centroidstart = (n!=0) ? twonm1 : SkipZero;
      
      offset = (twon-1)/twon1;

      // Place nodes at centroid of n level squares.
      for (int k=centroidstart;k<twon;++k)
      {
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       for (int p=0;p<3;++p)
	       {
		  Vectors_[CurrentPOS_][p] = ((i/twon)-offset)*(*LatticeVec_)[0][p];
		  Vectors_[CurrentPOS_][p] += ((j/twon)-offset)*(*LatticeVec_)[1][p];
		  Vectors_[CurrentPOS_][p] += ((k/twon)-offset)*(*LatticeVec_)[2][p];
	       }
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
	       for (int p=0;p<3;++p)
	       {
		  Vectors_[CurrentPOS_][p] = ((i/twon1) - onehalf)*(*LatticeVec_)[0][p];
		  Vectors_[CurrentPOS_][p] += ((j/twon) - onehalf)*(*LatticeVec_)[1][p];
		  Vectors_[CurrentPOS_][p] += ((k/twon) - offset)*(*LatticeVec_)[2][p];
	       }
	       ++CurrentPOS_;
	    }
	 }
	 
	 for (int j=0;j<twon;++j)
	 {
	    for (int i=0;i<=twon;++i)
	    {
	       for (int p=0;p<3;++p)
	       {
		  Vectors_[CurrentPOS_][p] = ((i/twon) - onehalf)*(*LatticeVec_)[0][p];
		  Vectors_[CurrentPOS_][p] += ((j/twon) - offset)*(*LatticeVec_)[1][p];
		  Vectors_[CurrentPOS_][p] += ((k/twon) - offset)*(*LatticeVec_)[2][p];
	       }
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
	       for (int p=0;p<3;++p)
	       {
		  Vectors_[CurrentPOS_][p] = ((i/twon1) - onehalf)*(*LatticeVec_)[0][p];
		  Vectors_[CurrentPOS_][p] += ((j/twon) - offset)*(*LatticeVec_)[1][p];
		  Vectors_[CurrentPOS_][p] += ((k/twon) - onehalf)*(*LatticeVec_)[2][p];
	       }
	       ++CurrentPOS_;
	    }
	 }

	 for (int j=0;j<=twon;++j)
	 {
	    for (int i=0;i<twon;++i)
	    {
	       for (int p=0;p<3;++p)
	       {
		  Vectors_[CurrentPOS_][p] = ((i/twon) - offset)*(*LatticeVec_)[0][p];
		  Vectors_[CurrentPOS_][p] += ((j/twon) - onehalf)*(*LatticeVec_)[1][p];
		  Vectors_[CurrentPOS_][p] += ((k/twon) - onehalf)*(*LatticeVec_)[2][p];
	       }
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
