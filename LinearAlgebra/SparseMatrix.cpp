#include "SparseMatrix.h"
#include "Matrix.h"
#include "Vector.h"
#include "Vector3D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>

// Global IDString
char SparseMatrixID[]="$Id: SparseMatrix.cpp,v 1.6 2008/06/02 20:50:14 elliott Exp $";

SparseMatrix::SparseMatrix(const Matrix& A)
{
   //This counts the number of nonzero entries
   
   int i,j=0;
   int k=0;
   int count=0;
   
   Rows_=A.Rows();
   Cols_=A.Cols();
   
   for(i=0; i<Rows_;i++)
   {
      for(j=0;j<Cols_;j++)
      {
	 if(A[i][j]!=0)
	 {
	    ++count;
	 }
      }
   }
   NoNonZero_ = count;
	
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for(i=0; i<Rows_;i++)
   {
      for(j=0;j<Cols_;j++)
      {
	 if(A[i][j]!=0)
	 {
	    Row_id_[k] = i;
	    Column_id_[k] = j;
	    Nonzero_entry_[k] = A[i][j];
	    k=k+1;
	 }
      }
   } 		
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(const SparseMatrix& A)
{
   Rows_=A.Rows();
   Cols_=A.Cols();
   NoNonZero_=A.NoNonZero_;
	
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];
	
   for(register int i=0; i<NoNonZero_; i++)
   {
      Row_id_[i] = A.Row_id_[i];
      Column_id_[i] = A.Column_id_[i];
      Nonzero_entry_[i] = A.Nonzero_entry_[i];
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(int NoNonZero, int Rows, int Cols)
{
   Rows_=Rows;
   Cols_=Cols;
   NoNonZero_= NoNonZero;
	
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for(register int i=0; i<NoNonZero_;i++)
   {
      Row_id_[i] = i;
      Column_id_[i] = i;
      Nonzero_entry_[i] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(const Matrix& A,int NoEntries)
{
   int i,j,k;

   NoNonZero_ = NoEntries;
   Rows_ = A.Rows();
   Cols_ = A.Cols();
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];
	
   k=0;
   for(i=0; i<Rows_;i++)
   {
      for(j=0;j<Cols_;j++)
      {
	 if(A[i][j]!=0)
	 {
	    Row_id_[k] = i;
	    Column_id_[k] = j;
	    Nonzero_entry_[k] = A[i][j];
	    k=k+1;
	 }
      }
   } 
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::~SparseMatrix()
{
   if (!IsNull())
   {
      delete [] Row_id_;
      delete [] Column_id_;
      delete [] Nonzero_entry_;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix operator+(const SparseMatrix& A,const SparseMatrix&B)
{
   if (A.Rows_ != B.Rows_ || A.Cols_ != B.Cols_ ||  A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator+() Diff Size Matrices or Null Matrix!!!"
	   << "\n";
      exit(-1);
   }
	
   Matrix C(A.Rows_, A.Cols_, 0);

   for(register int i=0; i < A.Rows_; i++)
   {
      for(register int j=0; j<A.Cols_; j++)
      {
	 for (register int l=0 ; l<A.NoNonZero_; l++)
	 {
	    if (i == A.Row_id_[l] && j == A.Column_id_[l])
	    {
	       C[i][j] = C[i][j]+A.Nonzero_entry_[l];
	    }
	 }
			
	 for (register int m=0; m<B.NoNonZero_; m++)
	 {
	    if (i == B.Row_id_[m] && j == B.Column_id_[m])
	    {
	       C[i][j] = C[i][j] +B.Nonzero_entry_[m];
	    }
	 }
			
      }
   }		
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix operator-(const SparseMatrix& A)
{
   SparseMatrix B(A);
	
   for(register int i=0;i<A.NoNonZero_;i++)
   {
      B.Nonzero_entry_[i]= -A.Nonzero_entry_[i];
   }
	
   return B;
	
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix operator-(const SparseMatrix& A,const SparseMatrix&B)
{
   if (A.Rows_ != B.Rows_ || A.Cols_ != B.Cols_ ||  A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator+() Diff Size Matrices or Null Matrix!!!"
	   << "\n";
      exit(-1);
   }
	
   Matrix C(A.Rows_, A.Cols_, 0);

   for(register int i=0; i < A.Rows_; i++)
   {
      for(register int j=0; j<A.Cols_; j++)
      {
	 for (register int l=0 ; l<A.NoNonZero_; l++)
	 {
	    if (i == A.Row_id_[l] && j == A.Column_id_[l])
	    {
	       C[i][j] = C[i][j]+A.Nonzero_entry_[l];
	    }
	 }
			
	 for (register int m=0; m<B.NoNonZero_; m++)
	 {
	    if (i == B.Row_id_[m] && j == B.Column_id_[m])
	    {
	       C[i][j] = C[i][j] -B.Nonzero_entry_[m];
	    }
	 }
			
      }
   }		
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix operator*(const double A, const SparseMatrix& B)
{
   SparseMatrix C(B);
	
   for(register int i=0;i<B.NoNonZero_;i++)
   {
      C.Nonzero_entry_[i]= A*B.Nonzero_entry_[i];
   }
	
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix operator*(const SparseMatrix& B, const double A)
{
   SparseMatrix C(B);
	
   for(register int i=0;i<B.NoNonZero_;i++)
   {
      C.Nonzero_entry_[i]= A*B.Nonzero_entry_[i];
   }
	
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix operator*(const SparseMatrix& A, const SparseMatrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In SparseMatrix Operator* : A.Cols!=B.Rows or Null Matrix"
	   <<"\n";
      exit(-1);
   }
	
   double entry;
	
	
   Matrix C(A.Rows_, B.Cols_, 0);
	
   for(register int i=0; i< A.Rows_; i++)
   {
      for(register int j=0; j< B.Cols_; j++)
      {
	 entry = 0;
			
	 for(register int k=0; k<A.NoNonZero_;k++)
	 {
	    if(A.Row_id_[k] == i)
	    {
	       for(register int l=0; l<B.NoNonZero_; l++)
	       {
		  if(B.Column_id_[l] == j)
		  {
		     if(B.Row_id_[l] == A.Column_id_[k])
		     {
			entry = entry + (A.Nonzero_entry_[k]*B.Nonzero_entry_[l]);
		     }
		  }
	       }
	    }
	 }
			
	 C[i][j] = entry;
      }		
   }
	
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix operator*(const SparseMatrix& A, const Matrix& B)
{	
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In SparseMatrix Operator* : A.Cols!=B.Rows or Null Matrix"
	   <<"\n";
      exit(-1);
   }
	
   Matrix C(A.Rows(), B.Cols(), 0);
	
   double entry;
	
   for (register int i=0; i<A.Rows(); i++)
   {
      for(register int j=0; j<B.Cols(); j++)
      {
	 entry = 0;
			
	 for(register int k=0; k<B.Rows();k++)
	 {
	    for(register int l=0; l<A.NoNonZero();l++)
	    {
	       if(A.Row_id_[l] == i && A.Column_id_[l] == k)
	       {
		  entry= entry +(A.Nonzero_entry_[l] * B[k][j]);
	       }
	    }	
			
	 }
			
	 C[i][j]=entry;
			
      }
   }
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix operator*(const Matrix& A, const SparseMatrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In SparseMatrix Operator* : A.Cols!=B.Rows or Null Matrix"
	   <<"\n";
      exit(-1);
   }
	
   Matrix C(A.Rows(), B.Cols(), 0);
	
   double entry;
	
   for (register int i=0; i<A.Rows(); i++)
   {
      for(register int j=0; j<B.Cols(); j++)
      {
	 entry = 0;
			
	 for(register int k=0; k< A.Cols();k++)
	 {
	    for (register int l=0;l<B.NoNonZero();l++)
	    {
	       if(B.Column_id_[l] == j & B.Row_id_[l] == k)
	       {
		  entry = entry + (A[i][k] * B.Nonzero_entry_[l]);
	       }
	    }
			
	 }
	 C[i][j]= entry;
      }
   }
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Vector operator*(const SparseMatrix& A, const Vector& B)
{
	
   if (A.Cols()!=B.Cols_ || A.IsNull() || B.Cols_==0)
   {
      cerr << "Error In Vector Operator* : A.Cols!=B.Cols or Null Matrix or Vector"
	   <<"\n";
      exit(-1);
   }
	
	
   Vector C(A.Rows(),0);
	
   double entry;
	
   for(register int i=0; i<A.Rows(); i++)
   {
      entry = 0;
		
      for (register int j=0; j<A.NoNonZero(); j++)
      {
	 if(A.Row_id_[j]==i)
	 {
	    entry = entry + (A.Nonzero_entry_[j] * B[A.Column_id_[j]]);
	 }
      }
		
      C[i] = entry;
   }
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Vector3D operator*(const SparseMatrix& A,const Vector3D& B)
{
   if (A.Cols() != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(matrix,vec): Wrong Size" << "\n";
      exit(-1);
   }

   Vector3D z(0.0);
	
   double entry;
	
   for (register int i=0;i<V3DLEN;i++)
   {
      entry =0;
		
      for(register int j=0; j<V3DLEN; j++)
      {
	 for(register int k=0; k<A.NoNonZero(); k++)
	 {
	    if(A.Row_id_[k] == i && A.Column_id_[k] == j)
	    {
	       entry = entry + (A.Nonzero_entry_[k] * B[j]);
	    }
	 }
      }
      z[i] = entry;
   }

   return z;

}

////////////////////////////////////////////////////////////////////////////////////////////////

Vector3D operator*(const Vector3D& A,const SparseMatrix& B)
{
   if (B.Rows() != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(matrix,vec): Wrong Size" << "\n";
      exit(-1);
   }
	
   Vector3D z(0);
	
   double entry;
	
   for (register int i=0; i < V3DLEN ;i++)
   {
      entry = 0;
		
      for(register int j=0; j< V3DLEN; j++)
      {
	 for(register int k=0; k<B.NoNonZero();k++)
	 {
	    if(B.Row_id_[k] == j && B.Column_id_[k] == i)
	    {
	       entry = entry + (B.Nonzero_entry_[k] * A[j]);
	    }
	 }
      }
      z[i] = entry;
   }
   return z;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Vector operator*(const Vector& A, const SparseMatrix& B)
{
   if (B.Rows()!=A.Cols_ || B.IsNull() || A.Cols_==0)
   {
      cerr << "Error In Vector Operator* : A.Cols!=B.Rows or Null Matrix or Vector"
	   <<"\n";
      exit(-1);
   }
	
   Vector C(B.Cols(), 0);
	
   double entry;
	
   for(register int i=0; i< B.Cols(); i++)
   {
      entry = 0;
		
      for(register int j=0; j<A.Cols_ ;j++)
      {
	 for(register int k=0; k<B.NoNonZero(); k++)
	 {
	    if(B.Column_id_[k]==i && B.Row_id_[k] == j)
	    {
	       entry = entry + (A[j]*B.Nonzero_entry_[k]);
	    }
	 }
      }
      C[i] = entry;
   }

   return C;

}

////////////////////////////////////////////////////////////////////////////////////////////////
SparseMatrix& SparseMatrix::operator=(const Matrix& A)
{
//This counts the number of nonzero entries
   
   int i,j=0;
   int k=0;
   int count=0;
   
   Rows_=A.Rows();
   Cols_=A.Cols();
   
   for(i=0; i<Rows_;i++)
   {
      for(j=0;j<Cols_;j++)
      {
	 if(A[i][j]!=0)
	 {
	    ++count;
	 }
      }
   }
   NoNonZero_ = count;
	
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for(i=0; i<Rows_;i++)
   {
      for(j=0;j<Cols_;j++)
      {
	 if(A[i][j]!=0)
	 {
	    Row_id_[k] = i;
	    Column_id_[k] = j;
	    Nonzero_entry_[k] = A[i][j];
	    k=k+1;
	 }
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////
SparseMatrix SparseMatrix::Transpose() const
{
   SparseMatrix B(NoNonZero_,Cols_,Rows_);
	
   for(register int i=0; i<NoNonZero_; i++)
   {
      B.Row_id_[i] = Column_id_[i];
      B.Column_id_[i] = Row_id_[i];
      B.Nonzero_entry_[i]=Nonzero_entry_[i];
   }
	
   return B;	
}

////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix& SparseMatrix::SetSparseIdentity(int Size)
{
   if (!IsNull())
   {
      delete [] Row_id_;
      delete [] Column_id_;
      delete [] Nonzero_entry_;
   }
   
   NoNonZero_ = Size;
   Rows_ = Size;
   Cols_ = Size;
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];
   
   for(register int i=0; i<Size;i++)
   {
      Row_id_[i] = i;
      Column_id_[i] = i;
      Nonzero_entry_[i] = 1;
   }
   
   return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Matrix ReverseSparse(const SparseMatrix& A)
{
   Matrix B(A.Rows(),A.Cols(),0);
	
   for (register int i=0; i<A.Rows(); i++)
   {
      for(register int j=0; j<A.Cols(); j++)
      {
	 for(register int k=0; k<A.NoNonZero(); k++)
	 {
	    if(A.Row_id_[k] == i && A.Column_id_[k] == j)
	    {
	       B[i][j] = A.Nonzero_entry_[k];
	    }
	 }
      }
   }
   return B;
} 

////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& out,const SparseMatrix& A)
{
   int W=out.width();
   int NoNonZero = A.NoNonZero();
   
   out << "\n";

   for (register int i=0;i<NoNonZero;i++)
   {
      out << "Row id = " << setw(W) << A.Row_id_[i] 
	  << "Column id = " << setw(W) << A.Column_id_[i] 
	  << "Entry = " << setw(W) << A.Nonzero_entry_[i]
	  << "\n";
   }

   out << "\n";

   return out;
}

char* SparseMatrix::Revision()
{
   return SparseMatrixID;
}
