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
char SparseMatrixID[] = "$Id: SparseMatrix.cpp,v 1.16 2011/04/16 02:45:42 elliott Exp $";

SparseMatrix::SparseMatrix(Matrix const& A, double const& tol)
{
   // This counts the number of nonzero entries

   int i, j = 0;
   int k = 0;
   int count = 0;

   Rows_ = A.Rows();
   Cols_ = A.Cols();

   for (i = 0; i < Rows_; i++)
   {
      for (j = 0; j < Cols_; j++)
      {
         if (fabs(A[i][j]) > tol)
         {
            ++count;
         }
      }
   }
   NoNonZero_ = count;

   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for (i = 0; i < Rows_; i++)
   {
      for (j = 0; j < Cols_; j++)
      {
         if (fabs(A[i][j]) > tol)
         {
            Row_id_[k] = i;
            Column_id_[k] = j;
            Nonzero_entry_[k] = A[i][j];
            k = k + 1;
         }
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(SparseMatrix const& A)
{
   Rows_ = A.Rows();
   Cols_ = A.Cols();
   NoNonZero_ = A.NoNonZero_;

   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for (register int i = 0; i < NoNonZero_; i++)
   {
      Row_id_[i] = A.Row_id_[i];
      Column_id_[i] = A.Column_id_[i];
      Nonzero_entry_[i] = A.Nonzero_entry_[i];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(int const& Rows, int const& Cols, int const& NoNonZero)
{
   Rows_ = Rows;
   Cols_ = Cols;
   NoNonZero_ = NoNonZero;

   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for (register int i = 0; i < NoNonZero_; i++)
   {
      Row_id_[i] = 0;
      Column_id_[i] = 0;
      Nonzero_entry_[i] = 0;
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(Matrix const& A, int const& NoEntries, double const& tol)
{
   int i, j, k;

   NoNonZero_ = NoEntries;
   Rows_ = A.Rows();
   Cols_ = A.Cols();
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   k = 0;
   for (i = 0; i < Rows_; i++)
   {
      for (j = 0; j < Cols_; j++)
      {
         if (fabs(A[i][j]) > tol)
         {
            Row_id_[k] = i;
            Column_id_[k] = j;
            Nonzero_entry_[k] = A[i][j];
            k = k + 1;
         }
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::~SparseMatrix()
{
   if (!IsNull())
   {
      delete[] Row_id_;
      delete[] Column_id_;
      delete[] Nonzero_entry_;
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void SparseMatrix::Resize(int const& Rows, int const& Cols, int const& NoNonZero)
{
   if ((Rows_ != Rows) || (Cols_ != Cols) || (NoNonZero_ != NoNonZero))
   {
      if (!IsNull())
      {
         delete[] Row_id_;
         delete[] Column_id_;
         delete[] Nonzero_entry_;
      }

      Rows_ = Rows;
      Cols_ = Cols;
      NoNonZero_ = NoNonZero;

      Row_id_ = new int[NoNonZero_];
      Column_id_ = new int[NoNonZero_];
      Nonzero_entry_ = new Elm[NoNonZero_];

      for (register int i = 0; i < NoNonZero_; ++i)
      {
         Row_id_[i] = 0;
         Column_id_[i] = 0;
         Nonzero_entry_[i] = 0;
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void SparseMatrix::GetNonZeroEntry(int const& EntryIndex, int& Row, int& Col, Elm& Entry) const
{
   if ((EntryIndex < 0) || (EntryIndex >= NoNonZero_))
   {
      cerr << "Error in SparseMatrix GetNonZeroEntry(). exiting...\n";
      exit(-32);
   }

   Row = Row_id_[EntryIndex];
   Col = Column_id_[EntryIndex];
   Entry = Nonzero_entry_[EntryIndex];
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void SparseMatrix::SetNonZeroEntry(int const& EntryIndex, int const& Row, int const& Col,
                                   Elm const& Entry)
{
   if ((EntryIndex < 0) || (EntryIndex >= NoNonZero_) || (Row < 0) || (Row >= Rows_) || (Col < 0)
       || (Col >= Cols_))
   {
      cerr << "Error in SparseMatrix SetNonZeroEntry(). exiting...\n";
      exit(-32);
   }

   Row_id_[EntryIndex] = Row;
   Column_id_[EntryIndex] = Col;
   Nonzero_entry_[EntryIndex] = Entry;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Add(Matrix& Y, SparseMatrix const& A, SparseMatrix const& B)
{
   if ((A.Rows_ != B.Rows_) || (A.Cols_ != B.Cols_) || (A.Rows_ != Y.Rows()) || (A.Cols_ != Y.Cols())
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Add() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, A.Cols_, 0.0); // initialize to zero
   for (register int l = 0; l < A.NoNonZero_; l++)
   {
      Y[A.Row_id_[l]][A.Column_id_[l]] += A.Nonzero_entry_[l];
   }
   for (register int m = 0; m < B.NoNonZero_; m++)
   {
      Y[B.Row_id_[m]][B.Column_id_[m]] += B.Nonzero_entry_[m];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Add(Matrix& Y, SparseMatrix const& A, Matrix const& B)
{
   if ((A.Rows_ != B.Rows()) || (A.Cols_ != B.Cols()) || (A.Rows_ != Y.Rows()) || (A.Cols_ != Y.Cols())
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Add() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, A.Cols_, 0.0); // initialize to zero
   for (register int l = 0; l < A.NoNonZero_; l++)
   {
      Y[A.Row_id_[l]][A.Column_id_[l]] += A.Nonzero_entry_[l];
   }
   for (register int m = 0; m < B.Rows(); m++)
   {
      for (register int n = 0; n < B.Cols(); ++n)
      {
         Y[m][n] += B[m][n];
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Add(Matrix& Y, Matrix const& A, SparseMatrix const& B)
{
   if ((A.Rows() != B.Rows_) || (A.Cols() != B.Cols_) || (A.Rows() != Y.Rows())
       || (A.Cols() != Y.Cols()) || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Add() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows(), A.Cols(), 0.0); // initialize to zero
   for (register int m = 0; m < A.Rows(); m++)
   {
      for (register int n = 0; n < A.Cols(); ++n)
      {
         Y[m][n] += A[m][n];
      }
   }
   for (register int l = 0; l < B.NoNonZero_; l++)
   {
      Y[B.Row_id_[l]][B.Column_id_[l]] += B.Nonzero_entry_[l];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix operator-(SparseMatrix const& A)
{
   SparseMatrix B(A);

   for (register int i = 0; i < A.NoNonZero_; i++)
   {
      B.Nonzero_entry_[i] = -A.Nonzero_entry_[i];
   }

   return B;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Subtract(Matrix& Y, SparseMatrix const& A, SparseMatrix const& B)
{
   if ((A.Rows_ != B.Rows_) || (A.Cols_ != B.Cols_) || (A.Rows_ != Y.Rows()) || (A.Cols_ != Y.Cols())
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Subtract() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, A.Cols_, 0.0); // initialize to zero
   for (register int l = 0; l < A.NoNonZero_; l++)
   {
      Y[A.Row_id_[l]][A.Column_id_[l]] += A.Nonzero_entry_[l];
   }
   for (register int m = 0; m < B.NoNonZero_; m++)
   {
      Y[B.Row_id_[m]][B.Column_id_[m]] -= B.Nonzero_entry_[m];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Subtract(Matrix& Y, SparseMatrix const& A, Matrix const& B)
{
   if ((A.Rows_ != B.Rows()) || (A.Cols_ != B.Cols()) || (A.Rows_ != Y.Rows()) || (A.Cols_ != Y.Cols())
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Subtract() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, A.Cols_, 0.0); // initialize to zero
   for (register int l = 0; l < A.NoNonZero_; l++)
   {
      Y[A.Row_id_[l]][A.Column_id_[l]] += A.Nonzero_entry_[l];
   }
   for (register int m = 0; m < B.Rows(); m++)
   {
      for (register int n = 0; n < B.Cols(); ++n)
      {
         Y[m][n] -= B[m][n];
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Subtract(Matrix& Y, Matrix const& A, SparseMatrix const& B)
{
   if ((A.Rows() != B.Rows_) || (A.Cols() != B.Cols_) || (A.Rows() != Y.Rows())
       || (A.Cols() != Y.Cols()) || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Subtract() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows(), A.Cols(), 0.0); // initialize to zero
   for (register int m = 0; m < A.Rows(); m++)
   {
      for (register int n = 0; n < A.Cols(); ++n)
      {
         Y[m][n] += A[m][n];
      }
   }
   for (register int l = 0; l < B.NoNonZero_; l++)
   {
      Y[B.Row_id_[l]][B.Column_id_[l]] -= B.Nonzero_entry_[l];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(SparseMatrix& Y, double const& A, SparseMatrix const& B)
{
   if ((B.Rows_ != Y.Rows_) || (B.Cols_ != Y.Cols_) || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Multiply() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   for (register int i = 0; i < B.NoNonZero_; i++)
   {
      Y.Nonzero_entry_[i] = A * B.Nonzero_entry_[i];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(SparseMatrix& Y, SparseMatrix const& A, double const& B)
{
   if ((A.Rows_ != Y.Rows_) || (A.Cols_ != Y.Cols_) || A.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Multiply() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   for (register int i = 0; i < A.NoNonZero_; i++)
   {
      Y.Nonzero_entry_[i] = B * A.Nonzero_entry_[i];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Matrix& Y, SparseMatrix const& A, SparseMatrix const& B)
{
   if ((Y.Rows() != A.Rows_) || (Y.Cols() != B.Cols_) || (A.Cols_ != B.Rows_) ||
       A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error in SparseMatrix Multiply() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, B.Cols_, 0.0); // Initialize to zero
   for (int i = 0; i < A.NoNonZero(); ++i)
   {
      for (int j = 0; j < B.NoNonZero(); ++j)
      {
         if (A.Column_id_[i] == B.Row_id_[j])
         {
            Y[A.Row_id_[i]][B.Column_id_[j]] += A.Nonzero_entry_[i] * B.Nonzero_entry_[j];
         }
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Matrix& Y, SparseMatrix const& A, Matrix const& B, SparseMatrix const& C)
{
   if ((Y.Rows() != A.Rows_) || (Y.Cols() != C.Cols_) || (A.Cols_ != B.Rows()) || (B.Cols() != C.Rows_)
       || A.IsNull() || B.IsNull() || C.IsNull() || Y.IsNull())
   {
      cerr << "Error In SparseMatrix Multiply : Wrong size Matricies or Null Matrix"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, C.Cols_, 0.0); // Initialize to zero
   for (int i = 0; i < A.NoNonZero_; ++i)
   {
      for (int j = 0; j < C.NoNonZero_; ++j)
      {
         Y[A.Row_id_[i]][C.Column_id_[j]]
            += A.Nonzero_entry_[i] * B[A.Column_id_[i]][C.Row_id_[j]] * C.Nonzero_entry_[j];
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Matrix& Y, SparseMatrix const& A, Matrix const& B)
{
   if ((Y.Rows() != A.Rows_) || (Y.Cols() != B.Cols()) || (A.Cols_ != B.Rows())
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error In SparseMatrix Multiply : Wrong size Matricies or Null Matrix"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, B.Cols(), 0.0); // Initialize to zero
   for (int i = 0; i < A.NoNonZero_; ++i)
   {
      for (int j = 0; j < B.Cols(); ++j)
      {
         Y[A.Row_id_[i]][j] += A.Nonzero_entry_[i] * B[A.Column_id_[i]][j];
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Matrix& Y, Matrix const& A, SparseMatrix const& B)
{
   if ((Y.Rows() != A.Rows()) || (Y.Cols() != B.Cols_) || (A.Cols() != B.Rows_)
       || A.IsNull() || B.IsNull() || Y.IsNull())
   {
      cerr << "Error In SparseMatrix Multiply : Wrong size Matricies or Null Matrix"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows(), B.Cols_, 0.0); // Initialize to zero
   for (int i = 0; i < A.Rows(); ++i)
   {
      for (int j = 0; j < B.NoNonZero_; ++j)
      {
         Y[i][B.Column_id_[j]] += A[i][B.Row_id_[j]] * B.Nonzero_entry_[j];
      }
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Vector& Y, SparseMatrix const& A, Vector const& B)
{
   if ((A.Cols_ != B.Dim()) || (Y.Dim() != A.Rows_) || A.IsNull() || (B.Dim() == 0) || (Y.Dim() == 0))
   {
      cerr << "Error In SparesMatrix Multiply : Wrong size or  Null Matrix or Vector"
           << "\n";
      exit(-1);
   }

   Y.Resize(A.Rows_, 0.0); // Initialize to zero
   for (int i = 0; i < A.NoNonZero_; ++i)
   {
      Y[A.Row_id_[i]] += A.Nonzero_entry_[i] * B[A.Column_id_[i]];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Vector& Y, Vector const& A, SparseMatrix const& B)
{
   if ((A.Dim() != B.Rows_) || (Y.Dim() != B.Cols_) || (A.Dim() == 0) || B.IsNull() || (Y.Dim() == 0))
   {
      cerr << "Error In SparesMatrix Multiply : Wrong size or  Null Matrix or Vector"
           << "\n";
      exit(-1);
   }

   Y.Resize(B.Cols_, 0.0); // Initialize to zero
   for (int i = 0; i < B.NoNonZero_; ++i)
   {
      Y[B.Column_id_[i]] += A[B.Row_id_[i]] * B.Nonzero_entry_[i];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Vector3D& Y, SparseMatrix const& A, Vector3D const& B)
{
   if ((A.Cols_ != V3DLEN) || (A.Rows_ != V3DLEN) || A.IsNull())
   {
      cerr << "Error In SparesMatrix Multiply : Wrong size or  Null Matrix"
           << "\n";
      exit(-1);
   }

   Y[0] = Y[1] = Y[2] = 0.0; // Initialize to zero
   for (int i = 0; i < A.NoNonZero_; ++i)
   {
      Y[A.Row_id_[i]] += A.Nonzero_entry_[i] * B[A.Column_id_[i]];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

void Multiply(Vector3D& Y, Vector3D const& A, SparseMatrix const& B)
{
   if ((V3DLEN != B.Cols_) || (B.Rows_ != V3DLEN) || B.IsNull())
   {
      cerr << "Error In SparesMatrix Multiply : Wrong size or  Null Matrix"
           << "\n";
      exit(-1);
   }

   Y[0] = Y[1] = Y[2] = 0.0;
   for (int i = 0; i < B.NoNonZero_; ++i)
   {
      Y[B.Column_id_[i]] += A[B.Row_id_[i]] * B.Nonzero_entry_[i];
   }
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix& SparseMatrix::operator=(Matrix const& A)
{
   // This counts the number of nonzero entries

   int i, j = 0;
   int k = 0;
   int count = 0;

   if (!IsNull())
   {
      delete[] Row_id_;
      delete[] Column_id_;
      delete[] Nonzero_entry_;
   }

   Rows_ = A.Rows();
   Cols_ = A.Cols();

   for (i = 0; i < Rows_; i++)
   {
      for (j = 0; j < Cols_; j++)
      {
         if (fabs(A[i][j]) < SPARSETOL)
         {
            ++count;
         }
      }
   }
   NoNonZero_ = count;

   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for (i = 0; i < Rows_; i++)
   {
      for (j = 0; j < Cols_; j++)
      {
         if (fabs(A[i][j]) < SPARSETOL)
         {
            Row_id_[k] = i;
            Column_id_[k] = j;
            Nonzero_entry_[k] = A[i][j];
            k = k + 1;
         }
      }
   }
   return *this;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix SparseMatrix::Transpose() const
{
   SparseMatrix B(Cols_, Rows_, NoNonZero_);

   for (register int i = 0; i < NoNonZero_; i++)
   {
      B.Row_id_[i] = Column_id_[i];
      B.Column_id_[i] = Row_id_[i];
      B.Nonzero_entry_[i] = Nonzero_entry_[i];
   }

   return B;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix& SparseMatrix::SetSparseIdentity(int const& Size)
{
   if (!IsNull())
   {
      delete[] Row_id_;
      delete[] Column_id_;
      delete[] Nonzero_entry_;
   }

   NoNonZero_ = Size;
   Rows_ = Size;
   Cols_ = Size;
   Row_id_ = new int[NoNonZero_];
   Column_id_ = new int[NoNonZero_];
   Nonzero_entry_ = new Elm[NoNonZero_];

   for (register int i = 0; i < Size; i++)
   {
      Row_id_[i] = i;
      Column_id_[i] = i;
      Nonzero_entry_[i] = 1;
   }

   return *this;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

Matrix ReverseSparse(SparseMatrix const& A)
{
   Matrix B(A.Rows(), A.Cols(), 0);

   for (register int i = 0; i < A.Rows(); i++)
   {
      for (register int j = 0; j < A.Cols(); j++)
      {
         for (register int k = 0; k < A.NoNonZero(); k++)
         {
            if ((A.Row_id_[k] == i) && (A.Column_id_[k] == j))
            {
               B[i][j] = A.Nonzero_entry_[k];
            }
         }
      }
   }
   return B;
}

// //////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& out, SparseMatrix const& A)
{
   int W = out.width();
   int NoNonZero = A.NoNonZero();

   out << "\n";

   for (register int i = 0; i < NoNonZero; i++)
   {
      out << "Row id = " << setw(W) << A.Row_id_[i]
      << "Column id = " << setw(W) << A.Column_id_[i]
      << "Entry = " << setw(W) << A.Nonzero_entry_[i]
      << "\n";
   }

   out << "\n";

   return out;
}

char const* const SparseMatrix::Revision()
{
   return SparseMatrixID;
}

