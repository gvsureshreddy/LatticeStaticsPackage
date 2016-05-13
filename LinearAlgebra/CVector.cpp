#include "CVector.h"
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector.h"
#include "Matrix.h"

// Global IDString
char CVectorID[] = "CVector";

// Private Functions...

// Public Functions...

CVector::CVector(int const Cols)
{
   Cols_ = Cols;

   if (Cols_ == 0)
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new CVector::Elm[Cols_];
   }

   return;
}

CVector::CVector(int const& Cols, CVector::Elm const& InitVal)
{
   Cols_ = Cols;

   if (Cols_ == 0)
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new CVector::Elm[Cols_];

      for (register int j = 0; j < Cols_; j++)
      {
         Elements_[j] = InitVal;
      }
   }

   return;
}

CVector::CVector(CVector const& A)
{
   Cols_ = A.Cols_;

   if (Cols_ == 0)
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new CVector::Elm[Cols_];
   }

   for (int i = 0; i < Cols_; i++)
   {
      Elements_[i] = A[i];
   }

   return;
}

CVector::CVector(CMatrix const& A)
{
   if (A.IsNull() || ((A.Rows() != 1) && (A.Cols() != 1)))
   {
      cerr << "Error in CVector::CVector(CMatrix& A) -- Null CMatrix or Non-CVector"
           << "\n";
      exit(-1);
   }

   if (A.Rows() == 1)
   {
      Cols_ = A.Cols();

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }

      memmove(Elements_, A[0], sizeof(CVector::Elm) * Cols_);
   }
   else
   {
      Cols_ = A.Rows();

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }

      memmove(Elements_, A[0], sizeof(CVector::Elm) * Cols_);
   }

   return;
}

CVector::CVector(Vector const& A)
{
   Cols_ = A.Dim();

   if (Cols_ == 0)
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new CVector::Elm[Cols_];
   }

   for (int i = 0; i < Cols_; i++)
   {
      Elements_[i] = CVector::Elm(A[i]);
   }

   return;
}

CVector::CVector(Matrix const& A)
{
   if (A.IsNull() || ((A.Rows() != 1) && (A.Cols() != 1)))
   {
      cerr << "Error in CVector::CVector(Matrix& A) -- Null Matrix or Non-CVector"
           << "\n";
      exit(-1);
   }

   if (A.Rows() == 1)
   {
      Cols_ = A.Cols();

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }

      for (int i = 0; i < Cols_; i++)
      {
         Elements_[i] = CVector::Elm(A[0][i]);
      }
   }
   else
   {
      Cols_ = A.Rows();

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }

      for (int i = 0; i < Cols_; i++)
      {
         Elements_[i] = CVector::Elm(A[0][i]);
      }
   }

   return;
}

CVector::~CVector()
{
   delete[] Elements_;

   return;
}

CVector operator+(CVector const& A, CVector const& B)
{
   if ((A.Cols_ != B.Cols_) || (A.Cols_ == 0) || (B.Cols_ == 0))
   {
      cerr << "Error in CVector Operator+() Diff Size CVectors or Null CVector!!!"
           << "\n";
      exit(-1);
   }

   CVector C(A.Cols_);

   for (register int j = 0; j < A.Cols_; j++)
   {
      C[j] = A[j] + B[j];
   }

   return C;
}

CVector operator-(CVector const& A, CVector const& B)
{
   if ((A.Cols_ != B.Cols_) || (A.Cols_ == 0) || (B.Cols_ == 0))
   {
      cerr << "Error in CVector Operator-() Diff Size CVectors or Null CVector!!!"
           << "\n";
      exit(-1);
   }

   CVector C(A.Cols_);

   for (register int j = 0; j < A.Cols_; j++)
   {
      C[j] = A[j] - B[j];
   }

   return C;
}

CVector operator-(CVector const& A)
{
   CVector B(A.Cols_);

   for (register int i = 0; i < A.Cols_; i++)
   {
      B[i] = -A[i];
   }

   return B;
}

// Dot Product
CVector::Elm const operator*(CVector const& A, CVector const& B)
{
   if ((A.Cols_ == 0) || (B.Cols_ == 0) || (A.Cols_ != B.Cols_))
   {
      cerr << "Error in Dot Product -- Null CVector or different Dimensions."
           << "\n";
      exit(-1);
   }

   CVector::Elm sum = 0;

   for (register int i = 0; i < A.Cols_; i++)
   {
      sum += A[i] * (B[i].conj());
   }

   return sum;
}

CVector operator*(CMatrix const& A, CVector const& B)
{
   if ((A.Cols() != B.Cols_) || A.IsNull() || (B.Cols_ == 0))
   {
      cerr << "Error In CVector Operator* : A.Cols!=B.Cols or Null CMatrix or CVector"
           << "\n";
      exit(-1);
   }

   CVector C(A.Rows(), 0);

   for (register int i = 0; i < A.Rows(); i++)
   {
      for (register int k = 0; k < B.Cols_; k++)
      {
         C[i] += A[i][k] * B[k];
      }
   }

   return C;
}

CVector operator*(CVector const& A, CMatrix const& B)
{
   if ((B.Cols() != A.Cols_) || B.IsNull() || (A.Cols_ == 0))
   {
      cerr << "Error In CVector Operator* : A.Cols!=B.Cols or Null CMatrix or CVector"
           << "\n";
      exit(-1);
   }

   CVector C(A.Cols_, 0);

   for (register int i = 0; i < B.Cols(); i++)
   {
      for (register int k = 0; k < A.Cols_; k++)
      {
         C[i] += A[k] * B[k][i];
      }
   }

   return C;
}

CVector operator*(CVector::Elm const& A, CVector const& B)
{
   CVector C(B.Cols_);

   for (register int i = 0; i < B.Cols_; i++)
   {
      C[i] = A * B[i];
   }

   return C;
}

CVector operator*(CVector const& A, CVector::Elm const& B)
{
   CVector C(A.Cols_);

   for (register int i = 0; i < A.Cols_; i++)
   {
      C[i] = B * A[i];
   }

   return C;
}

CVector operator/(CVector const& A, CVector::Elm const& B)
{
   CVector C(A.Cols_);

   for (register int i = 0; i < A.Cols_; i++)
   {
      C[i] = A[i] / B;
   }

   return C;
}

#ifdef CHECK_BOUNDS
CVector::Elm& CVector::operator[](int const& i)
{
   if (i >= Cols_)
   {
      cerr << "Error in CVector::operator[]() -- Index Overflow"
           << "\n";
      exit(-1);
   }

   return Elements_[i];
}

CVector::Elm const& CVector::operator[](int const& i) const
{
   if (i >= Cols_)
   {
      cerr << "Error in CVector::operator[]() -- Index Overflow"
           << "\n";
      exit(-1);
   }

   return Elements_[i];
}
#endif

CVector& CVector::operator=(CVector const& B)
{
   if (Cols_ != B.Cols_)
   {
      cerr << "Error in CVector& operator=() : CVectors not same size"
           << "\n";
      exit(-1);
   }

   memmove(Elements_, B.Elements_, sizeof(CVector::Elm) * Cols_);

   return *this;
}

void CVector::Resize(int const& Cols)
{
   if (Cols != Cols_)
   {
      delete[] Elements_;

      Cols_ = Cols;

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }
   }

   return;
}

void CVector::Resize(int const& Cols, CVector::Elm const& InitVal)
{
   if (Cols != Cols_)
   {
      delete[] Elements_;

      Cols_ = Cols;

      if (Cols_ == 0)
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new CVector::Elm[Cols_];
      }

      for (register int j = 0; j < Cols_; j++)
      {
         Elements_[j] = InitVal;
      }
   }

   return;
}

CVector::Elm CVector::Norm() const
{
   return sqrt(*this * (*this));
}

CVector SolvePLU(CMatrix const& A, CVector const& B)
{
   CMatrix C(B.Cols_, 1);

   for (register int i = 0; i < B.Cols_; i++)
   {
      C[i][0] = B[i];
   }

   return SolvePLU(A, C);
}

ostream& operator<<(ostream& out, CVector const& A)
{
   int W = out.width();

   for (register int i = 0; i < A.Cols_; i++)
   {
      out << setw(W) << A[i];
   }

   return out;
}

istream& operator>>(istream& in, CVector& A)
{
   for (register int i = 0; i < A.Cols_; i++)
   {
      in >> A[i];
   }

   return in;
}

char const* const CVector::Revision()
{
   return CVectorID;
}
