#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>

// Global IDString
char MatrixID[] = "$Id: Matrix.cpp,v 1.33 2011/05/30 15:48:37 elliott Exp $";

// Private Methods...


// Computes sqrt(a^2 + b^2) without destructive underflow or overflow.
Matrix::Elm pythag(Matrix::Elm const& a, Matrix::Elm const& b);
Matrix::Elm pythag(Matrix::Elm const& a, Matrix::Elm const& b)
{
   Matrix::Elm absa, absb;
   absa = fabs(a);
   absb = fabs(b);
   if (absa > absb)
   {
      return absa * sqrt(1.0 + ((absb / absa) * (absb / absa)));
   }
   else
   {
      return (absb == 0.0 ?
              0.0 :
              absb * sqrt(1.0 + ((absa / absb) * (absa / absb))));
   }
}


// Returns matrix of size Rows_-1 x Cols_-1 with ith row and
//   jth column removed
Matrix Matrix::Minor(int const& i, int const& j) const
{
   Matrix A(Rows_ - 1, Cols_ - 1);

   if (!IsNull() || !A.IsNull())
   {
      for (register int a = 0; a < A.Rows_; a++)
      {
         for (register int b = 0; b < A.Cols_; b++)
         {
            if (a < i)
            {
               if (b < j)
               {
                  A[a][b] = Elements_[a][b];
               }
               else
               {
                  A[a][b] = Elements_[a][b + 1];
               }
            }
            else
            {
               if (b < j)
               {
                  A[a][b] = Elements_[a + 1][b];
               }
               else
               {
                  A[a][b] = Elements_[a + 1][b + 1];
               }
            }
         }
      }
   }

   return A;
}

// Public Methods...

int Matrix::MathematicaPrintFlag = 0;

Matrix::Matrix(int const& Rows, int const& Cols)
{
   Rows_ = Rows;
   Cols_ = Cols;

   if (IsNull())
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new Matrix::Elm*[Rows_];
      Elements_[0] = new Matrix::Elm[Rows_ * Cols_];

      for (register int i = 1; i < Rows_; i++)
      {
         Elements_[i] = Elements_[i - 1] + Cols_;
      }
   }
   return;
}

Matrix::Matrix(int const& Rows, int const& Cols, Matrix::Elm const& InitVal)
{
   Rows_ = Rows;
   Cols_ = Cols;

   if (IsNull())
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new Matrix::Elm*[Rows_];
      Elements_[0] = new Matrix::Elm[Rows_ * Cols_];

      for (register int i = 1; i < Rows_; i++)
      {
         Elements_[i] = Elements_[i - 1] + Cols_;
      }

      for (register int i = 0; i < Rows_; i++)
      {
         for (register int j = 0; j < Cols_; j++)
         {
            Elements_[i][j] = InitVal;
         }
      }
   }
   return;
}

Matrix::Matrix(Matrix const& A)
{
   Rows_ = A.Rows_;
   Cols_ = A.Cols_;

   if (IsNull())
   {
      Elements_ = 0;
   }
   else
   {
      Elements_ = new Matrix::Elm*[Rows_];
      Elements_[0] = new Matrix::Elm[Rows_ * Cols_];

      for (register int i = 1; i < Rows_; i++)
      {
         Elements_[i] = Elements_[i - 1] + Cols_;
      }

      memmove(Elements_[0], A[0], sizeof(Matrix::Elm) * Rows_ * Cols_);
   }

   return;
}

Matrix::~Matrix()
{
   if (!IsNull())
   {
      delete[] Elements_[0];
      delete[] Elements_;
   }

   return;
}

Matrix operator+(Matrix const& A, Matrix const& B)
{
   if ((A.Rows_ != B.Rows_) || (A.Cols_ != B.Cols_) || A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator+() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Matrix C(A.Rows_, A.Cols_);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         C[i][j] = A[i][j] + B[i][j];
      }
   }

   return C;
}

Matrix operator-(Matrix const& A)
{
   Matrix B(A.Rows_, A.Cols_);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (int j = 0; j < A.Cols_; j++)
      {
         B[i][j] = -A[i][j];
      }
   }

   return B;
}

Matrix operator-(Matrix const& A, Matrix const& B)
{
   if ((A.Rows_ != B.Rows_) || (A.Cols_ != B.Cols_) || A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator-() Diff Size Matrices or Null Matrix!!!"
           << "\n";
      exit(-1);
   }

   Matrix C(A.Rows_, A.Cols_);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         C[i][j] = A[i][j] - B[i][j];
      }
   }

   return C;
}

Matrix operator*(Matrix const& A, Matrix const& B)
{
   if ((A.Cols_ != B.Rows_) || A.IsNull() || B.IsNull())
   {
      cerr << "Error In Matrix Operator* : A.Cols!=B.Rows or Null Matrix"
           << "\n";
      exit(-1);
   }

   Matrix C(A.Rows_, B.Cols_, 0);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < B.Cols_; j++)
      {
         for (register int k = 0; k < A.Cols_; k++)
         {
            C[i][j] += A[i][k] * B[k][j];
         }
      }
   }

   return C;
}

Matrix operator*(Matrix::Elm const& A, Matrix const& B)
{
   Matrix C(B.Rows_, B.Cols_);

   for (register int i = 0; i < B.Rows_; i++)
   {
      for (register int j = 0; j < B.Cols_; j++)
      {
         C[i][j] = A * B[i][j];
      }
   }

   return C;
}

Matrix operator*(Matrix const& A, Matrix::Elm const& B)
{
   Matrix C(A.Rows_, A.Cols_);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         C[i][j] = B * A[i][j];
      }
   }

   return C;
}

Matrix operator/(Matrix const& A, Matrix::Elm const& B)
{
   if (B == 0)
   {
      cerr << "Divide By Zero Error in Matrix operator/()"
           << "\n";
      exit(-1);
   }

   Matrix C(A.Rows_, A.Cols_);

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         C[i][j] = A[i][j] / B;
      }
   }

   return C;
}

#ifdef CHECK_BOUNDS
Matrix::Elm* const Matrix::operator[](int const& i)
{
   if (i >= Rows_)
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}

Matrix::Elm const* const Matrix::operator[](int const& i) const
{
   if (i >= Rows_)
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}
#endif

Matrix& Matrix::operator=(Matrix const& B)
{
   if ((Rows_ != B.Rows_) || (Cols_ != B.Cols_))
   {
      cerr << "Error in Matrix& operator=() : Matricies not same size"
           << "\n";
      exit(-1);
   }

   memmove(Elements_[0], B[0], sizeof(Matrix::Elm) * Rows_ * Cols_);

   return *this;
}

Matrix& Matrix::SetIdentity(int const& Size)
{
   if ((Rows_ != Size) || (Cols_ != Size))
   {
      Resize(Size, Size);
   }

   for (register int i = 0; i < Size; i++)
   {
      for (register int j = 0; j < Size; j++)
      {
         if (i == j)
         {
            Elements_[i][i] = 1.0;
         }
         else
         {
            Elements_[i][j] = 0;
         }
      }
   }

   return *this;
}

Matrix Matrix::Transpose() const
{
   Matrix A(Cols_, Rows_);

   for (register int i = 0; i < Rows_; i++)
   {
      for (register int j = 0; j < Cols_; j++)
      {
         A[j][i] = (*this)[i][j];
      }
   }

   return A;
}

Matrix Matrix::Inverse() const
{
   if (!IsSquare() || IsNull())
   {
      cerr << "Error in Matrix::Inverse() : Non-Square or Null Matrix" << "\n";
      exit(-1);
   }

   Matrix B(Rows_, 1, 0), X(Rows_, 1), C(Rows_, Cols_);

   B[0][0] = 1.0;

   for (register int i = 0; i < Cols_; i++)
   {
#ifdef SOLVE_SVD
      X = SolveSVD(*this, B);
#else
      X = SolvePLU(*this, B);
#endif

      for (register int j = 0; j < Rows_; j++)
      {
         C[j][i] = X[j][0];
      }

      B[i][0] = 0;
      if (i != Cols_ - 1)
      {
         B[i + 1][0] = 1.0;
      }
   }

   return C;
}

void Matrix::Resize(int const& Rows, int const& Cols)
{
   if ((Rows != Rows_) || (Cols != Cols_))
   {
      if (Elements_ != 0)
      {
         delete[] Elements_[0];
         delete[] Elements_;
      }

      Rows_ = Rows;
      Cols_ = Cols;

      if (IsNull())
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new Matrix::Elm*[Rows_];
         Elements_[0] = new Matrix::Elm[Rows_ * Cols_];

         for (register int i = 1; i < Rows_; i++)
         {
            Elements_[i] = Elements_[i - 1] + Cols_;
         }
      }
   }

   return;
}

void Matrix::Resize(int const& Rows, int const& Cols, Matrix::Elm const& InitVal)
{
   if ((Rows != Rows_) || (Cols != Cols_))
   {
      if (Elements_ != 0)
      {
         delete[] Elements_[0];
         delete[] Elements_;
      }

      Rows_ = Rows;
      Cols_ = Cols;

      if (IsNull())
      {
         Elements_ = 0;
      }
      else
      {
         Elements_ = new Matrix::Elm*[Rows_];
         Elements_[0] = new Matrix::Elm[Rows_ * Cols_];

         for (register int i = 1; i < Rows_; i++)
         {
            Elements_[i] = Elements_[i - 1] + Cols_;
         }
      }
   }

   for (register int i = 0; i < Rows_; i++)
   {
      for (register int j = 0; j < Cols_; j++)
      {
         Elements_[i][j] = InitVal;
      }
   }

   return;
}

// Recursivly calculate determinent
Matrix::Elm Matrix::Det() const
{
   if (IsNull() || !IsSquare())
   {
      cerr << "Error in Matrix::Det() : Null or Non-Square Matrix" << "\n";
      exit(-1);
   }

   if (Rows_ == 1)
   {
      return Elements_[0][0];
   }
   else
   {
      Matrix::Elm det = 0;

      for (register int i = 0; i < Cols_; i++)
      {
         det += (1 - 2 * (i % 2)) * Elements_[0][i] * (Minor(0, i).Det());
      }

      return det;
   }
}


// Decompose PA=LU using scaled partial pivoting.
void PLU(Matrix const& A, Matrix& P, Matrix& L, Matrix& U)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in PLU -- Non-Square or Null Matrix to decompose..." << "\n";
      exit(-1);
   }

   P.Resize(A.Rows_, A.Cols_, 0);
   L.SetIdentity(A.Rows_);
   U.Resize(A.Rows_, A.Cols_, 0);

   Matrix Temp = A,
   S(A.Rows_, 1, 0);
   int* Ipivot;
   Ipivot = new int[A.Rows_];

   for (register int i = 0; i < A.Rows_; i++)
   {
      Ipivot[i] = i;
   }

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         if (fabs(Temp[i][j]) > S[i][0])
         {
            S[i][0] = fabs(Temp[i][j]);
         }
      }
   }

   for (register int i = 0; i < A.Rows_; i++)
   {
      Matrix::Elm temp1;
      temp1 = fabs(Temp[i][i] / S[i][0]);

      int k = i;
      for (register int j = i; j < A.Rows_; j++)
      {
         if (fabs(Temp[j][i]) > temp1)
         {
            temp1 = fabs(Temp[j][i] / S[j][0]);
            k = j;
         }
      }

      if (k > i)
      {
         Matrix::Elm* Switch;
         Switch = new Matrix::Elm[A.Rows_];
         for (register int j = i; j < A.Rows_; j++)
         {
            Switch[j] = Temp[i][j];
            Temp[i][j] = Temp[k][j];
            Temp[k][j] = Switch[j];
         }

         for (register int j = 0; j < i; j++)
         {
            Switch[j] = L[i][j];
            L[i][j] = L[k][j];
            L[k][j] = Switch[j];
         }

         delete[] Switch;

         temp1 = S[i][0];
         S[i][0] = S[k][0];
         S[k][0] = temp1;

         int tempi1 = Ipivot[i];
         Ipivot[i] = Ipivot[k];
         Ipivot[k] = tempi1;
      }

      for (register int j = i + 1; j < A.Rows_; j++)
      {
         L[j][i] = Temp[j][i] / Temp[i][i];
      }

      for (register int j = i + 1; j < A.Rows_; j++)
      {
         for (register int k = i + 1; k < A.Rows_; k++)
         {
            Temp[j][k] = Temp[j][k] - (L[j][i] * Temp[i][k]);
         }
      }
   }

   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = i; j < A.Rows_; j++)
      {
         U[i][j] = Temp[i][j];
      }
   }

   for (register int i = 0; i < A.Rows_; i++)
   {
      P[i][Ipivot[i]] = 1;
   }

   delete[] Ipivot;

   return;
}

// find QR factorization of A or A.Transpose()
//
// A   = Q*R  -- CalcTranspose = 0
// A^T = Q*R  -- CalcTranspose = 1
void QR(Matrix const& A, Matrix& Q, Matrix& R, int const& CalcTranspose)
{
   int i, j, k, m, n;
   Matrix::Elm c, s, tau, A1, A2;

   if (CalcTranspose)
   {
      m = A.Cols_;
      n = A.Rows_;
   }
   else
   {
      m = A.Rows_;
      n = A.Cols_;
   }

   // initialize R and Q
   R.Resize(m, n);
   for (i = 0; i < m; ++i)
   {
      for (j = 0; j < n; ++j)
      {
         if (CalcTranspose)
         {
            R[i][j] = A[j][i];
         }
         else
         {
            R[i][j] = A[i][j];
         }
      }
   }
   Q.SetIdentity(m);

   // Perform Givens rotations
   for (j = 0; j < m; j++)
   {
      for (i = m - 1; i > j; i--)
      {
         // calculate G such that G^T*R sets R[i][j] = 0
         if (fabs(R[i][j]) == 0.0)
         {
            c = 1.0; s = 0.0;
         }
         else if (fabs(R[i - 1][j]) == 0.0)
         {
            c = 0.0; s = 1.0;
         }
         else
         {
            if (fabs(R[i][j]) > fabs(R[i - 1][j]))
            {
               tau = -R[i - 1][j] / R[i][j];
               s = 1.0 / sqrt(1.0 + tau * tau);
               c = tau * s;
            }
            else
            {
               tau = -R[i][j] / R[i - 1][j];
               c = 1.0 / sqrt(1.0 + tau * tau);
               s = tau * c;
            }
         }

         for (k = j; k < n; ++k)
         {
            // perform G^T*R
            A1 = R[i - 1][k];
            A2 = R[i][k];
            R[i - 1][k] = A1 * c - A2 * s;
            R[i][k] = A1 * s + A2 * c;
         }

         for (k = 0; k < m; ++k)
         {
            // perform Q*G
            A1 = Q[k][i - 1];
            A2 = Q[k][i];
            Q[k][i - 1] = A1 * c - A2 * s;
            Q[k][i] = A1 * s + A2 * c;
         }
      }
   }
   return;
}

// Return the solution x for the linear system A*x = B
// using A=Q*R if A.Rows()==A.Cols()
// using A^{+} = (A^{T}*A)^{-1}*A^{T}, with A=Q*R if A.Rows() > A.Cols()
// using A^{+} = A^{T}*(A*A^{T})^{-1}, with A^{T}=Q*R if A.Rows() < A.Cols()
void SolveQR(Matrix const& Q, Matrix const& R, Matrix& x, Matrix const& B)
{
   if (Q.Cols() != R.Rows())
   {
      cerr << "Error in SolveQR(): incorrect Q and R matrix sizes.\n";
      exit(-1);
   }

   // x = A^{-1}*B
   if ((x.Rows() == B.Rows()) && (Q.Rows() == B.Rows()) &&
       (Q.Cols() == R.Rows()) && (R.Cols() == B.Rows()))
   {
      // A = Q*R
      //
      // R*x = Q^{T}*B

      // calculate y = Q^{T}*B
      Matrix y(B.Rows(), 1, 0.0);
      for (int i = 0; i < y.Rows(); ++i)
      {
         for (int j = 0; j < Q.Rows(); ++j)
         {
            y[i][0] += Q[j][i] * B[j][0];
         }
      }

      // solve R*x = y
      // note: no need to initialize x
      double sum;
      for (int i = R.Rows() - 1; i >= 0; --i)
      {
         sum = 0.0;
         for (int j = R.Cols() - 1; j > i; --j)
         {
            sum += R[i][j] * x[j][0];
         }
         x[i][0] = (y[i][0] - sum) / R[i][i];
      }
   }
   // x = A^{+}*B
   else if ((x.Rows() == R.Cols()) && (B.Rows() == Q.Cols()))
   {
      // A = Q*R
      //
      // A^{+} = (A^{T}*A)^{-1}*A^{T}
      // using QR gives (A^{T}*A) = (R^{T}*R)
      // so (R^{T}*R)*x = A^{T}*B

      // calculate z = A^{T}*B = R^{T}*Q^{T}*B
      Matrix y(Q.Cols(), 1, 0.0);
      for (int i = 0; i < y.Rows(); ++i)
      {
         for (int j = 0; j < Q.Rows(); ++j)
         {
            y[i][0] += Q[j][i] * B[j][0];
         }
      }
      Matrix z(R.Cols(), 1, 0.0);
      for (int i = 0; i < z.Rows(); ++i)
      {
         for (int j = 0; j < R.Rows(); ++j)
         {
            z[i][0] += R[j][i] * y[j][0];
         }
      }

      // solve (R^{T}*R)*x = z
      //
      // solve R^{T}*y = z
      // note: no need to initialize y
      double sum;
      for (int i = 0; i < y.Rows(); ++i)
      {
         sum = 0.0;
         for (int j = 0; j < i; ++j)
         {
            sum += R[j][i] * y[j][0];
         }
         y[i][0] = (z[i][0] - sum) / R[i][i];
      }

      // solve R*x = y
      // note: no need to initialize x
      for (int i = R.Rows() - 1; i >= 0; --i)
      {
         sum = 0.0;
         for (int j = i + 1; j < R.Cols(); ++j)
         {
            sum += R[i][j] * x[j][0];
         }
         x[i][0] = (y[i][0] - sum) / R[i][i];
      }
   }
   else if ((x.Rows() == Q.Rows()) && (B.Rows() == R.Cols()))
   {
      // A^{T} = Q*R
      //
      // A^{+} = A^{T}*(A*A^{T})^{-1}
      // using QR gives (A*A^{T}) = (R^{T}*R)
      // x = A^{+}*B = A^{T}*(R^{T}*R)^{-1}*B

      // slove (R^{T}*R)*y = B
      //
      // solve R^{T}*z = B
      // note: no need to initialize z
      Matrix z(R.Rows(), 1);
      double sum;
      for (int i = 0; i < R.Rows(); ++i)
      {
         sum = 0.0;
         for (int j = 0; j < i; ++j)
         {
            sum += R[j][i] * z[j][0];
         }
         z[i][0] = (B[i][0] - sum) / R[i][i];
      }
      // solve R*y = z
      // note: no need to initialize y
      Matrix y(R.Cols(), 1);
      for (int i = R.Rows() - 1; i >= 0; --i)
      {
         sum = 0.0;
         for (int j = i + 1; j < R.Cols(); ++j)
         {
            sum += R[i][j] * y[j][0];
         }
         y[i][0] = (z[i][0] - sum) / R[i][i];
      }

      // calculate x = A^{T}*y = Q*R*y
      //
      // calculate z = R*y
      z.Resize(R.Rows(), 1, 0.0); // initialize x
      for (int i = 0; i < z.Rows(); ++i)
      {
         for (int j = 0; j < R.Cols(); ++j)
         {
            z[i][0] += R[i][j] * y[j][0];
         }
      }
      // calculate x = Q*z
      x.Resize(Q.Rows(), 1, 0.0); // initialize z
      for (int i = 0; i < x.Rows(); ++i)
      {
         for (int j = 0; j < Q.Cols(); ++j)
         {
            x[i][0] += Q[i][j] * z[j][0];
         }
      }
   }
   else
   {
      cerr << "Error in SolveQR(): x and B not compatible with size of A.\n";
      exit(-1);
   }
}

// Perform Broyden's update on QR factorization of a matrix (secant eq: Ax = y)
// it is expected that norm(x) == 1.0
// Anew = A + (y-Ax)x^T
void BroydenQRUpdate(Matrix& Q, Matrix& R, Matrix const& y, Matrix const& x)
{
   // A = Q*R so Q^{T}*A = R
   //
   // then Q^{T}*Anew = R + (Q^{T}*y - R*x)*x^{T}
   if ((y.Rows() != Q.Rows()) || (R.Cols() != x.Rows()))
   {
      cerr << "Error in BroydenQRUpdate(): y or x not compatible with A\n";
      exit(-1);
   }
   if (Q.Cols() != R.Rows())
   {
      cerr << "Error in BroydenQRUpdate(): Q and R sizes not compatible\n";
   }

   // following alg 16.3.3 of Allgower and Georg (except A=Q*R not A=Q^T*R)
   Matrix u(Q.Cols(), 1, 0.0);
   for (int i = 0; i < u.Rows(); ++i)
   {
      for (int j = 0; j < Q.Rows(); ++j)
      {
         u[i][0] += Q[j][i] * y[j][0];
      }
      for (int j = 0; j < R.Cols(); ++j)
      {
         u[i][0] -= R[i][j] * x[j][0];
      }
   }

   // Perform Givens rotations so as to render u[i] == 0 for i>0
   double r, c, s;
   double A1, A2;
   for (int i = u.Rows() - 2; i >= 0; --i)
   {
      c = u[i][0];
      s = u[i + 1][0];
      if (s != 0.0)
      {
         r = sqrt(c * c + s * s);
         c = c / r;
         s = s / r;

         // Rotate R
         for (int j = 0; j < R.Cols(); ++j)
         {
            A1 = c * R[i][j] + s * R[i + 1][j];
            A2 = -s * R[i][j] + c * R[i + 1][j];
            R[i][j] = A1;
            R[i + 1][j] = A2;
         }
         // Rotate Q
         for (int j = 0; j < Q.Rows(); ++j)
         {
            A1 = c * Q[j][i] + s * Q[j][i + 1];
            A2 = -s * Q[j][i] + c * Q[j][i + 1];
            Q[j][i] = A1;
            Q[j][i + 1] = A2;
         }
         // Rotate u
         A1 = c * u[i][0] + s * u[i + 1][0];
         A2 = -s * u[i][0] + c * u[i + 1][0];
         u[i][0] = A1;
         u[i + 1][0] = A2;
      }
   }

   // Now we have Hessenberg form
   //
   // Next perform update with u[i] and x
   // then do Givens rotations to bring into uppertriangular form.

   // do the update
   for (int i = 0; i < R.Cols(); ++i)
   {
      R[0][i] += u[0][0] * x[i][0];
   }

   // Perform Givens rotations to fixup R
   for (int i = 0; i < R.Rows() - 1; ++i)
   {
      c = R[i][i];
      s = R[i + 1][i];
      if (s != 0.0)
      {
         r = sqrt(c * c + s * s);
         c = c / r;
         s = s / r;

         // rotate R
         for (int j = 0; j < R.Cols(); ++j)
         {
            A1 = c * R[i][j] + s * R[i + 1][j];
            A2 = -s * R[i][j] + c * R[i + 1][j];
            R[i][j] = A1;
            R[i + 1][j] = A2;
         }
         // Rotate Q
         for (int j = 0; j < Q.Rows(); ++j)
         {
            A1 = c * Q[j][i] + s * Q[j][i + 1];
            A2 = -s * Q[j][i] + c * Q[j][i + 1];
            Q[j][i] = A1;
            Q[j][i + 1] = A2;
         }
      }
   }
}

// Singular Value Decomposition of A -- Algorithm from Numerical Recipies
//
// return value - condition number of A
// A - mxn matrix to decompose
// U - mxn "column-orthogonal" matrix
// W - nxn diagonal matrix (singular values)
// V - nxn orthogonal matrix
//
// A = U*W*V.Transpose();
//
Matrix::Elm SVD(Matrix const& A, Matrix& U, Matrix& W, Matrix& V,
                Matrix::Elm const& MaxCond, int const& PrintFlag)
{
   // Initialize U = A
   U.Resize(A.Rows_, A.Cols_);
   U = A;

   // Initialize W and V
   W.Resize(A.Cols_, A.Cols_, 0.0);
   V.Resize(A.Cols_, A.Cols_);

   // allocate temp storage space
   Matrix::Elm* temp;
   temp = new Matrix::Elm[A.Cols_];

   // define local variables
   int flag,
       l = 0,
       nm = 0;
   Matrix::Elm anorm,
               c, f, g, h, s, scale, x, y, z;

   g = scale = anorm = 0.0;

   // Householder reduction to bidiagonal form.
   for (int i = 0; i < A.Cols_; i++)
   {
      l = i + 1;
      temp[i] = scale * g;
      g = s = scale = 0.0;

      if (i < A.Rows_)
      {
         for (int k = i; k < A.Rows_; k++)
         {
            scale += fabs(U[k][i]);
         }
         if (scale)
         {
            for (int k = i; k < A.Rows_; k++)
            {
               U[k][i] /= scale;
               s += U[k][i] * U[k][i];
            }
            f = U[i][i];
            g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
            h = f * g - s;
            U[i][i] = f - g;
            for (int j = l; j < A.Cols_; j++)
            {
               s = 0.0;
               for (int k = i; k < A.Rows_; k++)
               {
                  s += U[k][i] * U[k][j];
               }
               f = s / h;
               for (int k = i; k < A.Rows_; k++)
               {
                  U[k][j] += f * U[k][i];
               }
            }

            for (int k = i; k < A.Rows_; k++)
            {
               U[k][i] *= scale;
            }
         }
      }

      W[i][i] = scale * g;

      g = s = scale = 0.0;
      if ((i < A.Rows_) && (i != A.Cols_ - 1))
      {
         for (int k = l; k < A.Cols_; k++)
         {
            scale += fabs(U[i][k]);
         }
         if (scale)
         {
            for (int k = l; k < A.Cols_; k++)
            {
               U[i][k] /= scale;
               s += U[i][k] * U[i][k];
            }
            f = U[i][l];
            g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
            h = f * g - s;
            U[i][l] = f - g;
            for (int k = l; k < A.Cols_; k++)
            {
               temp[k] = U[i][k] / h;
            }
            for (int j = l; j < A.Rows_; j++)
            {
               s = 0.0;
               for (int k = l; k < A.Cols_; k++)
               {
                  s += U[j][k] * U[i][k];
               }
               for (int k = l; k < A.Cols_; k++)
               {
                  U[j][k] += s * temp[k];
               }
            }
            for (int k = l; k < A.Cols_; k++)
            {
               U[i][k] *= scale;
            }
         }
      }

      anorm = (anorm > (fabs(W[i][i]) + fabs(temp[i])) ?
               anorm :
               (fabs(W[i][i]) + fabs(temp[i])));
   }

   // Accumulation of right-hand transformations.
   for (int i = A.Cols_ - 1; i >= 0; i--)
   {
      if (i < int(A.Cols_) - 1)
      {
         if (g)
         {
            for (int j = l; j < A.Cols_; j++)
            {
               // Double division to avoid possible underflow
               V[j][i] = (U[i][j] / U[i][l]) / g;
            }
            for (int j = l; j < A.Cols_; j++)
            {
               s = 0.0;
               for (int k = l; k < A.Cols_; k++)
               {
                  s += U[i][k] * V[k][j];
               }
               for (int k = l; k < A.Cols_; k++)
               {
                  V[k][j] += s * V[k][i];
               }
            }
         }
         for (int j = l; j < A.Cols_; j++)
         {
            V[i][j] = V[j][i] = 0.0;
         }
      }
      V[i][i] = 1.0;
      g = temp[i];
      l = i;
   }

   // Accumulation of left-hand transformations.
   for (int i = (A.Cols_ < A.Rows_ ? A.Cols_ : A.Rows_) - 1; i >= 0; i--)
   {
      l = i + 1;
      g = W[i][i];
      for (int j = l; j < A.Cols_; j++)
      {
         U[i][j] = 0.0;
      }
      if (g)
      {
         g = 1.0 / g;
         for (int j = l; j < A.Cols_; j++)
         {
            s = 0.0;
            for (int k = l; k < A.Rows_; k++)
            {
               s += U[k][i] * U[k][j];
            }
            f = (s / U[i][i]) * g;
            for (int k = i; k < A.Rows_; k++)
            {
               U[k][j] += f * U[k][i];
            }
         }
         for (int j = i; j < A.Rows_; j++)
         {
            U[j][i] *= g;
         }
      }
      else
      {
         for (int j = i; j < A.Rows_; j++)
         {
            U[j][i] = 0.0;
         }
      }
      ++U[i][i];
   }

   // Diagonalization of the bidiagonal form: Loop over singular values, and
   // -- over allowed iterations.
   for (int k = A.Cols_ - 1; k >= 0; k--)
   {
      for (int its = 0; its < 30; its++)
      {
         flag = 1;
         // Test for splitting
         for (l = k; l >= 0; l--)
         {
            // Note that temp[0] is always zero
            nm = l - 1;
            if (fabs(temp[l]) + anorm == anorm)
            {
               flag = 0;
               break;
            }
            if (fabs(W[nm][nm]) + anorm == anorm)
            {
               break;
            }
         }
         if (flag)
         {
            // Cancellation of temp[l], if l > 0
            c = 0.0;
            s = 1.0;
            for (int i = l; i <= k; i++)
            {
               f = s * temp[i];
               temp[i] = c * temp[i];
               if (fabs(f) + anorm == anorm)
               {
                  break;
               }
               g = W[i][i];
               h = pythag(f, g);
               W[i][i] = h;
               h = 1.0 / h;
               c = g * h;
               s = -f * h;
               for (int j = 0; j < A.Rows_; j++)
               {
                  y = U[j][nm];
                  z = U[j][i];
                  U[j][nm] = y * c + z * s;
                  U[j][i] = z * c - y * s;
               }
            }
         }
         z = W[k][k];

         // Convergence
         if (l == k)
         {
            // Singular value is made nonnegative.
            if (z < 0.0)
            {
               W[k][k] = -z;
               for (int j = 0; j < A.Cols_; j++)
               {
                  V[j][k] = -V[j][k];
               }
            }
            break;
         }
         if (its == 29)
         {
            cerr << "no convergence in 30 SVD iterations" << "\n";
            exit(-1);
         }
         // Shift from bottom 2-by-2 minor
         x = W[l][l];
         nm = k - 1;
         y = W[nm][nm];
         g = temp[nm];
         h = temp[k];
         f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
         g = pythag(f, 1.0);
         f = ((x - z) * (x + z) + h * ((y / (f + (f >= 0.0 ? fabs(g) : -fabs(g)))) - h)) / x;
         // Next QR transformation:
         c = s = 1.0;
         for (int j = l; j <= nm; j++)
         {
            int i;
            i = j + 1;
            g = temp[i];
            y = W[i][i];
            h = s * g;
            g = c * g;
            z = pythag(f, h);
            temp[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g * c - x * s;
            h = y * s;
            y *= c;
            for (int jj = 0; jj < A.Cols_; jj++)
            {
               x = V[jj][j];
               z = V[jj][i];
               V[jj][j] = x * c + z * s;
               V[jj][i] = z * c - x * s;
            }
            z = pythag(f, h);
            // Rotation can be arbitrary if z = 0
            W[j][j] = z;
            if (z)
            {
               z = 1.0 / z;
               c = f * z;
               s = h * z;
            }
            f = c * g + s * y;
            x = c * y - s * g;
            for (int jj = 0; jj < A.Rows_; jj++)
            {
               y = U[jj][j];
               z = U[jj][i];
               U[jj][j] = y * c + z * s;
               U[jj][i] = z * c - y * s;
            }
         }
         temp[l] = 0.0;
         temp[k] = f;
         W[k][k] = x;
      }
   }

   delete[] temp;

   // Condition number stuff...
   // Remember the singular values are >= 0.0
   Matrix::Elm
      ConditionNumber,
      wmax = 0.0,
      wmin;
   for (int j = 0; j < A.Cols_; j++)
   {
      if (W[j][j] > wmax)
      {
         wmax = W[j][j];
      }
   }
   wmin = wmax;
   for (int j = 0; j < A.Cols_; j++)
   {
      if (W[j][j] < wmin)
      {
         wmin = W[j][j];
      }
   }

   ConditionNumber = wmax / wmin;
   if (PrintFlag)
   {
      cerr << "SVD: Condition Number is : " << ConditionNumber << "\n";
   }

   // Fix up any singular values that are "too small"
   for (int j = 0; j < A.Cols_; j++)
   {
      if (W[j][j] < wmax / MaxCond)
      {
         W[j][j] = 0.0;
         cerr << "SVD: Explicitly set Singular Value #" << j
              << " to 0.0  !!!" << "\n";
      }
   }

   return ConditionNumber;
}

Matrix SymEigVal(Matrix A, Matrix* const B, int const& MaxItr, double const& Tol)
{
   int count = 0,
       converged = 0;
   Matrix EigVals(1, A.Cols_);
   double theta, c, s, cc, ss, cs, aij1, aii1, ajj1, aki1, akj1, tmp;
   const double PIby4 = 0.25 * acos(-1.0);

   if (B != 0)
   {
      B->SetIdentity(A.Cols_);
   }

   while ((count < MaxItr) && (!converged))
   {
      for (int i = 0; i < A.Cols_; i++)
      {
         for (int j = i + 1; j < A.Cols_; j++)
         {
            if (fabs(A[i][j]) < Tol)
            {
               A[i][j] = A[j][i] = 0.0;
               continue;
            }

            if (fabs(A[i][i] - A[j][j]) > Tol)
            {
               theta = 0.5 * atan(2.0 * A[i][j] / (A[i][i] - A[j][j]));
            }
            else
            {
               theta = PIby4 * (A[i][j] / fabs(A[i][j]));
            }

            if (fabs(theta) > PIby4)
            {
               theta -= 2.0 * PIby4 * (theta / fabs(theta));
            }

            c = cos(theta);
            s = sin(theta);
            cc = c * c;
            ss = s * s;
            cs = c * s;
            aij1 = A[i][j];
            aii1 = A[i][i];
            ajj1 = A[j][j];


            A[i][i] = aii1 * cc + 2.0 * aij1 * cs + ajj1 * ss;
            A[j][j] = aii1 * ss - 2.0 * aij1 * cs + ajj1 * cc;
            A[i][j] = A[j][i] = 0.0;

            for (int k = 0; k < A.Cols_; k++)
            {
               if (B != 0)
               {
                  tmp = (*B)[k][i] * c + (*B)[k][j] * s;
                  (*B)[k][j] = -(*B)[k][i] * s + (*B)[k][j] * c;
                  (*B)[k][i] = tmp;
               }


               if ((k == i) || (k == j))
               {
                  continue;
               }

               aki1 = A[k][i];
               akj1 = A[k][j];

               A[i][k] = A[k][i] = aki1 * c + akj1 * s;
               A[j][k] = A[k][j] = -aki1 * s + akj1 * c;
            }
         }
      }
      count++;

      converged = 1;
      for (int i = 0; i < A.Cols_; i++)
      {
         for (int j = i + 1; j < A.Cols_; j++)
         {
            if (fabs(A[i][j]) > Tol)
            {
               converged = 0;
            }
         }
      }
   }

   if (!converged)
   {
      cerr << "Error: SymEigVal(): Failed - No convergence!" << "\n";
      exit(-1);
   }

   for (int i = 0; i < A.Cols_; i++)
   {
      EigVals[0][i] = A[i][i];
   }

   return EigVals;
}

void Cholesky(Matrix const& A, Matrix& U, Matrix& D)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in Cholesky() -- Non-Square or Null Matrix" << "\n";
      exit(-1);
   }

   U.SetIdentity(A.Rows_);
   D.Resize(A.Rows_, A.Cols_, 0);


   for (register int i = 0; i < A.Rows_; i++)
   {
      D[i][i] = A[i][i];
      if (i > 0)
      {
         for (register int k = 0; k < i; k++)
         {
            D[i][i] -= D[k][k] * U[k][i] * U[k][i];
         }
      }
      for (register int j = i + 1; j < A.Rows_; j++)
      {
         U[i][j] = A[i][j];
         if (i > 0)
         {
            for (register int k = 0; k < i; k++)
            {
               U[i][j] -= D[k][k] * U[k][i] * U[k][j];
            }
         }
         U[i][j] /= D[i][i];
      }
   }

   return;
}

Matrix SolvePLU(Matrix const& A, Matrix const& B)
{
   if (!A.IsSquare() || A.IsNull() || (A.Cols_ != B.Rows_))
   {
      cerr << "Error in Solve() - Non-Square Matrix, Null Matrix, or system of != dimension"
           << "\n";
      exit(-1);
   }

   Matrix
      P, L, U; // PLU() resizes P,L,& U thus do not waste time initializing them.

   PLU(A, P, L, U);

   Matrix::Elm* Y;
   Y = new Matrix::Elm[B.Rows_];
   Matrix Temp = P * B;

   Y[0] = Temp[0][0];
   for (register int i = 1; i < B.Rows_; i++)
   {
      Y[i] = Temp[i][0];
      for (register int j = 0; j < i; j++)
      {
         Y[i] -= Y[j] * L[i][j];
      }
   }

   Matrix X(B.Rows_, 1, 0);

   X[B.Rows_ - 1][0] = Y[B.Rows_ - 1] / U[A.Rows_ - 1][A.Cols_ - 1];
   for (register int i = A.Rows_ - 2; i >= 0; i--)
   {
      X[i][0] = Y[i];
      for (register int j = A.Rows_ - 1; j > i; j--)
      {
         X[i][0] -= X[j][0] * U[i][j];
      }
      X[i][0] = X[i][0] / U[i][i];
   }

   delete[] Y;

   return X;
}

Matrix SolvePLU(Matrix const& P, Matrix const& L, Matrix const& U, Matrix const& B)
{
   Matrix::Elm* Y;
   Y = new Matrix::Elm[B.Rows_];
   Matrix Temp = P * B;

   Y[0] = Temp[0][0];
   for (register int i = 1; i < B.Rows_; i++)
   {
      Y[i] = Temp[i][0];
      for (register int j = 0; j < i; j++)
      {
         Y[i] -= Y[j] * L[i][j];
      }
   }

   Matrix X(B.Rows_, 1, 0);

   X[B.Rows_ - 1][0] = Y[B.Rows_ - 1] / U[P.Rows_ - 1][P.Cols_ - 1];
   for (register int i = P.Rows_ - 2; i >= 0; i--)
   {
      X[i][0] = Y[i];
      for (register int j = P.Rows_ - 1; j > i; j--)
      {
         X[i][0] -= X[j][0] * U[i][j];
      }
      X[i][0] = X[i][0] / U[i][i];
   }

   delete[] Y;

   return X;
}

Matrix SolveSVD(Matrix const& A, Matrix const& B,
                Matrix::Elm const& MaxCond, int const& PrintFlag)
{
   // SVD resizes U,W,V so don't bother now
   Matrix
      U, W, V;
   Matrix x(B.Rows_, B.Cols_);

   SVD(A, U, W, V, MaxCond, PrintFlag);

   int jj, j, i;
   Matrix::Elm s, * tmp;

   // Allocate temp space
   tmp = new Matrix::Elm[A.Cols_];

   // Claculate U.Transpose()*B
   for (j = 0; j < A.Cols_; j++)
   {
      s = 0.0;
      // Nonzero result only if W[j][j] is nonzero
      if (W[j][j])
      {
         for (i = 0; i < A.Rows_; i++)
         {
            s += U[i][j] * B[i][0];
         }
         // This is the divide by W[j][j];
         s /= W[j][j];
      }
      tmp[j] = s;
   }

   // Matrix multiply by V to get answer
   for (j = 0; j < A.Cols_; j++)
   {
      s = 0.0;
      for (jj = 0; jj < A.Cols_; jj++)
      {
         s += V[j][jj] * tmp[jj];
      }
      x[j][0] = s;
   }

   // release temp space
   delete[] tmp;

   return x;
}

ostream& operator<<(ostream& out, Matrix const& A)
{
   int W = out.width();

   out << "\n";

   if (Matrix::MathematicaPrintFlag)
   {
      out << setw(0) << "{{";
   }
   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         out << setw(W) << A[i][j];
         if ((Matrix::MathematicaPrintFlag) && (j != (A.Cols_ - 1)))
         {
            out << ",";
         }
      }

      if (Matrix::MathematicaPrintFlag)
      {
         if (i != (A.Rows_ - 1))
         {
            out << "},\n {";
         }
         else
         {
            out << "}}";
         }
      }
      else
      {
         out << "\n";
      }
   }

   out << "\n";

   return out;
}

istream& operator>>(istream& in, Matrix& A)
{
   for (register int i = 0; i < A.Rows_; i++)
   {
      for (register int j = 0; j < A.Cols_; j++)
      {
         in >> A[i][j];
      }
   }

   return in;
}

char const* const Matrix::Revision()
{
   return MatrixID;
}
