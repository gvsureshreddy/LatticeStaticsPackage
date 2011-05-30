#ifndef RSE__Vector
#define RSE__Vector

#ifndef RSE__LinearAlgebraBuildDate
#define RSE__LinearAlgebraBuildDate
char* LinearAlgebraBuildDate();
#endif

#include "Matrix.h"
#include "Vector3D.h"

using namespace std;

class SparseMatrix;

class Vector
{
protected:
   typedef double Elm;
   int Cols_;
   Elm* Elements_;

public:
   // Constructor...
   // Precond. Vector object has been declared
   // Receive. Cols,Initial Value
   // Output. None
   // Postcond. Vector of size 1xCols allocated and
   // each element set to Inital Value (or not set at all)
   // Devaults: Cols-0,Initial Value= (un initialized)
   Vector(int const& Cols = 0);
   Vector(int const& Cols, Elm const& InitVal);
   Vector(Vector const& A);
   Vector(Vector3D const& A);
   Vector(Matrix const& A);
   friend Vector3D::Vector3D(Vector const& A);

   // Deconstructor...
   // Release dynamic memory.
   ~Vector();

   // Algebraic Operators...

   friend Vector& operator+(Vector& A)
   {
      return A;
   }

   friend Vector operator+(Vector const& A, Vector const& B);
   friend Vector operator-(Vector const& A, Vector const& B);
   friend Vector operator-(Vector const& A);

   // Dot Product
   friend Elm const operator*(Vector const& A, Vector const& B);

   // Cross Product (Cols_==3 Only)
   friend Vector operator%(Vector const& A, Vector const& B);
   friend Vector operator*(Matrix const& A, Vector const& B);
   friend Vector operator*(Vector const& A, Matrix const& B);
   friend Vector operator*(Elm const& A, Vector const& B);
   friend Vector operator*(Vector const& A, Elm const& B);
   friend Vector operator/(Vector const& A, Elm const& B);
   friend Vector operator*(SparseMatrix const& A, Vector const& B);
   friend Vector operator*(Vector const& A, SparseMatrix const& B);


   // Element Access methods
#ifdef CHECK_BOUNDS
   // With Bounds checking!!!
   Elm& operator[](int const& i);
   Elm const& operator[](int const& i) const;
#else
   // Without Bounds Checking!!!
   inline Elm& operator[](int const& i)
   {
      return Elements_[i];
   }

   inline Elm const& operator[](int const& i) const
   {
      return Elements_[i];
   }
#endif

   // Assignment Operatons

   Vector& operator=(Vector const& B);
   Vector& operator+=(Vector const& B)
   {
      return *this = *this + B;
   }

   Vector& operator-=(Vector const& B)
   {
      return *this = *this - B;
   }

   Vector& operator*=(Elm const& B)
   {
      return *this = *this * B;
   }

   Vector& operator/=(Elm const& B)
   {
      return *this = *this / B;
   }

   // Destructively Resize Vector
   // No change if size dosen't change
   void Resize(int const& Cols = 0);
   void Resize(int const& Cols, Elm const& InitVal);

   // Operations & Etc...
   int const& Dim() const
   {
      return Cols_;
   }

   // Standard IR^n Norm
   Matrix::Elm Norm() const;

   // Operations & Etc...

   // Solves Ax=B with pseudoinverse of A
   // using A=Q*R if A.Rows()==A.Cols()
   // using A^{+} = (A^{T}*A)^{-1}*A^{T}, with A=Q*R if A.Rows() > A.Cols()
   // using A^{+} = A^{T}*(A*A^{T})^{-1}, with A^{T}=Q*R if A.Rows() < A.Cols()
   friend void SolveQR(Matrix const& Q, Matrix const& R, Vector& x, Vector const& B);

   // Uses PLU decomposition with Forward and Backwards substitution
   friend Vector SolvePLU(Matrix const& A, Vector const& B);

   // Uses SVD decomposition.
   friend Vector SolveSVD(Matrix const& A, Vector const& B, Elm const& MaxCond = MAXCONDITION,
                          int const& PrintFlag = 0);

   // Perform Broyden update on QR factorization of A
   // it is expected that x.Norm() == 0
   friend void BroydenQRUpdate(Matrix& Q, Matrix& R, Vector const& y, Vector const& x);

   // Output/Input Function
   friend ostream& operator<<(ostream& out, Vector const& A);
   friend istream& operator>>(istream& in, Vector& A);

   static char const* const Revision();
};

#endif
