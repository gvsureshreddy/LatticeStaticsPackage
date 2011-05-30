#ifndef RSE__CMatrix
#define RSE__CMatrix

#ifndef RSE__LinearAlgebraBuildDate
#define RSE__LinearAlgebraBuildDate
char* LinearAlgebraBuildDate();
#endif


#include <iostream>
#include <MyComplexDouble.h>
#include <Matrix.h>
#include <Vector.h>

using namespace std;

class Vector;
class Vector3D;

class CMatrix
{
protected:
   typedef MyComplexDouble Elm;

   Elm** Elements_;
   int Rows_;
   int Cols_;

   // Used by Det()
   CMatrix Minor(int const& i, int const& j) const;

public:
   // Flag to print in mathematica format
   static int MathematicaPrintFlag;

   // Constructor...
   // Precond. CMatrix object has been declared
   // Receive. Rows,Cols,Initial Value
   // Output.  None
   // Postcondition. CMatrix of size RowsXCols
   //   allocated and each element set to Initial Value
   // Defaults: Rows=0,Cols=0,Initial Value= (Uninitialized)

   CMatrix(int const& Rows = 0, int const& Cols = 0);
   CMatrix(int const& Rows, int const& Cols, Elm const& InitVal);
   CMatrix(CMatrix const& A);
   CMatrix(Matrix const& A);

   // Deconstructor...
   ~CMatrix();

   // Size Access...
   int const& Rows() const
   {
      return Rows_;
   }

   int const& Cols() const
   {
      return Cols_;
   }

   // Mathematical Operations...

   friend CMatrix& operator+(CMatrix& A)
   {
      return A;
   }

   friend CMatrix operator+(CMatrix const& A, CMatrix const& B);
   friend CMatrix operator-(CMatrix const& A);
   friend CMatrix operator-(CMatrix const& A, CMatrix const& B);
   friend CMatrix operator*(CMatrix const& A, CMatrix const& B);
   friend CMatrix operator*(Elm const& A, CMatrix const& B);
   friend CMatrix operator*(CMatrix const& A, Elm const& B);
   friend CMatrix operator*(double const& A, CMatrix const& B);
   friend CMatrix operator*(CMatrix const& A, double const& B);

   // Below are defined in corresponding class --------------------
   friend Vector operator*(CMatrix const& A, Vector const& B);
   friend Vector operator*(Vector const& A, CMatrix const& B);
   friend Vector3D operator*(CMatrix const& A, Vector3D const& B);
   friend Vector3D operator*(Vector3D const& A, CMatrix const& B);

   // -------------------------------------------------------------
   friend CMatrix operator/(CMatrix const& A, Elm const& B);

   // Element Access methods

#ifdef CHECK_BOUNDS
   // Note: Index checking on Rows but not on Columns....
   Elm const* const operator[](int const& i) const;
   Elm* const operator[](int const& i);
#else
   // Note: NO Index checking
   inline Elm const* const operator[](int const& i) const
   {
      return Elements_[i];
   }

   inline Elm* const operator[](int const& i)
   {
      return Elements_[i];
   }
#endif

   // Assignment Operations

   CMatrix& operator=(CMatrix const& B);
   CMatrix& operator+=(CMatrix const& B)
   {
      return *this = *this + B;
   }

   CMatrix& operator-=(CMatrix const& B)
   {
      return *this = *this - B;
   }

   CMatrix& operator*=(CMatrix const& B)
   {
      return *this = *this * B;
   }

   CMatrix& operator*=(Elm const& B)
   {
      return *this = *this * B;
   }

   // Misc. CMatrix Operatons

   CMatrix& SetIdentity(int const& Size = 0);
   CMatrix Transpose() const;
   CMatrix Conjugate() const;

   CMatrix ConjTrans() const
   {
      return (this->Transpose()).Conjugate();
   }

   CMatrix Inverse() const;
   int IsSquare() const
   {
      return Rows_ == Cols_;
   }

   int IsNull() const
   {
      return (Rows_ == 0 || Cols_ == 0);
   }

   // Destructively Resize CMatrix
   // No change if size does not change
   void Resize(int const& Rows = 0, int const& Cols = 0);
   void Resize(int const& Rows, int const& Cols, Elm const& InitVal);

   // Operations & Etc...

   // Deterimnent
   Elm Det() const;

   // Set P,L,U to the corresponding matricies of the PLU
   // decomposition of A
   friend void PLU(CMatrix const& A, CMatrix& P, CMatrix& L, CMatrix& U);

   // HermiteEigVal -- determine the eigenvalues of A
   // Diag(eigen values) = B.ConjTrans()*A*B
   //
   // Modified for Hermite matricies from symmetric jacobi in
   // "Matrix Computations"
   //
   // Returns vector containing eigenvalues
   // AND-- CMatrix of eigenvectors (as columns) if a pointer is passed
   // MaxItr - max number of iterations
   // Tol - tolerance for convergence
   //

   // Note: Assumes A is HERMITIAN
   friend Matrix HermiteEigVal(CMatrix A, CMatrix* const B = 0,
                               int const& MaxItr = 100, double const& Tol = 1.0e-13);

   // QR decomposition
   //
   // A   = Q*R  -- CalcTranspose = 0
   // A^T = Q*R  -- CalcTranspose = 1
   friend void QR(CMatrix const& A, CMatrix& Q, CMatrix& R, int const& CalcTranspose = 0);

   // Right Eigenvalues
   //
   // Find the right eigenvalues of a complex matrix using the QR algorithm
   friend CMatrix RightEigVals(CMatrix const& A, int const& MaxItr = 5000,
                               double const& Tol = 1.0e-13);

   // Cholesky Decomposition of CMatrix
   // A=U.ConjTrans()*D*U
   //
   // D - diagonal CMatrix of real values
   //
   // Assumes Hermitian CMatrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void Cholesky(CMatrix const& A, CMatrix& U, CMatrix& D);

   // Return solution x of the linear system A*x=B
   // Uses PLU decomposition and Forward and Backwards substitution
   friend CMatrix SolvePLU(CMatrix const& A, CMatrix const& B);
   friend CMatrix SolvePLU(CMatrix const& P, CMatrix const& L, CMatrix const& U,
                           CMatrix const& B);

   // Output/Input Functions
   friend ostream& operator<<(ostream& out, CMatrix const& A);
   friend istream& operator>>(istream& in, CMatrix& A);

   static char* Revision();
};

#endif
