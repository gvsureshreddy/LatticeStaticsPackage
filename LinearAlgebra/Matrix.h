#ifndef __Matrix
#define __Matrix

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif


#include <iostream>

using namespace std;

// Sentinal to allow for no initialization of data
#define SENTINAL -9999999.8888888887777788

// Maximum condition number for double precision
// above this a matrix is very ill-conditioned
//
// MAXCONDITION is used by SVD.
#define MAXCONDITION 10.0e12

class SparseMatrix;
class Vector;
class Vector3D;

class Matrix
{
public:
   typedef double Elm;

protected:

   Elm **Elements_;
   int Rows_;
   int Cols_;

   // Used by Det()
   Matrix Minor(int const& i,int const& j) const;

public:

   // Flag to print in mathematica format
   static int MathematicaPrintFlag;

   // Constructor...
   // Precond. Matrix object has been declared
   // Receive. Rows,Cols,Initial Value
   // Output.  None
   // Postcondition. Matrix of size RowsXCols
   //   allocated and each element set to Initial Value
   // Defaults: Rows=0,Cols=0,Initial Value= (Uninitialized)

   Matrix(int const& Rows=0,int const& Cols=0,Elm const& InitVal=SENTINAL);
   Matrix(Matrix const& A);

   // Deconstructor...
   ~Matrix();

   // Size Access...
   int const& Rows() const {return Rows_;}
   int const& Cols() const {return Cols_;}
   
   // Mathematical Operations...

   friend Matrix& operator+(Matrix& A) {return A;}
   friend Matrix operator+(Matrix const& A,Matrix const&B);
   friend Matrix operator-(Matrix const& A);
   friend Matrix operator-(Matrix const& A,Matrix const& B);
   friend Matrix operator*(Matrix const& A,Matrix const& B);
   friend Matrix operator*(Elm const& A,Matrix const& B);
   friend Matrix operator*(Matrix const& A,Elm const& B);
   friend Matrix operator*(SparseMatrix const& A,Matrix const& B);
   friend Matrix operator*(Matrix const& A,SparseMatrix const& B);
   // Below are defined in corresponding class --------------------
   friend Vector operator*(Matrix const& A,Vector const& B);
   friend Vector operator*(Vector const& A,Matrix const& B);
   friend Vector3D operator*(Matrix const& A,Vector3D const& B);
   friend Vector3D operator*(Vector3D const& A,Matrix const& B);
   // -------------------------------------------------------------
   friend Matrix operator/(Matrix const& A,Elm const& B);

   // Element Access methods

#ifdef CHECK_BOUNDS
   // Note: Index checking on Rows but not on Columns....
   Elm* const operator[](int const& i);
   Elm const* const operator[](int const& i) const;
#else
   // Note: NO Index checking
   inline Elm* const operator[](int const& i) {return Elements_[i];}
   inline Elm const* const operator[](int const& i) const {return Elements_[i];}
#endif
   
   // Assignment Operations

   Matrix& operator=(Matrix const& B);
   Matrix& operator+=(Matrix const& B) {return *this=*this+B;}
   Matrix& operator-=(Matrix const& B) {return *this=*this-B;}
   Matrix& operator*=(Matrix const& B) {return *this=*this*B;}
   Matrix& operator*=(Elm const& B)    {return *this=*this*B;}

   // Misc. Matrix Operatons
   
   Matrix& SetIdentity(int const& Size=0);
   Matrix Transpose() const;
   Matrix Inverse() const;
   int IsSquare() const {return Rows_==Cols_;}
   int IsNull() const {return (Rows_==0 || Cols_==0);}

   // Destructively Resize Matrix
   // No change if size does not change
   void Resize(int const& Rows=0,int const& Cols=0,Elm const& InitVal=SENTINAL);
   
   // Operations & Etc...

   // Deterimnent
   Elm Det() const;

   // Set P,L,U to the corresponding matricies of the PLU
   //   decomposition of A
   friend void PLU(Matrix const& A,Matrix& P,Matrix& L,Matrix& U);

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
   // each W[i][i] < fabs(MAX(W)) / MaxCond; will be set to 0.0
   // -- this most often reduces error when solving a system of equations
   // -- that is ill-conditioned (as compaired with a straight SVD or LU
   // -- decomposition).
   //
   // if (PrintFlag); then the condition number of A will be echoed on cerr
   //
   // WHENEVER a W[i][i] is explicitly set to 0.0 a message is echoed to cerr
   // NOTE: this situation may be detected by the calling program by compairing
   // -- the value of MaxCond with the returned condition number.
   //
   friend Elm SVD(Matrix const& A,Matrix& U,Matrix& W,Matrix& V,
                  Elm const& MaxCond=MAXCONDITION,int const& PrintFlag=0);

   // SymEigVal -- determine the eigenvalues of A
   // Diag(eigen values) = B.Transpose()*A*B
   //
   // Use Cyclic Jacobi Method -- Ref. "Numerical Analysis" by Pratel pg 440
   //
   // Returns vector containing eigenvalues
   // AND-- Matrix of eigenvectors (as columns) if a pointer is passed
   // MaxItr - max number of iterations
   // Tol - tolerance for convergence
   //
   // Note: Assumes A is SYMMETRIC
   friend Matrix SymEigVal(Matrix A,Matrix* const B=0,int const& MaxItr=100,
			   double const& Tol=1.0e-13);
   
   // Cholesky Decomposition of Matrix
   // A=U.Transpose()*D*U
   //
   // D - diagonal Matrix
   //
   // Assumes Symmetric Matrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void QR(Matrix const& A,Matrix& Q,Matrix& R,int const& CalcTranspose=0);
   // QR decomposition of A

   friend void Cholesky(Matrix const& A,Matrix& U,Matrix& D);

   // Return solution x of the linear system A*x=B
   // Uses PLU decomposition and Forward and Backwards substitution
   friend Matrix SolvePLU(Matrix const& A,Matrix const& B);
   
   // Return solution x of the linear system A*x=B
   // Uses SVD decomposition
   //
   // x = V*W.Inverse()*(U.Transpose()*B);
   // WHERE: W.Inverse() is actually calculated by hand and any
   // -- W[i][i] == 0.0 has inverse component 0.0
   friend Matrix SolveSVD(Matrix const& A,Matrix const& B,
			  Elm const& MaxCond=MAXCONDITION,
			  int const& PrintFlag=0);
   
   // Output/Input Functions
   friend ostream& operator<<(ostream& out,Matrix const& A);
   friend istream& operator>>(istream& in, Matrix& A);

   static char const* const Revision();
};

#endif
