#ifndef __Matrix
#define __Matrix

#include <iostream.h>

// Sentinal to allow for no initialization of data
#define SENTINAL -9999999.8888888887777788

// Maximum condition number for double precision
// above this a matrix is very ill-conditioned
//
// MAXCONDITION is used by SVD.
#define MAXCONDITION 10.0e12


class Vector;
class Vector3D;

class Matrix
{
protected:
   typedef double Elm;

   Elm **Elements_;
   unsigned Rows_;
   unsigned Cols_;

   // Used by Det()
   Matrix Minor(unsigned i,unsigned j) const;

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

   Matrix(unsigned Rows=0,unsigned Cols=0,Elm InitVal=SENTINAL);
   Matrix(const Matrix& A);

   // Deconstructor...
   ~Matrix();

   // Size Access...
   unsigned Rows() const {return Rows_;}
   unsigned Cols() const {return Cols_;}
   
   // Mathematical Operations...

   friend Matrix& operator+(Matrix& A) {return A;}
   friend Matrix operator+(const Matrix& A,const Matrix&B);
   friend Matrix operator-(const Matrix& A);
   friend Matrix operator-(const Matrix& A,const Matrix& B);
   friend Matrix operator*(const Matrix& A,const Matrix& B);
   friend Matrix operator*(const Elm& A,const Matrix& B);
   friend Matrix operator*(const Matrix& A,const Elm& B);
   // Below are defined in corresponding class --------------------
   friend Vector operator*(const Matrix& A,const Vector& B);
   friend Vector operator*(const Vector& A,const Matrix& B);
   friend Vector3D operator*(const Matrix& A,const Vector3D& B);
   friend Vector3D operator*(const Vector3D& A,const Matrix& B);
   // -------------------------------------------------------------
   friend Matrix operator/(const Matrix& A,const Elm& B);

   // Element Access methods

   // Note: Index checking on Rows but not on Columns....
   Elm* operator[](unsigned i);
   Elm* operator[](unsigned i) const;

   // Assignment Operations

   Matrix& operator=(const Matrix& B);
   Matrix operator+=(const Matrix& B) {return *this=*this+B;}
   Matrix operator-=(const Matrix& B) {return *this=*this-B;}
   Matrix operator*=(const Matrix& B) {return *this=*this*B;}
   Matrix operator*=(const Elm& B)    {return *this=*this*B;}

   // Misc. Matrix Operatons
   
   Matrix& SetIdentity(unsigned Size=0);
   Matrix Transpose() const;
   Matrix Inverse() const;
   int IsSquare() const {return Rows_==Cols_;}
   int IsNull() const {return (Rows_==0 || Cols_==0);}

   // Destructively Resize Matrix
   // No change if size does not change
   void Resize(unsigned Rows=0,unsigned Cols=0,Elm InitVal=SENTINAL);
   
   // Operations & Etc...

   // Deterimnent
   Elm Det() const;

   // Set P,L,U to the corresponding matricies of the PLU
   //   decomposition of A
   friend void PLU(const Matrix& A,Matrix& P,Matrix& L,Matrix& U);

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
   friend Elm SVD(const Matrix& A,Matrix& U,Matrix& W,Matrix& V,
		  const Elm MaxCond=MAXCONDITION,const int PrintFlag=0);

   // SymEigVal -- determine the eigenvalues of A
   // Use Cyclic Jacobi Method -- Ref. "Numerical Analysis" by Pratel pg 440
   //
   // Returns vector containing eigenvalues
   // MaxItr - max number of iterations
   // Tol - tolerance for convergence
   //
   // Note: Assumes A is SYMMETRIC
   friend Matrix SymEigVal(Matrix A,const int MaxItr=100,const double Tol=1.0e-12);
   
   // Cholesky Decomposition of Matrix
   // A=U.Transpose()*D*U
   //
   // D - diagonal Matrix
   //
   // Assumes Symmetric Matrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void Cholesky(const Matrix& A,Matrix& U,Matrix& D);

   // Return solution x of the linear system A*x=B
   // Uses PLU decomposition and Forward and Backwards substitution
   friend Matrix SolvePLU(const Matrix& A,const Matrix& B);
   
   // Return solution x of the linear system A*x=B
   // Uses SVD decomposition
   //
   // x = V*W.Inverse()*(U.Transpose()*B);
   // WHERE: W.Inverse() is actually calculated by hand and any
   // -- W[i][i] == 0.0 has inverse component 0.0
   friend Matrix SolveSVD(const Matrix& A,const Matrix& B,
			  const Elm MaxCond=MAXCONDITION,
			  const int PrintFlag=0);
   
   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const Matrix& A);
   friend istream& operator>>(istream& in, Matrix& A);

   static char* Revision();
};

#endif
