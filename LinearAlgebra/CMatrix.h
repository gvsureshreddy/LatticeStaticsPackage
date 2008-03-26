#ifndef __CMatrix
#define __CMatrix

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif


#include <iostream>
#include <MyComplexDouble.h>
#include <Matrix.h>

using namespace std;

// Sentinal to allow for no initialization of data
#define SENTINAL -9999999.8888888887777788

class Vector;
class Vector3D;

class CMatrix
{
protected:
   typedef MyComplexDouble Elm;

   Elm **Elements_;
   int Rows_;
   int Cols_;

   // Used by Det()
   CMatrix Minor(int i,int j) const;

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

   CMatrix(int Rows=0,int Cols=0,Elm InitVal=SENTINAL);
   CMatrix(const CMatrix& A);
   CMatrix(const Matrix& A);

   // Deconstructor...
   ~CMatrix();

   // Size Access...
   int Rows() const {return Rows_;}
   int Cols() const {return Cols_;}
   
   // Mathematical Operations...

   friend CMatrix& operator+(CMatrix& A) {return A;}
   friend CMatrix operator+(const CMatrix& A,const CMatrix&B);
   friend CMatrix operator-(const CMatrix& A);
   friend CMatrix operator-(const CMatrix& A,const CMatrix& B);
   friend CMatrix operator*(const CMatrix& A,const CMatrix& B);
   friend CMatrix operator*(const Elm& A,const CMatrix& B);
   friend CMatrix operator*(const CMatrix& A,const Elm& B);
   // Below are defined in corresponding class --------------------
   friend Vector operator*(const CMatrix& A,const Vector& B);
   friend Vector operator*(const Vector& A,const CMatrix& B);
   friend Vector3D operator*(const CMatrix& A,const Vector3D& B);
   friend Vector3D operator*(const Vector3D& A,const CMatrix& B);
   // -------------------------------------------------------------
   friend CMatrix operator/(const CMatrix& A,const Elm& B);

   // Element Access methods

#ifdef CHECK_BOUNDS
   // Note: Index checking on Rows but not on Columns....
   Elm* operator[](int i);
   Elm* operator[](int i) const;
#else
   // Note: NO Index checking
   Elm* operator[](int i) {return Elements_[i];}
   Elm* operator[](int i) const {return Elements_[i];}
#endif
   
   // Assignment Operations

   CMatrix& operator=(const CMatrix& B);
   CMatrix operator+=(const CMatrix& B) {return *this=*this+B;}
   CMatrix operator-=(const CMatrix& B) {return *this=*this-B;}
   CMatrix operator*=(const CMatrix& B) {return *this=*this*B;}
   CMatrix operator*=(const Elm& B)    {return *this=*this*B;}

   // Misc. CMatrix Operatons
   
   CMatrix& SetIdentity(int Size=0);
   CMatrix Transpose() const;
   CMatrix Conjugate() const;
   CMatrix ConjTrans() {return (this->Transpose()).Conjugate();}
   CMatrix Inverse() const;
   int IsSquare() const {return Rows_==Cols_;}
   int IsNull() const {return (Rows_==0 || Cols_==0);}

   // Destructively Resize CMatrix
   // No change if size does not change
   void Resize(int Rows=0,int Cols=0,Elm InitVal=SENTINAL);
   
   // Operations & Etc...

   // Deterimnent
   Elm Det() const;
   // Set P,L,U to the corresponding matricies of the PLU
   //   decomposition of A
   friend void PLU(const CMatrix& A,CMatrix& P,CMatrix& L,CMatrix& U);

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
   friend Matrix HermiteEigVal(CMatrix A,CMatrix *B=NULL,const int MaxItr=100,
			       const double Tol=1.0e-13);

   // Cholesky Decomposition of CMatrix
   // A=U.ConjTrans()*D*U
   //
   // D - diagonal CMatrix of real values
   //
   // Assumes Hermitian CMatrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void Cholesky(const CMatrix& A,CMatrix& U,CMatrix& D);

   // Return solution x of the linear system A*x=B
   // Uses PLU decomposition and Forward and Backwards substitution
   friend CMatrix SolvePLU(const CMatrix& A,const CMatrix& B);
   
   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const CMatrix& A);
   friend istream& operator>>(istream& in, CMatrix& A);

   static char* Revision();
};

#endif
