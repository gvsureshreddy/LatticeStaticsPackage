#ifndef __CMatrix
#define __CMatrix

#include <iostream.h>
#include <complex.h>

// Sentinal to allow for no initialization of data
#define SENTINAL -9999999.8888888887777788

class Vector;
class Vector3D;
class Matrix;

class CMatrix
{
protected:
   typedef complex<double> Elm;

   Elm **Elements_;
   unsigned Rows_;
   unsigned Cols_;

   // Used by Det()
   CMatrix Minor(unsigned i,unsigned j) const;

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

   CMatrix(unsigned Rows=0,unsigned Cols=0,Elm InitVal=SENTINAL);
   CMatrix(const CMatrix& A);
   CMatrix(const Matrix& A);

   // Deconstructor...
   ~CMatrix();

   // Size Access...
   unsigned Rows() const {return Rows_;}
   unsigned Cols() const {return Cols_;}
   
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
   Elm* operator[](unsigned i);
   Elm* operator[](unsigned i) const;
#else
   // Note: NO Index checking
   Elm* operator[](unsigned i) {return Elements_[i];}
   Elm* operator[](unsigned i) const {return Elements_[i];}
#endif
   
   // Assignment Operations

   CMatrix& operator=(const CMatrix& B);
   CMatrix operator+=(const CMatrix& B) {return *this=*this+B;}
   CMatrix operator-=(const CMatrix& B) {return *this=*this-B;}
   CMatrix operator*=(const CMatrix& B) {return *this=*this*B;}
   CMatrix operator*=(const Elm& B)    {return *this=*this*B;}

   // Misc. CMatrix Operatons
   
   CMatrix& SetIdentity(unsigned Size=0);
   CMatrix Transpose() const;
   //CMatrix Inverse() const;
   int IsSquare() const {return Rows_==Cols_;}
   int IsNull() const {return (Rows_==0 || Cols_==0);}

   // Destructively Resize CMatrix
   // No change if size does not change
   void Resize(unsigned Rows=0,unsigned Cols=0,Elm InitVal=SENTINAL);
   
   // Operations & Etc...

   // Deterimnent
   Elm Det() const;

   // Cholesky Decomposition of CMatrix
   // A=U.Transpose()*D*U
   //
   // D - diagonal CMatrix
   //
   // Assumes Symmetric CMatrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void Cholesky(const CMatrix& A,CMatrix& U,CMatrix& D);

   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const CMatrix& A);
   friend istream& operator>>(istream& in, CMatrix& A);

   static char* Revision();
};

#endif
