#ifndef __Matrix
#define __Matrix

#include <iostream.h>

// ***********************************************************************************
// $Log: Matrix.h,v $
// Revision 1.1  1999/08/13 14:53:49  elliottr
// Initial revision
//
// ***********************************************************************************

// Sentinal to allow for no initialization of data
#define SENTINAL -9999999.8888888887777788

class Vector;

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
   friend Vector operator*(const Matrix& A,const Vector& B);
   friend Vector operator*(const Vector& A,const Matrix& B);
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
   
   // Cholesky Decomposition of Matrix
   // A=U.Transpose()*D*U
   // Assumes Symmetric Matrix (thus uses only Upper Diagonal part of A
   // Note: will fail if A has EigenValue of 0.0
   friend void Cholesky(const Matrix& A,Matrix& U,Matrix& D);

   // Return solution x of the linear system A*x=B
   // Uses PLU decomposition and Forward and Backwards substitution
   friend Matrix Solve(const Matrix& A,const Matrix& B);

   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const Matrix& A);
   friend istream& operator>>(istream& in, Matrix& A);
};

#endif
