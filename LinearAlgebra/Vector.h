#ifndef __Vector
#define __Vector

#include "Matrix.h"

// ***********************************************************************************
// $Log: Vector.h,v $
// Revision 1.1  1999/08/13 14:54:07  elliottr
// Initial revision
//
// Revision 1.1  1999/07/26 16:12:11  elliottr
// Initial revision
//
// ***********************************************************************************

// Sentinal Value to allow conditional initialization of data
#define SENTINAL -9999999.8888888887777788

class Vector
{
protected:
   typedef double Elm;
   unsigned Cols_;
   Elm *Elements_;
   
public:

   // Constructor...
   // Precond. Vector object has been declared
   // Receive. Cols,Initial Value
   // Output. None
   // Postcond. Vector of size 1xCols allocated and
   // each element set to Inital Value (or not set at all)
   // Devaults: Cols-0,Initial Value= (un initialized)
   Vector(const unsigned& Cols=0,const Elm& InitVal=SENTINAL);
   Vector(const Vector& A);
   Vector(const Matrix& A);

   // Deconstructor...
   // Release dynamic memory.
   ~Vector();

   // Algebraic Operators...
   
   friend Vector& operator+(Vector& A) {return A;}
   friend Vector operator+(const Vector& A,const Vector& B);
   friend Vector operator-(const Vector& A,const Vector& B);
   friend Vector operator-(const Vector& A);
   // Dot Product
   friend Elm operator*(const Vector& A,const Vector& B);
   // Cross Product (Cols_==3 Only)
   friend Vector operator%(const Vector& A,const Vector& B);
   friend Vector operator*(const Matrix& A,const Vector& B);
   friend Vector operator*(const Vector& A,const Matrix& B);
   friend Vector operator*(const Elm& A,const Vector& B);
   friend Vector operator*(const Vector& A,const Elm& B);
   friend Vector operator/(const Vector& A,const Elm& B);

   // Element Access methods
   // With Bounds checking!!!
   Elm& operator[](const unsigned& i);
   const Elm operator[](const unsigned& i) const;

   // Assignment Operatons

   Vector& operator=(const Vector& B);
   Vector operator+=(const Vector& B) {return *this=*this+B;}
   Vector operator-=(const Vector& B) {return *this=*this-B;}
   Vector operator*=(const Vector& B) {return *this=*this*B;}
   Vector operator*=(const Elm& B) {return *this=*this*B;}

   // Destructively Resize Vector
   // No change if size dosen't change
   void Resize(const unsigned& Cols=0,const Elm& InitVal=SENTINAL);

   // Operations & Etc...
   unsigned Dim() const {return Cols_;}
   // Standard IR^n Norm
   Matrix::Elm Norm();

   // Operations & Etc...

   // Uses PLU decomposition with Forward and Backwards substitution
   friend Vector Solve(const Matrix& A,const Vector& B);

   // Output/Input Function
   friend ostream& operator<<(ostream& out,const Vector& A);
   friend istream& operator>>(istream& in,Vector& A);
};

#endif
