#ifndef __Vector
#define __Vector

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif

#include "Matrix.h"
#include "Vector3D.h"

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
   Vector(const Vector3D& A);
   Vector(const Matrix& A);
   friend Vector3D::Vector3D(const Vector& A);

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
#ifdef CHECK_BOUNDS
   // With Bounds checking!!!
   Elm& operator[](const unsigned& i);
   const Elm operator[](const unsigned& i) const;
#else
   // Without Bounds Checking!!!
   inline Elm& operator[](const unsigned& i) {return Elements_[i];}
   inline const Elm operator[](const unsigned& i) const {return Elements_[i];}
#endif
   
   // Assignment Operatons

   Vector& operator=(const Vector& B);
   Vector operator+=(const Vector& B) {return *this=*this+B;}
   Vector operator-=(const Vector& B) {return *this=*this-B;}
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
   friend Vector SolvePLU(const Matrix& A,const Vector& B);
   // Uses SVD decomposition.
   friend Vector SolveSVD(const Matrix& A,const Vector& B,const Elm MaxCond=MAXCONDITION,
			  const int PrintFlag=0);

   // Output/Input Function
   friend ostream& operator<<(ostream& out,const Vector& A);
   friend istream& operator>>(istream& in,Vector& A);

   static char* Revision();
};

#endif
