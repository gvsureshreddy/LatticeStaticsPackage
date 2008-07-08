#ifndef __CVector
#define __CVector

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif


#include "CMatrix.h"

using namespace std;

// Sentinal Value to allow conditional initialization of data
#define SENTINAL -9999999.8888888887777788

class Vector;
class Matrix;

class CVector
{
protected:
   typedef MyComplexDouble Elm;
   int Cols_;
   Elm *Elements_;
   
public:

   // Constructor...
   // Precond. CVector object has been declared
   // Receive. Cols,Initial Value
   // Output. None
   // Postcond. CVector of size 1xCols allocated and
   // each element set to Inital Value (or not set at all)
   // Devaults: Cols-0,Initial Value= (un initialized)
   CVector(int const& Cols=0,Elm const& InitVal=SENTINAL);
   CVector(CVector const& A);
   CVector(CMatrix const& A);
   CVector(Vector const& A);
   CVector(Matrix const& A);

   // Deconstructor...
   // Release dynamic memory.
   ~CVector();

   // Algebraic Operators...
   
   friend CVector& operator+(CVector& A) {return A;}
   friend CVector operator+(CVector const& A,CVector const& B);
   friend CVector operator-(CVector const& A,CVector const& B);
   friend CVector operator-(CVector const& A);
   // Dot Product
   friend Elm const operator*(CVector const& A,CVector const& B);
   friend CVector operator*(CMatrix const& A,CVector const& B);
   friend CVector operator*(CVector const& A,CMatrix const& B);
   friend CVector operator*(Elm const& A,CVector const& B);
   friend CVector operator*(CVector const& A,Elm const& B);
   friend CVector operator/(CVector const& A,Elm const& B);
   
   // Element Access methods
#ifdef CHECK_BOUNDS
   // With Bounds checking!!!
   Elm& operator[](int const& i);
   Elm const& operator[](int const& i) const;
#else
   // Without Bounds Checking!!!
   inline Elm& operator[](int const& i) {return Elements_[i];}
   inline Elm const& operator[](int const& i) const {return Elements_[i];}
#endif
   
   // Assignment Operatons

   CVector& operator=(CVector const& B);
   CVector& operator+=(CVector const& B) {return *this=*this+B;}
   CVector& operator-=(CVector const& B) {return *this=*this-B;}
   CVector& operator*=(Elm const& B) {return *this=(*this)*B;}

   // Destructively Resize CVector
   // No change if size dosen't change
   void Resize(int const& Cols=0,Elm const& InitVal=SENTINAL);

   // Operations & Etc...
   int const& Dim() const {return Cols_;}
   // Standard IC^n Norm
   Elm Norm() const;

   // Operations & Etc...

   // Uses PLU decomposition with Forward and Backwards substitution
   friend CVector SolvePLU(CMatrix const& A,CVector const& B);

   // Output/Input Function
   friend ostream& operator<<(ostream& out,CVector const& A);
   friend istream& operator>>(istream& in,CVector& A);

   static char const* const Revision();
};

#endif
