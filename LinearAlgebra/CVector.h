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
   typedef complex<double> Elm;
   unsigned Cols_;
   Elm *Elements_;
   
public:

   // Constructor...
   // Precond. CVector object has been declared
   // Receive. Cols,Initial Value
   // Output. None
   // Postcond. CVector of size 1xCols allocated and
   // each element set to Inital Value (or not set at all)
   // Devaults: Cols-0,Initial Value= (un initialized)
   CVector(const unsigned& Cols=0,const Elm& InitVal=SENTINAL);
   CVector(const CVector& A);
   CVector(const CMatrix& A);
   CVector(const Vector& A);
   CVector(const Matrix& A);

   // Deconstructor...
   // Release dynamic memory.
   ~CVector();

   // Algebraic Operators...
   
   friend CVector& operator+(CVector& A) {return A;}
   friend CVector operator+(const CVector& A,const CVector& B);
   friend CVector operator-(const CVector& A,const CVector& B);
   friend CVector operator-(const CVector& A);
   // Dot Product
   friend Elm operator*(const CVector& A,const CVector& B);
   friend CVector operator*(const CMatrix& A,const CVector& B);
   friend CVector operator*(const CVector& A,const CMatrix& B);
   friend CVector operator*(const Elm& A,const CVector& B);
   friend CVector operator*(const CVector& A,const Elm& B);
   friend CVector operator/(const CVector& A,const Elm& B);
   
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

   CVector& operator=(const CVector& B);
   CVector operator+=(const CVector& B) {return *this=*this+B;}
   CVector operator-=(const CVector& B) {return *this=*this-B;}
   CVector operator*=(const Elm& B) {return *this=(*this)*B;}

   // Destructively Resize CVector
   // No change if size dosen't change
   void Resize(const unsigned& Cols=0,const Elm& InitVal=SENTINAL);

   // Operations & Etc...
   unsigned Dim() const {return Cols_;}
   // Standard IC^n Norm
   Elm Norm();

   // Operations & Etc...

   // Uses PLU decomposition with Forward and Backwards substitution
   friend CVector SolvePLU(const CMatrix& A,const CVector& B);

   // Output/Input Function
   friend ostream& operator<<(ostream& out,const CVector& A);
   friend istream& operator>>(istream& in,CVector& A);

   static char* Revision();
};

#endif
