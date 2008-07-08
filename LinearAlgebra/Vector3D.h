#ifndef __Vector3D
#define __Vector3D

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif


#include <cmath>

using namespace std;

#define V3DLEN 3

class SparseMatrix;
class Matrix;
class Vector;

class Vector3D
{
protected:
   typedef double Elm;
   Elm Elements_[V3DLEN];
   
public:

   // Constructor...
   Vector3D() {}
   Vector3D(Elm const& InitVal);
   Vector3D(Elm const& x,Elm const& y,Elm const& z);
   Vector3D(Vector3D const& A);
   Vector3D(Vector const& A);


   // Deconstructor...
   ~Vector3D() {}

   // Algebraic Operators...
   
   inline friend Vector3D& operator+(Vector3D& A) {return A;}
   inline friend Vector3D operator+(Vector3D const& A,Vector3D const& B)
   {
      return Vector3D(A[0]+B[0],
		      A[1]+B[1],
		      A[2]+B[2]);
   }
   inline friend Vector3D operator-(Vector3D const& A,Vector3D const& B)
   {
      return Vector3D(A[0]-B[0],
		      A[1]-B[1],
		      A[2]-B[2]);
   }
   inline friend Vector3D operator-(Vector3D const& A)
   {return Vector3D(-A[0],-A[1],-A[2]);}
   // Dot Product
   inline friend Elm operator*(Vector3D const& A,Vector3D const& B)
   {
      return (A[0]*B[0] +
	      A[1]*B[1] +
	      A[2]*B[2]);
   }
   // Cross Product
   inline friend Vector3D operator%(Vector3D const& A,Vector3D const& B)
   {
      return Vector3D(A[1]*B[2]-A[2]*B[1],
		      -(A[0]*B[2]-A[2]*B[0]),
		      A[0]*B[1]-A[1]*B[0]);
   }
   // Scalar Products
   inline friend Vector3D operator*(Elm const& A,Vector3D const& B)
   {return Vector3D(A*B[0],A*B[1],A*B[2]);}
   inline friend Vector3D operator*(Vector3D const& A,Elm const& B)
   {return Vector3D(B*A[0],B*A[1],B*A[2]);}
   inline friend Vector3D operator/(Vector3D const& A,Elm const& B)
   {return Vector3D(A[0]/B,A[1]/B,A[2]/B);}

   // Matrix Products
   friend Vector3D operator*(Vector3D const& A,Matrix const& B);
   friend Vector3D operator*(Matrix const& A,Vector3D const& B);
   friend Vector3D operator*(SparseMatrix const& A,Vector3D const& B);
   friend Vector3D operator*(Vector3D const& A,SparseMatrix const& B);	
   
   // Element Access methods
   inline Elm& operator[](int const& i) {return Elements_[i];}
   inline Elm const& operator[](int const& i) const {return Elements_[i];}

   // Assignment Operatons

   Vector3D& operator=(Vector3D const& B);
   inline Vector3D& operator+=(Vector3D const& B) {return *this=*this+B;}
   inline Vector3D& operator-=(Vector3D const& B) {return *this=*this-B;}
   inline Vector3D& operator*=(Vector3D const& B) {return *this=*this*B;}
   inline Vector3D& operator*=(Elm const& B) {return *this=*this*B;}

   inline Elm Norm() const {return sqrt((*this)*(*this));}

   // Output/Input Function
   friend ostream& operator<<(ostream& out,Vector3D const& A);
   friend istream& operator>>(istream& in,Vector3D& A);

   static char const* const Revision();
};

#endif
