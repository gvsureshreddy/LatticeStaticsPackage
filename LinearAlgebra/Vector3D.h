#ifndef __Vector3D
#define __Vector3D

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif


#include <math.h>

#define V3DLEN 3

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
   Vector3D(const Elm& InitVal);
   Vector3D(const Elm& x,const Elm& y,const Elm& z);
   Vector3D(const Vector3D& A);
   Vector3D(const Vector& A);


   // Deconstructor...
   ~Vector3D() {}

   // Algebraic Operators...
   
   inline friend Vector3D& operator+(Vector3D& A) {return A;}
   inline friend Vector3D operator+(const Vector3D& A,const Vector3D& B)
   {
      return Vector3D(A.Elements_[0]+B.Elements_[0],
		      A.Elements_[1]+B.Elements_[1],
		      A.Elements_[2]+B.Elements_[2]);
   }
   inline friend Vector3D operator-(const Vector3D& A,const Vector3D& B)
   {
      return Vector3D(A.Elements_[0]-B.Elements_[0],
		      A.Elements_[1]-B.Elements_[1],
		      A.Elements_[2]-B.Elements_[2]);
   }
   inline friend Vector3D operator-(const Vector3D& A)
   {return Vector3D(-A.Elements_[0],-A.Elements_[1],-A.Elements_[2]);}
   // Dot Product
   inline friend Elm operator*(const Vector3D& A,const Vector3D& B)
   {
      return (A.Elements_[0]*B.Elements_[0] +
	      A.Elements_[1]*B.Elements_[1] +
	      A.Elements_[2]*B.Elements_[2]);
   }
   // Cross Product
   inline friend Vector3D operator%(const Vector3D& A,const Vector3D& B)
   {
      return Vector3D(A.Elements_[1]*B.Elements_[2]-A.Elements_[2]*B.Elements_[1],
		      -(A.Elements_[0]*B.Elements_[2]-A.Elements_[2]*B.Elements_[0]),
		      A.Elements_[0]*B.Elements_[1]-A.Elements_[1]*B.Elements_[0]);
   }
   // Scalar Products
   inline friend Vector3D operator*(const Elm& A,const Vector3D& B)
   {return Vector3D(A*B.Elements_[0],A*B.Elements_[1],A*B.Elements_[2]);}
   inline friend Vector3D operator*(const Vector3D& A,const Elm& B)
   {return Vector3D(B*A.Elements_[0],B*A.Elements_[1],B*A.Elements_[2]);}
   inline friend Vector3D operator/(const Vector3D& A,const Elm& B)
   {return Vector3D(A.Elements_[0]/B,A.Elements_[1]/B,A.Elements_[2]/B);}

   // Matrix Products
   friend Vector3D operator*(const Vector3D& A,const Matrix& B);
   friend Vector3D operator*(const Matrix& A,const Vector3D& B);
   
   // Element Access methods
   inline Elm& operator[](const unsigned& i) {return Elements_[i];}
      
   inline const Elm operator[](const unsigned& i) const {return Elements_[i];}

   // Assignment Operatons

   Vector3D& operator=(const Vector3D& B);
   inline Vector3D operator+=(const Vector3D& B) {return *this=*this+B;}
   inline Vector3D operator-=(const Vector3D& B) {return *this=*this-B;}
   inline Vector3D operator*=(const Vector3D& B) {return *this=*this*B;}
   inline Vector3D operator*=(const Elm& B) {return *this=*this*B;}

   inline Elm Norm() {return sqrt((*this)*(*this));}

   // Output/Input Function
   friend ostream& operator<<(ostream& out,const Vector3D& A);
   friend istream& operator>>(istream& in,Vector3D& A);

   static char* Revision();
};

#endif
