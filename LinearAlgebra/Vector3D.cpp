#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector3D.h"
#include "Vector.h"
#include "Matrix.h"

// Global IDString
char Vector3DID[]="$Id: Vector3D.cpp,v 1.4 2004/08/24 19:38:00 elliottr Exp $";

// Private Functions...

// Public Functions...

Vector3D::Vector3D(const Vector3D::Elm& InitVal)
{
   for (int i=0;i<V3DLEN;i++)
      Elements_[i]=InitVal;
}

Vector3D::Vector3D(const Vector3D::Elm& x,const Vector3D::Elm& y,const Vector3D::Elm& z)
{
   Elements_[0]=x;
   Elements_[1]=y;
   Elements_[2]=z;
}

Vector3D::Vector3D(const Vector3D& A)
{
   memmove(Elements_,A.Elements_,sizeof(Vector3D::Elm[V3DLEN]));

   return;
}

Vector3D::Vector3D(const Vector& A)
{
   if (A.Cols_ != V3DLEN)
   {
      cerr << "Vector3D: error: Constructor: Wrong Size" << endl;
      exit(-1);
   }
   
   memmove(Elements_,A.Elements_,sizeof(Vector3D::Elm[V3DLEN]));

   return;
}

Vector3D& Vector3D::operator=(const Vector3D& B)
{
   memmove(Elements_,B.Elements_,sizeof(Vector3D::Elm[V3DLEN]));

   return *this;
}

Vector3D operator*(const Vector3D& A,const Matrix& B)
{
   if (B.Rows_ != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(vec,matrix): Wrong Size" << endl;
      exit(-1);
   }

   Vector3D z(0.0);
   for (int i=0;i<V3DLEN;i++)
   {
      for (int j=0;j<V3DLEN;j++)
      {
	 z.Elements_[i] += A.Elements_[j]*B.Elements_[j][i];
      }
   }

   return z;
}

Vector3D operator*(const Matrix& A,const Vector3D& B)
{
   if (A.Cols_ != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(matrix,vec): Wrong Size" << endl;
      exit(-1);
   }

   Vector3D z(0.0);
   for (int i=0;i<V3DLEN;i++)
   {
      for (int j=0;j<V3DLEN;j++)
      {
	 z.Elements_[i] += A.Elements_[i][j]*B.Elements_[j];
      }
   }

   return z;
}

ostream& operator<<(ostream& out,const Vector3D& A)
{
   int W=out.width();

   for (register int i=0;i<V3DLEN;i++)
   {
      out << setw(W) << A.Elements_[i];
   }

   return out;
}

istream& operator>>(istream& in,Vector3D& A)
{
   for (register int i=0;i<V3DLEN;i++)
      in >> A.Elements_[i];

   return in;
}

char* Vector3D::Revision()
{
   return Vector3DID;
}
