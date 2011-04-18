#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "Vector3D.h"
#include "Vector.h"
#include "Matrix.h"

// Global IDString
char Vector3DID[] = "$Id: Vector3D.cpp,v 1.11 2011/04/18 16:31:46 elliott Exp $";

// Private Functions...

// Public Functions...

Vector3D::Vector3D(Vector3D::Elm const& InitVal)
{
   for (int i = 0; i < V3DLEN; i++)
   {
      Elements_[i] = InitVal;
   }
}

Vector3D::Vector3D(Vector3D::Elm const& x, Vector3D::Elm const& y, Vector3D::Elm const& z)
{
   Elements_[0] = x;
   Elements_[1] = y;
   Elements_[2] = z;
}

Vector3D::Vector3D(Vector3D const& A)
{
   memmove(Elements_, A.Elements_, sizeof(Vector3D::Elm[V3DLEN]));

   return;
}

Vector3D::Vector3D(Vector const& A)
{
   if (A.Cols_ != V3DLEN)
   {
      cerr << "Vector3D: error: Constructor: Wrong Size" << "\n";
      exit(-1);
   }

   memmove(Elements_, A.Elements_, sizeof(Vector3D::Elm[V3DLEN]));

   return;
}

Vector3D& Vector3D::operator=(Vector3D const& B)
{
   memmove(Elements_, B.Elements_, sizeof(Vector3D::Elm[V3DLEN]));

   return *this;
}

Vector3D operator*(Vector3D const& A, Matrix const& B)
{
   if (B.Rows_ != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(vec,matrix): Wrong Size" << "\n";
      exit(-1);
   }

   Vector3D z(0.0);
   for (int i = 0; i < V3DLEN; i++)
   {
      for (int j = 0; j < V3DLEN; j++)
      {
         z[i] += A[j] * B[j][i];
      }
   }

   return z;
}

Vector3D operator*(Matrix const& A, Vector3D const& B)
{
   if (A.Cols_ != V3DLEN)
   {
      cerr << "Vector3D: error: operator*(matrix,vec): Wrong Size" << "\n";
      exit(-1);
   }

   Vector3D z(0.0);
   for (int i = 0; i < V3DLEN; i++)
   {
      for (int j = 0; j < V3DLEN; j++)
      {
         z[i] += A[i][j] * B[j];
      }
   }

   return z;
}

ostream& operator<<(ostream& out, Vector3D const& A)
{
   int W = out.width();

   for (register int i = 0; i < V3DLEN; i++)
   {
      out << setw(W) << A[i];
   }

   return out;
}

istream& operator>>(istream& in, Vector3D& A)
{
   for (register int i = 0; i < V3DLEN; i++)
   {
      in >> A[i];
   }

   return in;
}

char const* const Vector3D::Revision()
{
   return Vector3DID;
}

