#include <string.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <Vector3D.h>

// Global IDString
char Vector3DID[]="$Id: Vector3D.cpp,v 1.1 2000/09/24 06:54:38 elliottr Exp $";

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

Vector3D& Vector3D::operator=(const Vector3D& B)
{
   memmove(Elements_,B.Elements_,sizeof(Vector3D::Elm[V3DLEN]));

   return *this;
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
