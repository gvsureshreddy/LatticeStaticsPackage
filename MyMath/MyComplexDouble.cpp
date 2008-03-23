#include "MyComplexDouble.h"

// Global IDString
char MyComplexDoubleID[]="$Id: MyComplexDouble.cpp,v 1.4 2008/03/23 02:25:51 elliott Exp $";

// Private Functions...

// Public Functions...

MyComplexDouble& MyComplexDouble::operator=(const MyComplexDouble& A)
{
   Re_ = A.Re_;
   Im_ = A.Im_;
   
   return *this;
}

MyComplexDouble operator-(const MyComplexDouble& A)
{
   return MyComplexDouble(-A.Re_,-A.Im_);
}

MyComplexDouble operator+(const MyComplexDouble& A,const MyComplexDouble& B)
{
   return MyComplexDouble(A.Re_+B.Re_,A.Im_+B.Im_);
}

MyComplexDouble operator-(const MyComplexDouble& A,const MyComplexDouble& B)
{
   return MyComplexDouble(A.Re_-B.Re_,A.Im_-B.Im_);
}

MyComplexDouble operator*(const MyComplexDouble& A,const MyComplexDouble& B)
{
   return MyComplexDouble(A.Re_*B.Re_-A.Im_*B.Im_, A.Im_*B.Re_+A.Re_*B.Im_);
}

MyComplexDouble operator*(const double& A,const MyComplexDouble& B)
{
   return MyComplexDouble(A*B.Re_, A*B.Im_);
}

MyComplexDouble operator*(const MyComplexDouble& A,const double& B)
{
   return MyComplexDouble(A.Re_*B, A.Im_*B);
}

MyComplexDouble operator/(const MyComplexDouble& A,const MyComplexDouble& B)
{
   if ((B.Re_ == 0.0) && (B.Im_ == 0.0))
   {
      cerr << "Can not divide by complex (0,0)!" << "\n";
      exit(-1);
   }

   double det = B.Re_*B.Re_ + B.Im_*B.Im_;
   return MyComplexDouble(
      (A.Re_*B.Re_ + A.Im_*B.Im_)/det,
      (A.Im_*B.Re_ - A.Re_*B.Im_)/det);
}

MyComplexDouble operator/(const MyComplexDouble& A,const double& B)
{
   if (B == 0.0)
   {
      cerr << "Can not divide by zero in MyComplexDouble operator/" << "\n";
      exit(-1);
   }

   return MyComplexDouble(A.Re_/B,A.Im_/B);
}

MyComplexDouble sqrt(const MyComplexDouble& A)
{
   double
      x=A.Re_,
      y=A.Im_;

   if (x == double())
   {
      double t = sqrt(abs(y)/2.0);
      return MyComplexDouble(t,y<double()?-t:t);
   }
   else
   {
      double t = sqrt(2.0*(abs(A) + abs(x)));
      double u = t/2.0;
      return x>double()
	 ? MyComplexDouble(u,y/t)
	 : MyComplexDouble(abs(y)/t, y < double() ? -u : u);
   }
}

ostream& operator<<(ostream& out, const MyComplexDouble& A)
{
   int W=out.width();
   out.width(0);

   out << "(" << setw(W) << A.Re_
       << "," << setw(W) << A.Im_
       << ")";

   return out;
}

istream & operator>>(istream& in,MyComplexDouble& A)
{
   char ch;
   double re=0.0,im=0.0;

   in >> ch;
   if (ch =='(')
   {
      in >> re >> ch;
      if (ch == ',')
      {
	 in >> im >> ch;
	 if (ch == ')')
	 {
	    A.Re_=re;
	    A.Im_=im;
	 }
	 else
	    in.setstate(ios_base::failbit);
      }
      else if (ch == ')')
      {
	 A.Re_ = re;
	 A.Im_ = 0.0;
      }
      else
      {
	 in.setstate(ios_base::failbit);
      }
   }
   else
   {
      in.putback(ch);
      in >> re;
      A.Re_=re;
      A.Im_=0.0;
   }

   return in;
}

char* MyComplexDouble::Revision()
{
   return MyComplexDoubleID;
}
   
