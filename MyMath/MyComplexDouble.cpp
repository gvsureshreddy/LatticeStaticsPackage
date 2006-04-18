#include "MyComplexDouble.h"

// Global IDString
char MyComplexDoubleID[]="$Id: MyComplexDouble.cpp,v 1.2 2006/04/18 14:26:51 elliott Exp $";

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
      cerr << "Can not divide by complex (0,0)!" << endl;
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
      cerr << "Can not divide by zero in MyComplexDouble operator/" << endl;
      exit(-1);
   }

   return MyComplexDouble(A.Re_/B,A.Im_/B);
}

MyComplexDouble sqrt(MyComplexDouble A)
{
   double y,x;
   if (fabs(A.imag()) <= 2.0*EPS*fabs(A.real())) 
   {
      if (A.Re_ >= 0.0)
      {
	 return MyComplexDouble(sqrt(A.Re_),0.0);
      }
      else
      {
	 return MyComplexDouble(0.0,sqrt(-A.Re_));
      }
   }
   else
   {
      y = sqrt((A.mod() - A.Re_)/2.0);
      x = A.Im_/(2.0*y);
      return MyComplexDouble(x,y);
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
   
