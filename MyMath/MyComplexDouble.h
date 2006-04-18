#ifndef __MyComplexDouble
#define __MyComplexDouble

#ifndef __MyMathBuildDate
#define __MyMathBuildDate
char *MyMathBuildDate();
#endif

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class MyComplexDouble
{
protected:
   typedef double Elm;
   Elm Re_;
   Elm Im_;

public:

   // Constructors
   MyComplexDouble(double Re=0.0, double Im=0.0): Re_(Re), Im_(Im) {}
   MyComplexDouble(const MyComplexDouble &CDB) {Re_=CDB.Re_; Im_=CDB.Im_;}

   // Destructor
   ~MyComplexDouble() {}

   // Data Acces methods
   Elm real() { return Re_;}
   Elm imag() { return Im_;}
   Elm mod() { return sqrt(Re_*Re_ + Im_*Im_);}
   Elm arg() { return (Re_>0.0)?atan(Re_/Im_):4.0*atan(1.0)-atan(Re_/Im_);}

   // Comparison methods...

   int operator==(const MyComplexDouble &right) const
   {return ((Re_==right.Re_) && (Im_==right.Im_));}
   int operator!=(const MyComplexDouble &right) const
   {return ((Re_!=right.Re_) || (Im_!=right.Im_));}

   // Algebraic Operators...

   MyComplexDouble conj() { return MyComplexDouble(Re_,-Im_);}
   const friend MyComplexDouble& operator+(const MyComplexDouble& A) {return A;}
   friend MyComplexDouble operator-(const MyComplexDouble& A);
   friend MyComplexDouble operator+(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator-(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator*(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator*(const double& A,const MyComplexDouble& B);
   friend MyComplexDouble operator*(const MyComplexDouble& A,const double& B);
   friend MyComplexDouble operator/(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator/(const MyComplexDouble& A,const double& B);
   friend double abs(MyComplexDouble A) {return A.mod();}
   friend MyComplexDouble sqrt(const MyComplexDouble& A);

   // Elementary Functions...
   friend MyComplexDouble exp(const MyComplexDouble &A)
   {return MyComplexDouble(exp(A.Re_)*cos(A.Im_),exp(A.Re_)*sin(A.Im_));}

   // Assignment Operators
   MyComplexDouble& operator=(const MyComplexDouble& B);
   MyComplexDouble operator+=(const MyComplexDouble& B) {return *this = *this+B;}
   MyComplexDouble operator-=(const MyComplexDouble& B) {return *this = *this-B;}
   MyComplexDouble operator*=(const MyComplexDouble& B) {return *this = *this*B;}
   MyComplexDouble operator*=(const double& B) {return *this = *this*B;}
   MyComplexDouble operator/=(const MyComplexDouble& B) {return *this = *this/B;}
   MyComplexDouble operator/=(const double& B) {return *this = *this/B;}

   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const MyComplexDouble& A);
   friend istream& operator>>(istream& in,MyComplexDouble& A);

   static char* Revision();
};

#endif



   
