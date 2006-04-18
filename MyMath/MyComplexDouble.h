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

#define EPS 2.0e-12

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
   

   // Algebraic Operators...

   const friend MyComplexDouble& operator+(const MyComplexDouble& A) {return A;}
   friend MyComplexDouble operator-(const MyComplexDouble& A);
   friend MyComplexDouble operator+(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator-(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator*(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator/(const MyComplexDouble& A,const MyComplexDouble& B);
   friend MyComplexDouble operator/(const MyComplexDouble& A,const double& B);
   friend double abs(MyComplexDouble& A) {return A.mod();}
   friend MyComplexDouble sqrt(MyComplexDouble A);

   // Assignment Operators
   MyComplexDouble& operator=(const MyComplexDouble& B);
   MyComplexDouble operator+=(const MyComplexDouble& B) {return *this = *this+B;}
   MyComplexDouble operator-=(const MyComplexDouble& B) {return *this = *this-B;}
   MyComplexDouble operator*=(const MyComplexDouble& B) {return *this = *this*B;}
   MyComplexDouble operator/=(const MyComplexDouble& B) {return *this = *this/B;}

   // Output/Input Functions
   friend ostream& operator<<(ostream& out,const MyComplexDouble& A);
   friend istream& operator>>(istream& in,MyComplexDouble& A);

   static char* Revision();
};

#endif



   
