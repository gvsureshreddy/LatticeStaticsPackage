#ifndef RSE__MyComplexDouble
#define RSE__MyComplexDouble

#ifndef RSE__MyMathBuildDate
#define RSE__MyMathBuildDate
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
   MyComplexDouble(double const& Re=0.0, double const& Im=0.0): Re_(Re), Im_(Im) {}
   MyComplexDouble(MyComplexDouble const& CDB) {Re_=CDB.Re_; Im_=CDB.Im_;}

   // Destructor
   ~MyComplexDouble() {}

   // Data Acces methods
   Elm const& real() const { return Re_;}
   Elm const& imag() const { return Im_;}
   Elm mod() const { return sqrt(Re_*Re_ + Im_*Im_);}
   Elm arg() const { return (Re_>0.0)?atan(Re_/Im_):4.0*atan(1.0)-atan(Re_/Im_);}

   // Comparison methods...

   int operator==(MyComplexDouble const& right) const
   {return ((Re_==right.Re_) && (Im_==right.Im_));}
   int operator!=(MyComplexDouble const& right) const
   {return ((Re_!=right.Re_) || (Im_!=right.Im_));}

   // Algebraic Operators...

   MyComplexDouble conj() const { return MyComplexDouble(Re_,-Im_);}
   friend MyComplexDouble& operator+(MyComplexDouble& A) {return A;}
   friend MyComplexDouble operator-(MyComplexDouble const& A);
   friend MyComplexDouble operator+(MyComplexDouble const& A,
                                    MyComplexDouble const& B);
   friend MyComplexDouble operator-(MyComplexDouble const& A,
                                    MyComplexDouble const& B);
   friend MyComplexDouble operator*(MyComplexDouble const& A,
                                    MyComplexDouble const& B);
   friend MyComplexDouble operator*(double const& A,MyComplexDouble const& B);
   friend MyComplexDouble operator*(MyComplexDouble const& A,double const& B);
   friend MyComplexDouble operator/(MyComplexDouble const& A,
                                    MyComplexDouble const& B);
   friend MyComplexDouble operator/(MyComplexDouble const& A,double const& B);
   friend MyComplexDouble operator/(double const& A,MyComplexDouble const& B)
   {return MyComplexDouble(A)/B;}
   friend double abs(MyComplexDouble const& A) {return A.mod();}
   friend MyComplexDouble sqrt(MyComplexDouble const& A);

   // Elementary Functions...
   friend MyComplexDouble exp(MyComplexDouble const& A)
   {return MyComplexDouble(exp(A.Re_)*cos(A.Im_),exp(A.Re_)*sin(A.Im_));}

   // Assignment Operators
   MyComplexDouble& operator=(MyComplexDouble const& B);
   MyComplexDouble& operator+=(MyComplexDouble const& B) {return *this = *this+B;}
   MyComplexDouble& operator-=(MyComplexDouble const& B) {return *this = *this-B;}
   MyComplexDouble& operator*=(MyComplexDouble const& B) {return *this = *this*B;}
   MyComplexDouble& operator*=(double const& B) {return *this = *this*B;}
   MyComplexDouble& operator/=(MyComplexDouble const& B) {return *this = *this/B;}
   MyComplexDouble& operator/=(double const& B) {return *this = *this/B;}

   // Output/Input Functions
   friend ostream& operator<<(ostream& out,MyComplexDouble const& A);
   friend istream& operator>>(istream& in,MyComplexDouble& A);

   static char const* const Revision();
};

#endif
