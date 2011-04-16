#ifndef RSE__Math
#define RSE__Math

#ifndef RSE__MyMathBuildDate
#define RSE__MyMathBuildDate
char* MyMathBuildDate();
#endif


#include <iostream>
#include "MyComplexDouble.h"

using namespace std;

// This is to be a math library for misc routines
// that I have needed.  It is likely that many of them
// will come from Numerical Receipes and various other
// places.


// Routine to find roots of a polynomial
// By Laguerre's method
// Taken from Numerical Recipes
#define MAXDEGREE 100
#define EPS 2.0e-12
//
// Coeff[0..Degree] define the poly sum(i=0)(m) Coeff[i]x^i
void PolyRootsLaguerre(MyComplexDouble const* const Coeff, int const& Degree,
                       MyComplexDouble* const Roots, int const& Polish);
// Lagurre's method to find a root of a polynomial near X
int Laguerre(MyComplexDouble const* const Coeff, int const& Degree,
             MyComplexDouble* const X);

// Multiply polynomials A and B with resulting poly of degree DegA+DegB
void PolyMult(MyComplexDouble const* const A, int const& DegA, MyComplexDouble const* const B,
              int const& DegB, MyComplexDouble* const Result);

#endif

