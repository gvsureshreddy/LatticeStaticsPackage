#ifndef __Math
#define __Math

#ifndef __MyMathBuildDate
#define __MyMathBuildDate
char *MyMathBuildDate();
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
void PolyRootsLaguerre(MyComplexDouble Coeff[],int Degree,MyComplexDouble Roots[],
		       int Polish);
// Lagurre's method to find a root of a polynomial near X
int Laguerre(MyComplexDouble Coeff[],int Degree,MyComplexDouble *X);

// Multiply polynomials A and B with resulting poly of degree DegA+DegB
void PolyMult(MyComplexDouble A[],int DegA,MyComplexDouble B[],int DegB,
	      MyComplexDouble Result[]);

#endif
