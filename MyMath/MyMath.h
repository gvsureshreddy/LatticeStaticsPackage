#ifndef __Math
#define __Math

#include <iostream.h>
#include <complex.h>

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
void PolyRootsLaguerre(double_complex Coeff[],int Degree,double_complex Roots[],
		       int Polish);
// Lagurre's method to find a root of a polynomial near X
int Laguerre(double_complex Coeff[],int Degree,double_complex *X);

// Multiply polynomials A and B with resulting poly of degree DegA+DegB
void PolyMult(double_complex A[],int DegA,double_complex B[],int DegB,
	      double_complex Result[]);

#endif
