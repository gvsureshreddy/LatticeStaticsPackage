#include <iomanip.h>
#include "Math.h"

void PolyRootsLaguerre(double_complex Coeff[],int Degree,double_complex Roots[],
		       int Polish)
{
   int i,its,j,jj;
   double_complex x,b,c,ad[MAXDEGREE];

   // Copy of coefficients for successive deflation.
   for (j=0;j<=Degree;j++) ad[j] = Coeff[j];
   // Loop over each root to be found.
   for (j=Degree;j>=1;j--)
   {
      // start at zero to favor convergence to smallest
      // remaining root
      x = double_complex(0.0,0.0);
      // find the root
      its=Laguerre(ad,j,&x);
      if (fabs(x.imag()) <= 2.0*EPS*fabs(x.real())) x=double_complex(x.real());
      Roots[j-1] = x;
      b = ad[j];  // forward deflation.
      for (jj=j-1;jj>=0;jj--)
      {
	 c = ad[jj];
	 ad[jj] = b;
	 b = x*b+c;
      }
   }

   if (Polish)
   {
      for (j=1;j<=Degree;j++)
	 its=Laguerre(Coeff,Degree,&Roots[j-1]);
   }

   for (j=1;j<=Degree-1;j++)
   {
      x=Roots[j];
      for (i=j-1;i>=0;i--)
      {
	 if (Roots[i].real() <= x.real()) break;
	 Roots[i+1]=Roots[i];
      }
      Roots[i+1]=x;
   }
}


#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define EPSS 1.0e-14
#define FMAX(a,b) (a >= b ? a : b)

int Laguerre(double_complex Coeff[],int Degree,double_complex *X)
{
   int iter,j;
   double abx,abp,abm,err;
   double_complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
   static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

   for (iter=1;iter<=MAXIT;iter++)
   {
      b = Coeff[Degree];
      err = abs(b);
      d=f = double_complex(0.0,0.0);
      abx = abs(*X);

      // Efficient computation of the polynomial and its first two derivatives
      for (j=Degree-1;j>=0;j--)
      {
	 f = (*X)*f + d;
	 d = (*X)*d + b;
	 b = (*X)*b + Coeff[j];
	 err = abs(b) + abx*err;
      }

      // Estimate of roundoff error in evaluating polynomial
      err *= EPSS;

      // we are on the root
      if (abs(b) <= err) return iter;
      // the generic case: use Laguerre's formula

      g = d/b;
      g2 = g*g;
      h = g2 - (2.0*(f/b));
      sq = sqrt((Degree-1)*(Degree*h - g2));
      gp = g + sq;
      gm = g - sq;
      abp = abs(gp);
      abm = abs(gm);
      if (abp < abm) gp = gm;
      dx = (( FMAX(abp,abm) > 0.0 ? double_complex(Degree,0.0)/gp
	      : (1.0 + abx)*(double_complex(cos(double(iter)),sin(double(iter))))));
      x1 = *X - dx;

      if (X->real() == x1.real() && X->imag() == x1.imag()) return iter;
      if (iter % MT)
	 (*X) = x1;
      else
	 (*X) = (*X) - ((double(iter)/MT)*dx);
   }

   cerr << "too many iteration in Laguerre" << endl;
   return iter;
}
      
void PolyMult(double_complex A[],int DegA,double_complex B[],int DegB,
	      double_complex Result[])
{
   int DegC = DegA + DegB;

   for (int i=0;i<=DegC;i++) Result[i] = double_complex(0.0,0.0);
   for (int i=0;i<=DegA;i++)
      for (int j=0;j<=DegB;j++)
      {
	 Result[i+j] += A[i]*B[j];
      }
}