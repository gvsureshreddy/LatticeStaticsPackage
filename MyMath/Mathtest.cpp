#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "MyMath.h"

int main()
{
   const int m=6;
   complex<double> Coeff[m+1],tst;

   Coeff[0] = complex<double>(6.0,0.0);
   Coeff[1] = complex<double>(5.0,0.0);
   Coeff[2] = complex<double>(1.0,0.0);
   Coeff[3] = complex<double>(12.0,0.0);
   Coeff[4] = complex<double>(8.0,0.0);
   Coeff[5] = complex<double>(0.5,0.0);
   Coeff[6] = complex<double>(1.5,0.0);

   complex<double> Roots[m+1];

   PolyRootsLaguerre(Coeff,2,Roots,1);

   for (int i=0;i<2;i++)
      cout << setw(20) << Roots[i] << endl;

   tst = -Roots[0];
   for (int i=1;i<2;i++)
      tst *= -Roots[i];

   // this should be Coeff[0]
   cout << tst*Coeff[2] << endl;
   cout << endl;
   
   PolyRootsLaguerre(Coeff,m,Roots,1);

   for (int i=0;i<m;i++)
      cout << setw(20) << Roots[i] << endl;

   tst = -Roots[0];
   for (int i=1;i<m;i++)
      tst *= -Roots[i];

   // this should be Coeff[0]
   cout << tst*Coeff[m] << endl;

   cout << endl << endl;

   complex<double> A[3],B[3],R[5];
   A[0] = complex<double>(-4.0,0.0);
   A[1] = complex<double>(0.0,0.0);
   A[2] = complex<double>(1.0,0.0);
   B[0] = complex<double>(-21.0,0.0);
   B[1] = complex<double>(4.0,0.0);
   B[2] = complex<double>(1.0,0.0);

   PolyMult(A,2,B,2,R);
   for (int i=0;i<=4;i++)
      cout << setw(20) << R[i];
   cout << endl;

   PolyRootsLaguerre(R,4,Roots,1);
   for (int i=0;i<4;i++)
      cout << setw(20) << Roots[i] << endl;

}
