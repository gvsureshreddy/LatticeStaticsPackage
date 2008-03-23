#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "MyMath.h"
#include "MyComplexDouble.h"

int main()
{
   const int m=6;
   MyComplexDouble Coeff[m+1],tst;

   Coeff[0] = MyComplexDouble(6.0,0.0);
   Coeff[1] = MyComplexDouble(5.0,0.0);
   Coeff[2] = MyComplexDouble(1.0,0.0);
   Coeff[3] = MyComplexDouble(12.0,0.0);
   Coeff[4] = MyComplexDouble(8.0,0.0);
   Coeff[5] = MyComplexDouble(0.5,0.0);
   Coeff[6] = MyComplexDouble(1.5,0.0);

   for (int i=0;i<m+1;++i)
   {
      cout << "Coeff[" << i << "]=" << setw(20) << Coeff[i] << "\n";
   }
   
   MyComplexDouble Roots[m+1];

   PolyRootsLaguerre(Coeff,2,Roots,1);

   for (int i=0;i<2;i++)
      cout << setw(20) << Roots[i] << "\n";

   tst = -Roots[0];
   for (int i=1;i<2;i++)
      tst *= -Roots[i];

   // this should be Coeff[0]
   cout << tst*Coeff[2] << "\n";
   cout << "\n";
   
   PolyRootsLaguerre(Coeff,m,Roots,1);
   for (int i=0;i<m;i++)
      cout << setw(20) << Roots[i] << "\n";

   tst = -Roots[0];
   for (int i=1;i<m;i++)
      tst *= -Roots[i];

   // this should be Coeff[0]
   cout << tst*Coeff[m] << "\n";

   cout << "\n" << "\n";

   MyComplexDouble A[3],B[3],R[5];
   A[0] = MyComplexDouble(-4.0,0.0);
   A[1] = MyComplexDouble(0.0,0.0);
   A[2] = MyComplexDouble(1.0,0.0);
   B[0] = MyComplexDouble(-21.0,0.0);
   B[1] = MyComplexDouble(4.0,0.0);
   B[2] = MyComplexDouble(1.0,0.0);

   PolyMult(A,2,B,2,R);
   for (int i=0;i<=4;i++)
      cout << setw(20) << R[i];
   cout << "\n";

   PolyRootsLaguerre(R,4,Roots,1);
   for (int i=0;i<4;i++)
      cout << setw(20) << Roots[i] << "\n";

}
