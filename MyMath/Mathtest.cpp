#include <iostream.h>
#include <stdlib.h>
#include <iomanip.h>

#include "Math.h"

int main()
{
   const int m=6;
   double_complex Coeff[m+1],tst;

   Coeff[0] = double_complex(6.0,0.0);
   Coeff[1] = double_complex(5.0,0.0);
   Coeff[2] = double_complex(1.0,0.0);
   Coeff[3] = double_complex(12.0,0.0);
   Coeff[4] = double_complex(8.0,0.0);
   Coeff[5] = double_complex(0.5,0.0);
   Coeff[6] = double_complex(1.5,0.0);

   double_complex Roots[m+1];

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

   double_complex A[3],B[3],R[5];
   A[0] = double_complex(-4.0,0.0);
   A[1] = double_complex(0.0,0.0);
   A[2] = double_complex(1.0,0.0);
   B[0] = double_complex(-21.0,0.0);
   B[1] = double_complex(4.0,0.0);
   B[2] = double_complex(1.0,0.0);

   PolyMult(A,2,B,2,R);
   for (int i=0;i<=4;i++)
      cout << setw(20) << R[i];
   cout << endl;

   PolyRootsLaguerre(R,4,Roots,1);
   for (int i=0;i<4;i++)
      cout << setw(20) << Roots[i] << endl;

}
