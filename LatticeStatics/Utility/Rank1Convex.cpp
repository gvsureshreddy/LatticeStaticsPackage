#include "UtilityFunctions.h"
#include <iomanip.h>

int main(int argc, char *argv[])
{
   if (argc != 2)
   {
      cerr << "Usage: " << argv[0] << " Dx (deg)" << endl;
      exit(1);
   }
   
   double pi = 4.0*atan(1.0);
   double Dx = strtod(argv[1],NULL);
   Matrix A(6,6);

   cin >> A;

   int rk;

//   rk = Rank1Convex3D(A,Dx*pi/180);
   rk = FullScanRank1Convex3D(A,Dx*pi/180);

   cout << "Rank1Convex:" << setw(20) << rk << endl;

   return 0;
}
