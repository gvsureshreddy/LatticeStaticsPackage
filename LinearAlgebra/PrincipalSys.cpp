#include <cstdlib>
#include <iomanip>
#include "Matrix.h"
#include "Vector.h"


int main(int argc, char* argv[])
{
   if (argc != 7)
   {
      cerr << "Usage: " << argv[0] << " A11 A22 A33 A12 A13 A23" << "\n";
      exit(-1);
   }

   int const size = 3;
   Matrix A(size, size);
   A[0][0] = atof(argv[1]);
   A[1][1] = atof(argv[2]);
   A[2][2] = atof(argv[3]);
   A[0][1] = A[1][0] = atof(argv[4]);
   A[0][2] = A[2][0] = atof(argv[5]);
   A[1][2] = A[2][1] = atof(argv[6]);

   Matrix EigVec(size, size), EigVals = SymEigVal(A, &EigVec);

   double tmp;
   if (EigVals[0][1] > EigVals[0][0])
   {
      tmp = EigVals[0][0];
      EigVals[0][0] = EigVals[0][1];
      EigVals[0][1] = tmp;
      for (int i = 0; i < size; ++i)
      {
         tmp = EigVec[i][0];
         EigVec[i][0] = EigVec[i][1];
         EigVec[i][1] = tmp;
      }
   }
   if (EigVals[0][2] > EigVals[0][0])
   {
      tmp = EigVals[0][0];
      EigVals[0][0] = EigVals[0][2];
      EigVals[0][2] = tmp;
      for (int i = 0; i < size; ++i)
      {
         tmp = EigVec[i][0];
         EigVec[i][0] = EigVec[i][2];
         EigVec[i][2] = tmp;
      }
   }
   if (EigVals[0][2] > EigVals[0][1])
   {
      tmp = EigVals[0][1];
      EigVals[0][1] = EigVals[0][2];
      EigVals[0][2] = tmp;
      for (int i = 0; i < size; ++i)
      {
         tmp = EigVec[i][1];
         EigVec[i][1] = EigVec[i][2];
         EigVec[i][2] = tmp;
      }
   }


   Vector X(size, 0.0),
   Z(size, 0.0),
   Xp(size, 0.0),
   Yp(size, 0.0),
   Zp(size, 0.0),
   LineofNodes(size, 0.0);

   X[0] = 1.0;
   Z[2] = 1.0;
   for (int i = 0; i < size; ++i)
   {
      Xp[i] = EigVec[i][0];
      Yp[i] = EigVec[i][1];
      Zp[i] = EigVec[i][2];
   }
   LineofNodes = Yp % Zp;

   for (int i = 0; i < size; ++i)
      cout << setw(20) << EigVals[0][i];

   cout << setw(20) << acos(X * LineofNodes);
   cout << setw(20) << acos(Z * Zp);
   cout << setw(20) << acos(Xp * LineofNodes);
   cout << "\n";

   return 0;
}

