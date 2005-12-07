#include <cstdlib>
#include <iomanip>
#include "Matrix.h"
#include "Vector.h"


int main(int argc, char *argv[])
{
   if (argc != 10)
   {
      cerr << "Usage: " << argv[0] << " A11 A12 A13 A21 A22 A23 A31 A32 A33" << endl;
      exit(-1);
   }

   int size=3;
   Matrix A(size,size);
   for (int i=0;i<size;++i)
      for (int j=0;j<size;++j)
      {
	 A[i][j] = atof(argv[size*i + j + 1]);
      }

   Matrix EigVec(size,size),EigVals=SymEigVal(A,&EigVec);

   Vector X(size,0.0),
      Z(size,0.0),
      Xp(size,0.0),
      Yp(size,0.0),
      Zp(size,0.0),
      LineofNodes(size,0.0);

   X[0] = 1.0;
   Z[2] = 1.0;
   for (int i=0;i<size;++i)
   {
      Xp[i] = EigVec[i][0];
      Yp[i] = EigVec[i][1];
      Zp[i] = EigVec[i][2];
   }
   LineofNodes = Yp % Zp;

   for (int i=0;i<size;++i)
      cout << setw(20) << EigVals[0][i];

   cout << setw(20) << acos(X*LineofNodes);
   cout << setw(20) << acos(Z*Zp);
   cout << setw(20) << acos(Xp*LineofNodes);
   cout << endl;
   
   return 0;
}
