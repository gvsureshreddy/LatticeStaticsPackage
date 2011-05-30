#include <cstdlib>
#include <iomanip>
#include "CMatrix.h"


int main(int argc, char* argv[])
{
   if (argc != 2)
   {
      cerr << "Usage: " << argv[0] << " size_int" << "\n";
      exit(-1);
   }

   int const size = atoi(argv[1]);
   CMatrix Z(size, size);

   cin >> Z;

   CMatrix EigVec(size, size);

   cout << setw(20)
        << HermiteEigVal(Z, &EigVec)
        << setw(20)
        << EigVec;

   return 0;
}
