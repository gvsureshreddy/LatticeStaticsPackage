#include <stdlib.h>
#include <iomanip.h>
#include "Matrix.h"


int main(int argc, char *argv[])
{
   if (argc != 2)
   {
      cerr << "Usage: " << argv[0] << " size_int" << endl;
      exit(-1);
   }

   int size=atoi(argv[1]);
   Matrix Z(size,size);

   cin >> Z;

   Matrix EigVec(size,size);

   cout << setw(20)
	<< SymEigVal(Z,&EigVec)
	<< setw(20)
	<< EigVec;

   return 0;
}
