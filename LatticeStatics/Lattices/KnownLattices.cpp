#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile,const char *prefix)
{
   const int NoLats = 1;
   const char *Lattices[]={"MultiLatticeTPP"};
   
   switch (GetStringParameter(prefix,"MainLatticeType",datafile,Lattices,NoLats))
   {
      case 0:
      {
	 return new MultiLatticeTPP(datafile,prefix);
      }
      break;
      case -1:
	 {
	    cerr << "Unknown Lattice Type " << endl;
	    exit(-1);
	 }
   }

   return NULL;
}
