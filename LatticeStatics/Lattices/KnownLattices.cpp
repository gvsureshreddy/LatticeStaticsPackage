#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile,const char *prefix)
{
   const int NoLats = 3;
   const char *Lattices[]={"SquarePressTempPairPotLat",
		           "TrianglePressTempPairPotLat",
			   "MultiLatticeTPP"};
   
   switch (GetStringParameter(prefix,"MainLatticeType",datafile,Lattices,NoLats))
   {
      case 0:
      {
	 return new SquarePressTempPairPotLat(datafile,prefix);
      }
      break;
      case 1:
      {
	 return new TrianglePressTempPairPotLat(datafile,prefix);
      }
      break;
      case 2:
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
