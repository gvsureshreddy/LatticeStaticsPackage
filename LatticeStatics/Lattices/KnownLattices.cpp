#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile)
{
   const int NoLats = 6;
   const char *Lattices[]={"SquarePressTempPairPotLat",
		           "TrianglePressTempPairPotLat",
		           "NiTiPressTempPairPotLat",
		           "NiTiShuffleTPPLat",
			   "NiTiRadii",
			   "NiTi15TPPLat"};
   
   switch (GetStringParameter("^MainLatticeType",datafile,Lattices,NoLats))
   {
      case 0:
      {
	 return new SquarePressTempPairPotLat(datafile);
      }
      break;
      case 1:
      {
	 return new TrianglePressTempPairPotLat(datafile);
      }
      break;
      case 2:
      {
	 return new NiTiPressTempPairPotLat(datafile);
      }
      break;
      case 3:
      {
	 return new NiTiShuffleTPPLat(datafile);
      }
      break;
      case 4:
      {
	 return new NiTiRadii(datafile);
      }
      break;
      case 5:
      {
	 return new NiTi15TPPLat(datafile);
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
