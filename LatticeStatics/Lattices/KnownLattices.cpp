#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile)
{
   const int NoLats = 9;
   const char *Lattices[]={"SquarePressTempPairPotLat",
		           "TrianglePressTempPairPotLat",
			   "NiTiPressTempPairPotLat",
		           "NiTi6TPPLat",
		           "NiTiShuffle1TPPLat",
		           "NiTiShuffle2TPPLat",
		           "NiTiShuffle3TPPLat",
			   "NiTi15TPPLat",
			   "NiTi9TPPLat"};
   
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
	 return new NiTi6TPPLat(datafile);
      }
      break;
      case 4:
      {
	 return new NiTiShuffle1TPPLat(datafile);
      }
      break;
      case 5:
      {
	 return new NiTiShuffle2TPPLat(datafile);
      }
      break;
      case 6:
      {
	 return new NiTiShuffle3TPPLat(datafile);
      }
      break;
      case 7:
      {
	 return new NiTi15TPPLat(datafile);
      }
      break;
      case 8:
      {
	 return new NiTi9TPPLat(datafile);
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
