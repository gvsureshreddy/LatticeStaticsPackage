#include "KnownLattices.h"
#include "SquarePressTempPairPotLat.h"
#include "TrianglePressTempPairPotLat.h"
#include "NiTiPressTempPairPotLat.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile)
{
   char tmp[LINELENGTH];

   enum lattices {SquarePressTempPairPot,TrianglePressTempPairPot,
		  NiTiPressTempPairPot};
   lattices lattice;

   GetParameter("^MainLatticeType",datafile,"%s",tmp);
   if ((!strcmp("SquarePressTempPairPotLat",tmp))
       || (!strcmp("squarepresstemppairpotlat",tmp)))
      lattice = SquarePressTempPairPot;
   else if ((!strcmp("TrianglePressTempPairPotLat",tmp))
       || (!strcmp("trianglepresstemppairpotlat",tmp)))
      lattice = TrianglePressTempPairPot;
   else if ((!strcmp("NiTiPressTempPairPotLat",tmp))
       || (!strcmp("nitipresstemppairpotlat",tmp)))
      lattice = NiTiPressTempPairPot;
   else
   {
      cerr << "Unknown Lattice Type : " << tmp << endl;
      exit(-1);
   }

   switch (lattice)
   {
      case SquarePressTempPairPot:
      {
	 return new SquarePressTempPairPotLat(datafile);
      }
      break;
      case TrianglePressTempPairPot:
      {
	 return new TrianglePressTempPairPotLat(datafile);
      }
      break;
      case NiTiPressTempPairPot:
      {
	 return new NiTiPressTempPairPotLat(datafile);
      }
      break;
   }

   return NULL;
}
