#include "KnownLattices.h"
#include "SquarePressTempPairPotLat.h"
#include "TrianglePressTempPairPotLat.h"
#include "NiTiPressTempPairPotLat.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile)
{
   FILE *pipe;
   char command[LINELENGTH];

   enum lattices {SquarePressTempPairPot,TrianglePressTempPairPot,
		  NiTiPressTempPairPot};
   lattices lattice;

   char lat[]="^MainLatticeType";
   SetPerlCommand(command,datafile,lat);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(lat);
   if ((!strcmp("SquarePressTempPairPotLat",command))
       || (!strcmp("squarepresstemppairpotlat",command)))
      lattice = SquarePressTempPairPot;
   else if ((!strcmp("TrianglePressTempPairPotLat",command))
       || (!strcmp("trianglepresstemppairpotlat",command)))
      lattice = TrianglePressTempPairPot;
   else if ((!strcmp("NiTiPressTempPairPotLat",command))
       || (!strcmp("nitipresstemppairpotlat",command)))
      lattice = NiTiPressTempPairPot;
   else
   {
      cerr << "Unknown Lattice Type : " << command << endl;
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
