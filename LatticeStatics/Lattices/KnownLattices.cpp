#include "KnownLattices.h"
#include "SquarePressTempPairPotLat.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile)
{
   FILE *pipe;
   char command[LINELENGTH];

   enum lattices {SquarePressTempPairPot};
   lattices lattice;

   char lat[]="^MainLatticeType";
   SetPerlCommand(command,datafile,lat);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(lat);
   if ((!strcmp("SquarePressTempPairPotLat",command))
       || (!strcmp("squarepresstemppairpotlat",command)))
      lattice = SquarePressTempPairPot;
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
   }
}
