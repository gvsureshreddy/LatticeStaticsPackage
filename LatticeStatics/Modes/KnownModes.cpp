#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile)
{
   FILE *pipe;
   char command[LINELENGTH];

   enum modes {UniDefTemp2DExp};
   modes mode;

   char mod[]="^MainModeType";
   SetPerlCommand(command,datafile,mod);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(mod);
   if ((!strcmp("UniDefTemp2DExpand",command))
       || (!strcmp("unideftemp2dexpand",command)))
      mode = UniDefTemp2DExp;
   else
   {
      cerr << "Unknown Mode Type : " << command << endl;
      exit(-1);
   }

   switch (mode)
   {
      case UniDefTemp2DExp:
      {
	 return new UniDefTemp2DExpand(Lat);
      }
      break;
   }
}

