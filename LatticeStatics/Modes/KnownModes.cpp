#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile)
{
   FILE *pipe;
   char command[LINELENGTH];

   enum modes {UniDefTemp2DExp,UniDefTemp3DExp,UniDefTemp3Dnacl};
   modes mode;

   char mod[]="^MainModeType";
   SetPerlCommand(command,datafile,mod);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(mod);
   if ((!strcmp("UniDefTemp2DExpand",command))
       || (!strcmp("unideftemp2dexpand",command)))
      mode = UniDefTemp2DExp;
   else if ((!strcmp("UniDefTemp3DExpand",command))
            || (!strcmp("unideftemp3dexpand",command)))
      mode = UniDefTemp3DExp;
   else if ((!strcmp("UniDefTemp3DNaCl",command))
            || (!strcmp("unideftemp3dnacl",command)))
      mode = UniDefTemp3Dnacl;
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
      case UniDefTemp3DExp:
      {
         return new UniDefTemp3DExpand(Lat);
      }
      break;
      case UniDefTemp3Dnacl:
      {
	 return new UniDefTemp3DNaCl(Lat);
      }
   }

   return NULL;
}

