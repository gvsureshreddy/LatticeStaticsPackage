#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile)
{
   char tmp[LINELENGTH];

   enum modes {UniDefTemp2DExp,UniDefTemp3DExp,UniDefTemp3Dnacl};
   modes mode;

   GetParameter("^MainModeType",datafile,"%s",tmp);
   if ((!strcmp("UniDefTemp2DExpand",tmp))
       || (!strcmp("unideftemp2dexpand",tmp)))
      mode = UniDefTemp2DExp;
   else if ((!strcmp("UniDefTemp3DExpand",tmp))
            || (!strcmp("unideftemp3dexpand",tmp)))
      mode = UniDefTemp3DExp;
   else if ((!strcmp("UniDefTemp3DNaCl",tmp))
            || (!strcmp("unideftemp3dnacl",tmp)))
      mode = UniDefTemp3Dnacl;
   else
   {
      cerr << "Unknown Mode Type : " << tmp << endl;
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

