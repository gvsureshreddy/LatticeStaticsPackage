#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,const char *datafile,const char *prefix)
{
   const int NoModes = KNOWNMODES;
   const char *Modes[]={"MultiMode"};

   switch (GetStringParameter(prefix,"MainModeType",datafile,Modes,NoModes))
   {
      case 0:
      {
	 return new MultiMode(Lat,datafile,prefix);
      }
      break;
      case -1:
      {
	 cerr << "Unknown Mode Type" << endl;
	 exit(-1);
      }
      break;
   }

   return NULL;
}

