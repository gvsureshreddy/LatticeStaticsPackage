#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile,const char *prefix)
{
   const int NoModes = KNOWNMODES;
   const char *Modes[]={"UniDefTemp2DExpand",
			"ML_Expand",
			"ML_NaCl",
			"ML_Rhombo",
			"ML_B19"};

   switch (GetStringParameter(prefix,"MainModeType",datafile,Modes,NoModes))
   {
      case 0:
      {
	 return new UniDefTemp2DExpand(Lat);
      }
      break;
      case 1:
      {
	 return new ML_Expand(Lat);
      }
      break;
      case 2:
      {
	 return new ML_NaCl(Lat);
      }
      break;
      case 3:
      {
	 return new ML_Rhombo(Lat);
      }
      break;
      case 4:
      {
	 return new ML_B19(Lat);
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

