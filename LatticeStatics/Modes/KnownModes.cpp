#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile)
{
   const int NoModes = 5;
   const char *Modes[]={"UniDefTemp2DExpand",
		        "UniDefTemp3DExpand",
		        "UniDefTemp3DNaCl",
			"UniDefTemp3DRhombo",
			"UniDefTemp3DOrtho"};

   switch (GetStringParameter("^MainModeType",datafile,Modes,NoModes))
   {
      case 0:
      {
	 return new UniDefTemp2DExpand(Lat);
      }
      break;
      case 1:
      {
         return new UniDefTemp3DExpand(Lat);
      }
      break;
      case 2:
      {
	 return new UniDefTemp3DNaCl(Lat);
      }
      break;
      case 3:
      {
	 return new UniDefTemp3DRhombo(Lat);
      }
      break;
      case 4:
      {
	 return new UniDefTemp3DOrtho(Lat);
      }
      break;
      case -1:
      {
	 cerr << "Unknown Mode Type" << endl;
	 exit(-1);
      }
   }

   return NULL;
}

