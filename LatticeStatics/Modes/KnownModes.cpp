#include "KnownModes.h"

#include "UtilityFunctions.h"

LatticeMode *InitializeMode(Lattice *Lat,char *datafile)
{
   const int NoModes = KNOWNMODES;
   const char *Modes[]={"UniDefTemp2DExpand",
		        "UniDefTemp3DExpand",
		        "UniDefTemp3DNaCl",
			"UniDefTemp3DTetrag",
			"UniDefTemp3DRhombo",
			"UniDefTemp3DOrtho",
			"UniDefTemp3D3PMono",
			"UniDefTemp3D3MMono",
			"UniDefTemp3D1PMono",
			"UniDefTemp3D1MMono",
                        "UniDefTemp3DFull",
			"ExpandShuffle",
			"NaClShuffle",
			"RhomboShuffle",
			"OrthoShuffle",
			"Mono3PShuffle",
			"Mono3MShuffle",
			"FullShuffle"};

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
	 return new UniDefTemp3DTetrag(Lat);
      }
      case 4:
      {
	 return new UniDefTemp3DRhombo(Lat);
      }
      break;
      case 5:
      {
	 return new UniDefTemp3DOrtho(Lat);
      }
      break;
      case 6:
      {
	 return new UniDefTemp3D3PMono(Lat);
      }
      break;
      case 7:
      {
	 return new UniDefTemp3D3MMono(Lat);
      }
      break;
      case 8:
      {
	 return new UniDefTemp3D1PMono(Lat);
      }
      break;
      case 9:
      {
	 return new UniDefTemp3D1MMono(Lat);
      }
      break;
      case 10:
      {
         return new UniDefTemp3DFull(Lat);
      }
      break;
      case 11:
      {
	 return new ExpandShuffle(Lat);
      }
      case 12:
      {
	 return new NaClShuffle(Lat);
      }      
      case 13:
      {
	 return new RhomboShuffle(Lat);
      }
      case 14:
      {
	 return new OrthoShuffle(Lat);
      }
      case 15:
      {
	 return new Mono3PShuffle(Lat);
      }
      case 16:
      {
	 return new Mono3MShuffle(Lat);
      }
      case 17:
      {
	 return new FullShuffle(Lat);
      }
      case -1:
      {
	 cerr << "Unknown Mode Type" << endl;
	 exit(-1);
      }
   }

   return NULL;
}

