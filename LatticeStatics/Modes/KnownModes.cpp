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
			"Expand6",
			"ExpandShuffle",
			"Shuffle3NaCl",
			"Shuffle3Rhombo",
			"Shuffle3Ortho",
			"Shuffle3Mono3P",
			"Shuffle3Mono3M",
			"FullShuffle1",
			"FullShuffle2",
			"FullShuffle3",
			"Expand15",
			"Full15",
			"Ortho15Shuff",
			"Expand9",
			"NaCl9"};

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
	 return new Expand6(Lat);
      }
      break;
      case 12:
      {
	 return new ExpandShuffle(Lat);
      }
      break;
      case 13:
      {
	 return new Shuffle3NaCl(Lat);
      }
      break;
      case 14:
      {
	 return new Shuffle3Rhombo(Lat);
      }
      break;
      case 15:
      {
	 return new Shuffle3Ortho(Lat);
      }
      break;
      case 16:
      {
	 return new Shuffle3Mono3P(Lat);
      }
      break;
      case 17:
      {
	 return new Shuffle3Mono3M(Lat);
      }
      break;
      case 18:
      {
	 return new FullShuffle1(Lat);
      }
      break;
      case 19:
      {
	 return new FullShuffle2(Lat);
      }
      break;
      case 20:
      {
	 return new FullShuffle3(Lat);
      }
      break;
      case 21:
      {
	 return new Expand15(Lat);
      }
      break;
      case 22:
      {
	 return new Full15(Lat);
      }
      break;
      case 23:
      {
	 return new Ortho15Shuff(Lat);
      }
      break;
      case 24:
      {
	 return new Expand9(Lat);
      }
      break;
      case 25:
      {
	 return new NaCl9(Lat);
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

