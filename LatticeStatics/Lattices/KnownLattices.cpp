#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(char *datafile,const char *prefix,int Echo,int Width,int Debug)
{
   const int NoLats = 2;
   const char *Lattices[]={"MultiLatticeTPP","MultiChainTPP"};
   
   switch (GetStringParameter(prefix,"MainLatticeType",datafile,Lattices,NoLats))
   {
      case 0:
      {
	 return new MultiLatticeTPP(datafile,prefix,Echo,Width,Debug);
      }
      case 1:
      {
	 return new MultiChainTPP(datafile,prefix,Echo,Width,Debug);
      }
      case -1:
	 {
	    cerr << "Unknown Lattice Type " << endl;
	    exit(-1);
	 }
   }

   return NULL;
}
