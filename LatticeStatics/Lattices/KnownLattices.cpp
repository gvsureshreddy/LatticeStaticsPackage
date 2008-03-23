#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice *InitializeLattice(PerlInput &Input,int Echo,int Width,int Debug)
{
   const char *Lat = Input.getString("Lattice","Type");

   if (!strcmp("MultiLatticeTPP",Lat))
   {
      return new MultiLatticeTPP(Input,Echo,Width,Debug);
   }
   else if (!strcmp("MultiChainTPP",Lat))
   {
      return new MultiChainTPP(Input,Echo,Width,Debug);
   }
   else if (!strcmp("MultiChainTTPP",Lat))
   {
      return new MultiChainTTPP(Input,Echo,Width,Debug);
   }
   else
   {
      cerr << "Unknown Lattice Type " << "\n";
      exit(-1);
   }
   
   return NULL;
}
