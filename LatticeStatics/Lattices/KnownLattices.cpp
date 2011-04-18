#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice* const InitializeLattice(PerlInput& Input, int const& Echo, int const& Width,
                                 int const& Debug)
{
   char const* const Lat = Input.getString("Lattice", "Type");

   if (!strcmp("MultiLatticeTPP", Lat))
   {
      return new MultiLatticeTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("MultiChainTPP", Lat))
   {
      return new MultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("MultiChainTTPP", Lat))
   {
      return new MultiChainTTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("SCLDQMultiChainTPP", Lat))
   {
      return new SCLDQMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("SCLDCMultiChainTPP", Lat))
   {
      return new SCLDCMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("QHQMultiChainTPP", Lat))
   {
      return new QHQMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("QHCMultiChainTPP", Lat))
   {
      return new QHCMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("TwoBarTruss", Lat))
   {
      return new TwoBarTruss(Input, Echo, Width);
   }
   else if (!strcmp("TwoBarTrussExternal", Lat))
   {
      return new TwoBarTrussExternal(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternal", Lat))
   {
      return new DFTExternal(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternalOld", Lat))
   {
      return new DFTExternalOld(Input, Echo, Width);
   }
   else if (!strcmp("QC", Lat))
   {
      return new QC(Input, Echo, Width);
   }
   else
   {
      cerr << "Unknown Lattice Type " << "\n";
      exit(-1);
   }
   Input.EndofInputSection();

   return 0;
}

