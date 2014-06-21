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
#ifdef USE_KIM
   else if (!strcmp("MultiLatticeKIM", Lat))
   {
      return new MultiLatticeKIM(Input, Echo, Width, Debug);
   }
#endif
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
   else if (!strcmp("FourBarTruss", Lat))
   {
      return new FourBarTruss(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternal", Lat))
   {
      return new DFTExternal(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternalOld", Lat))
   {
      return new DFTExternalOld(Input, Echo, Width);
   }
#ifdef USE_QC
   else if (!strcmp("QC", Lat))
   {
      return new QC(Input, Echo, Width);
   }
#endif
#ifdef USE_FEAP
   else if (!strcmp("FEAP", Lat))
   {
      return new FEAP(Input, Echo, Width);
   }
#endif
   else
   {
      cerr << "Unknown Lattice Type " << "\n";
      exit(-1);
   }
   Input.EndofInputSection();

   return 0;
}
