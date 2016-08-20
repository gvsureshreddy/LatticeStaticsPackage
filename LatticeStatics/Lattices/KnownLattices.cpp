#include "KnownLattices.h"
#include "UtilityFunctions.h"

Lattice* const InitializeLattice(PerlInput const& Input, int const& Echo,
                                 int const& Width, int const& Debug)
{
  char const* const LatticeType = Input.getString("Lattice", "Type");
  return InitializeLattice(LatticeType, Input, Echo, Width, Debug);
}


Lattice* const InitializeLattice(char const* const LatticeType,
                                 PerlInput const& Input, int const& Echo,
                                 int const& Width, int const& Debug)
{
   if (!strcmp("DisplacementControl", LatticeType))
   {
      return new DisplacementControl(Input, Echo, Width);
   }
   else if (!strcmp("MultiLatticeTPP", LatticeType))
   {
      return new MultiLatticeTPP(Input, Echo, Width, Debug);
   }
#ifdef USE_KIM
   else if (!strcmp("MultiLatticeKIM", LatticeType))
   {
      return new MultiLatticeKIM(Input, Echo, Width);
   }
#endif
   else if (!strcmp("MultiChainTPP", LatticeType))
   {
      return new MultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("MultiChainTTPP", LatticeType))
   {
      return new MultiChainTTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("SCLDQMultiChainTPP", LatticeType))
   {
      return new SCLDQMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("SCLDCMultiChainTPP", LatticeType))
   {
      return new SCLDCMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("QHQMultiChainTPP", LatticeType))
   {
      return new QHQMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("QHCMultiChainTPP", LatticeType))
   {
      return new QHCMultiChainTPP(Input, Echo, Width, Debug);
   }
   else if (!strcmp("TwoBarTruss", LatticeType))
   {
      return new TwoBarTruss(Input, Echo, Width);
   }
   else if (!strcmp("TwoBarTrussExternal", LatticeType))
   {
      return new TwoBarTrussExternal(Input, Echo, Width);
   }
   else if (!strcmp("FourBarTruss", LatticeType))
   {
      return new FourBarTruss(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternal", LatticeType))
   {
      return new DFTExternal(Input, Echo, Width);
   }
   else if (!strcmp("DFTExternalOld", LatticeType))
   {
      return new DFTExternalOld(Input, Echo, Width);
   }
#ifdef USE_QC
   else if (!strcmp("QC", LatticeType))
   {
      return new QC(Input, Echo, Width);
   }
#endif
#ifdef USE_FEAP
   else if (!strcmp("FEAP", LatticeType))
   {
      return new FEAP(Input, Echo, Width);
   }
#endif
#ifdef USE_NEO
   else if (!strcmp("NeoHookean2D", LatticeType))
   {
      return new NeoHookean2D(Input, Echo, Width);
   }
#endif
   else
   {
      cerr << "Unknown Lattice Type " << LatticeType << "\n";
      exit(-1);
   }
   Input.EndofInputSection();

   return 0;
}
