#include "KnownRestrictions.h"
#include <cstdlib>
#include <cstring>

Restriction* const InitializeRestriction(Lattice* const Lat, PerlInput const& Input)
{
   const char* Restrict = Input.getString("Restriction", "Type");

   if (!strcmp("RestrictToTranslatedSubSpace", Restrict))
   {
      return new RestrictToTranslatedSubSpace(Lat, Input);
   }
   else if (!strcmp("NoRestriction", Restrict))
   {
      return new NoRestriction(Lat, Input);
   }
   else
   {
      cerr << "Unknown Restriction Type" << "\n";
      exit(-1);
   }
   Input.EndofInputSection();

   return 0;
}
