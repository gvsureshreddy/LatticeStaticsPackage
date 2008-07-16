#include "KnownRestrictions.h"

Restriction* const InitializeRestriction(Lattice* const Lat,PerlInput const& Input)
{
   const char *Restrict = Input.getString("Restriction","Type");

   if (!strcmp("RestrictToSubSpaceOld",Restrict))
   {
      return new RestrictToSubSpaceOld(Lat,Input);
   }
   else
   {
      cerr << "Unknown Restriction Type" << "\n";
      exit(-1);
   }
   Input.EndofInputSection();
   
   return 0;
}

