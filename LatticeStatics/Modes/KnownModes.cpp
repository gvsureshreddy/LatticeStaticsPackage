#include "KnownModes.h"

LatticeMode* const InitializeMode(Lattice* const Lat,PerlInput const& Input)
{
   const char *Mode = Input.getString("Mode","Type");

   if (!strcmp("MultiMode",Mode))
   {
      return new MultiMode(Lat,Input);
   }
   else
   {
      cerr << "Unknown Mode Type" << "\n";
      exit(-1);
   }
   
   return 0;
}

