#include "KnownModes.h"

LatticeMode *InitializeMode(Lattice *Lat,PerlInput &Input)
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
   
   return NULL;
}

