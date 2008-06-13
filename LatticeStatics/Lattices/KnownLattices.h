#ifndef __KnownLattices
#define __KnownLattices

#include "PerlInput.h"
#include "MultiLatticeTPP.h"
#include "MultiChainTPP.h"
#include "MultiChainTTPP.h"
#include "TwoBarTruss.h"

Lattice *InitializeLattice(PerlInput &Input,int Echo,int Width=20,int Debug=0);

#endif
