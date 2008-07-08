#ifndef __KnownLattices
#define __KnownLattices

#include "PerlInput.h"
#include "MultiLatticeTPP.h"
#include "MultiChainTPP.h"
#include "MultiChainTTPP.h"
#include "TwoBarTruss.h"

Lattice* const InitializeLattice(PerlInput& Input,int const& Echo,int const& Width=20,
                                 int const& Debug=0);

#endif
