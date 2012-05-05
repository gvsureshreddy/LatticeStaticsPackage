#ifndef RSE__KnownLattices
#define RSE__KnownLattices

#include "PerlInput.h"
#include "MultiLatticeTPP.h"
#include "MultiChainTPP.h"
#include "MultiChainTTPP.h"
#include "SCLDQMultiChainTPP.h"
#include "SCLDCMultiChainTPP.h"
#include "QHQMultiChainTPP.h"
#include "QHCMultiChainTPP.h"
#include "TwoBarTruss.h"
#include "TwoBarTrussExternal.h"
#include "FourBarTruss.h"
#include "DFTExternal.h"
#include "DFTExternalOld.h"
#include "QC.h"
#include "FEAP.h"

Lattice* const InitializeLattice(PerlInput& Input, int const& Echo, int const& Width = 20,
                                 int const& Debug = 0);

#endif
