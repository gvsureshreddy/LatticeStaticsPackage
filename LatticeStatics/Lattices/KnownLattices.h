#ifndef RSE__KnownLattices
#define RSE__KnownLattices

#include "PerlInput.h"
#ifdef USE_KIM
#include "MultiLatticeKIM.h"
#endif
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
#ifdef USE_QC
#include "QC.h"
#endif
#ifdef USE_FEAP
#include "FEAP.h"
#endif

Lattice* const InitializeLattice(PerlInput& Input, int const& Echo, int const& Width = 20,
                                 int const& Debug = 0);

#endif
