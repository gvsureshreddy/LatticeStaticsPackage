#include "UniDefTemp2DExpand.h"
#include "UniDefTemp3DExpand.h"
#include "UniDefTemp3DNaCl.h"
#include "UniDefTemp3DTetrag.h"
#include "UniDefTemp3DRhombo.h"
#include "UniDefTemp3DOrtho.h"
#include "UniDefTemp3D3PMono.h"
#include "UniDefTemp3D3MMono.h"
#include "UniDefTemp3D1PMono.h"
#include "UniDefTemp3D1MMono.h"
#include "UniDefTemp3DFull.h"
#include "Expand6.h"
#include "ExpandShuffle.h"
#include "Shuffle3NaCl.h"
#include "Shuffle3Rhombo.h"
#include "Shuffle3Ortho.h"
#include "Shuffle3Mono3P.h"
#include "Shuffle3Mono3M.h"
#include "FullShuffle1.h"
#include "FullShuffle2.h"
#include "FullShuffle3.h"
#include "Expand15.h"
#include "Full15.h"
#include "Ortho15Shuff.h"
#include "Expand9.h"

#define KNOWNMODES 25

LatticeMode *InitializeMode(Lattice *Lat,char *datafile);
