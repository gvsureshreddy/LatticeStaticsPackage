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
#include "ExpandShuffle.h"
#include "NaClShuffle.h"
#include "RhomboShuffle.h"
#include "OrthoShuffle.h"
#include "Mono3PShuffle.h"
#include "Mono3MShuffle.h"
#include "FullShuffle.h"
#include "Expand15.h"
#include "Full15.h"

#define KNOWNMODES 20

LatticeMode *InitializeMode(Lattice *Lat,char *datafile);
