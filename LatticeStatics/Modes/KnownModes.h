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
#include "OrthoShuffle.h"
#include "RhomboShuffle.h"
#include "FullShuffle.h"

#define KNOWNMODES 15

LatticeMode *InitializeMode(Lattice *Lat,char *datafile);
