#ifndef __UtilityFunctions
#define __UtilityFunctions

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <Matrix.h>
#include <Vector.h>
#include <MyMath.h>
#include <SparseMatrix.h>

#include "CBKinematics.h"

using namespace std;

// Utility function definitions
char kbhitWait();
int EnterDebugMode();

unsigned FullScanRank1Convex3D(CBKinematics *CBK,Matrix K,double dx);
unsigned FullScanRank1Convex2D(Matrix K,double dx);
unsigned Rank1Convex3D(CBKinematics *CBK,Matrix K,double dx);
unsigned Rank1Convex2D(Matrix K,double dx);

Matrix TranslationProjection1D(int NoAtoms);
Matrix TranslationProjection3D(int Fsize,int NoAtoms);

#endif
