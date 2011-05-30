#ifndef RSE__UtilityFunctions
#define RSE__UtilityFunctions

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

int FullScanRank1Convex3D(CBKinematics const* const CBK, Matrix const& K, double const& dx);
int FullScanRank1Convex2D(Matrix const& K, double const& dx);
int Rank1Convex3D(CBKinematics const* const CBK, Matrix const& K, double const& dx);
int Rank1Convex2D(Matrix const& K, double const& dx);

Matrix TranslationProjection1D(int const& NoAtoms);
Matrix TranslationProjection3D(int const& Fsize, int const& NoAtoms);

#endif
