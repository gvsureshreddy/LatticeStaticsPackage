#ifndef __UtilityFunctions
#define __UtilityFunctions

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Matrix.h>
#include <MyMath.h>

#define LINELENGTH 600

// Utility function definitions
void GetParameter(const char *tag,const char *datafile,const char *scanffmt,
		  void *parameter);
int GetStringParameter(const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb);
void SetPerlCommand(char *string,const char *datafile,const char *tag);
void Errfun(const char *string);
FILE *OpenPipe(const char *command,const char *mode);

unsigned FullScanRank1Convex3D(Matrix K,double dx);
unsigned FullScanRank1Convex2D(Matrix K,double dx);
unsigned Rank1Convex3D(Matrix K,double dx);
unsigned Rank1Convex2D(Matrix K,double dx);

static const double UtilityALT[3][3][3] = {0.0, 0.0, 0.0,
					   0.0, 0.0, 1.0,
					   0.0, -1.0, 0.0,
					   0.0, 0.0, -1.0,
					   0.0, 0.0, 0.0,
					   1.0, 0.0, 0.0,
					   0.0, 1.0, 0.0,
					   -1.0, 0.0, 0.0,
					   0.0, 0.0, 0.0};

inline double Alt(int i,int j,int k) {return UtilityALT[i][j][k];}


char *builddate();

#endif
