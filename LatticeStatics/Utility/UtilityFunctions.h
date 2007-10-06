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

using namespace std;

#define LINELENGTH 600

// Utility function definitions
char kbhitWait();
int EnterDebugMode();

extern char UTILITYechocommand[LINELENGTH];
extern char *UTILITYechocmd;
int GetParameter(const char *prefix,const char *tag,const char *datafile,
		 const char *scanffmt,void *parameter,int DispErr=1);
int GetVectorParameter(const char *prefix,const char *tag,
		       const char *datafile,Vector *V,int DispErr=1);
int GetIntVectorParameter(const char *prefix,const char *tag,const char *datafile,
			  int N,int *Vec,int DispErr=1);
int GetMatrixParameter(const char *prefix,const char *tag,const char *datafile,Matrix *M,
		       int DispErr=1);
int GetStringParameter(const char *prefix,const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb,int DispErr=1);
void SetPerlCommand(char *string,const char *datafile,const char *prefix,const char *tag);
void Errfun(const char *string);
FILE *OpenPipe(const char *command,const char *mode);

unsigned FullScanRank1Convex3D(Matrix K,double dx);
unsigned FullScanRank1Convex2D(Matrix K,double dx);
unsigned Rank1Convex3D(Matrix K,double dx);
unsigned Rank1Convex2D(Matrix K,double dx);

Matrix TranslationProjection1D(int NoAtoms);
Matrix TranslationProjection3D(int Fsize,int NoAtoms);

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
