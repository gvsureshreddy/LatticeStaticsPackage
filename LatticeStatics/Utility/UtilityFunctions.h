#ifndef __UtilityFunctions
#define __UtilityFunctions

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINELENGTH 600

// Utility function definitions
void GetParameter(const char *tag,const char *datafile,const char *scanffmt,
		  void *parameter);
int GetStringParameter(const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb);
void SetPerlCommand(char *string,const char *datafile,const char *tag);
void Errfun(const char *string);
FILE *OpenPipe(const char *command,const char *mode);

#endif
