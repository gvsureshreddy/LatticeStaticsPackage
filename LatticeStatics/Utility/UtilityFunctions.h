#ifndef __UtilityFunctions
#define __UtilityFunctions

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>

#define LINELENGTH 257

// Utility function definitions
void SetPerlCommand(char *string,char *datafile,char *tag);
void Errfun(char *string);
FILE *OpenPipe(char *command,char *mode);

#endif
