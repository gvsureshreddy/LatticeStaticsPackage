#include "UtilityFunctions.h"

void SetPerlCommand(char *string,char *datafile,char *tag)
{
   char format[]=
      {"perl -e \"\\$found=1; while (<>){"\
       "if (/%s/) {\\$found=0; split('=',\\$_); print eval(\\$_[1]);}} "\
       "exit \\$found; \" %s"};
      sprintf(string,format,tag,datafile);
}

FILE *OpenPipe(char *command,char *mode)
{
   FILE *pipe;

   pipe=popen(command,mode);
   if (!pipe)
   {
      cerr << "popen failed! -- " << endl;
      exit(-1);
   }

   return pipe;
}

void Errfun(char *string)
{
   cerr << "Error -- Unable to find : "
	<< string << endl;
   exit(-1);
}
