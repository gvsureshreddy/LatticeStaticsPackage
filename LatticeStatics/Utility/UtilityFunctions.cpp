#include "UtilityFunctions.h"

void SetPerlCommand(char *string,char *datafile,char *tag)
{
   char format[]=
      {"perl -e \"\\$R=findref(\\$ARGV[1],\\$ARGV[0]); print \\$R;"\
       "sub findref {my(\\$tag,\\$df) = @_; my(\\$fnd); \\$fnd=1; "\
       "open(R,\\$df); while (<R>) {if (/\\$tag/) {\\$fnd=0; "\
       "\\$_=deref(\\$_,\\$df); split('=',\\$_); return eval(\\$_[1]);}} "\
       "close(R); if (\\$fnd == 1) {exit \\$fnd;}} sub deref "\
       "{my(\\$fld,\\$df)=@_; while (\\$fld =~ m/\\<([^>]+)>/g) "\
       "{\\$v=findref(\\\"^\\$1\\\",\\$df); \\$fld =~ s/<\\$1>/\\$v/} "\
       "return \\$fld;}\" %s %s"};

      sprintf(string,format,datafile,tag);
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
