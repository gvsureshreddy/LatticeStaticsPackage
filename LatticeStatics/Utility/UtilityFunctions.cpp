#include "UtilityFunctions.h"

void GetParameter(const char *tag,const char *datafile,const char *scanffmt,
		  void *parameter)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,tag);
   pipe = OpenPipe(command,"r");
   fscanf(pipe,scanffmt,parameter);
   if (pclose(pipe)) Errfun(tag);
}

int GetStringParameter(const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb)
{
   int i;
   char strng[LINELENGTH];
   GetParameter(tag,datafile,"%s",strng);
   
   for (i=numb-1;i>=0;i--)
   {
      if (!strcasecmp(strng,choices[i]))
      {
	 return i;
      }
   }

   return i;
}

void SetPerlCommand(char *string,const char *datafile,const char *tag)
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

FILE *OpenPipe(const char *command,const char *mode)
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

void Errfun(const char *string)
{
   cerr << "Error -- Unable to find : "
	<< string << endl;
   exit(-1);
}
