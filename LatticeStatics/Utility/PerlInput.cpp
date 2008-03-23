#include <cstdlib>
#include <iostream>
#include "PerlInput.h"

static PerlInterpreter *my_perl = NULL;

void PerlInput::Initialize()
{
   char *args[] = {"","-e ''"};
   if (my_perl != NULL)
   {
      cerr << "PerlInput Error: Can only have one instance of PerlInput!\n";
      exit(-3);
   }
   
   PERL_SYS_INIT3(2,args,NULL);
   my_perl = perl_alloc();
   perl_construct(my_perl);
   PL_exit_flags|=PERL_EXIT_DESTRUCT_END;
   perl_parse(my_perl,NULL,2,args, (char **)NULL);
   perl_run(my_perl);
}

void PerlInput::Readfile(const char *datafile)
{
   char tmp[256];
   SV *retval;
   
   sprintf(tmp,"do '%s';",datafile);
   retval = eval_pv(tmp,TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }
}

void PerlInput::Readfile(const char *datafile,const char *prefix)
{
   char templt[]
      = "$___dmy = \"\"; open(___INPUT,%s) or die \"can't open file for input from perl.\";"
      "while(<___INPUT>){ if ( /^%s/ ) {$_ =~ s/^%s//;"
      "$___dmy = $___dmy . $_;}}; close(___INPUT);";
   char temp[512];
   sprintf(temp,templt,datafile,prefix,prefix);
   
   SV *retval;
   retval = eval_pv(temp,TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }

   retval = eval_pv("eval $___dmy;",TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }
}

PerlInput::PerlInput(char *datafile)
{
   Initialize();
   Readfile(datafile);
}

PerlInput::~PerlInput()
{
   perl_destruct(my_perl);
   perl_free(my_perl);
   PERL_SYS_TERM();
   
   my_perl = NULL;
}

PerlInput::HashStruct PerlInput::getHash(char *HashName)
{
   int Errno = -7;

   HashStruct Hash;
   Hash.Ptr = get_hv(HashName,FALSE);
   if (Hash.Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << HashName << " does not exist.\n";
      exit(Errno);
   }

   strcpy(Hash.Name,HashName);
   return Hash;
}

SV *getScalar(PerlInput::HashStruct Hash,char *ParamName,
              int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);

SV *getScalar(PerlInput::HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   int Errno = -5;
   
   SV **ParamValPtr = hv_fetch(Hash.Ptr,ParamName,strlen(ParamName),FALSE);
   if (ParamValPtr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "} does not exist.\n";
      exit(Errno);
   }
   
   SV *ParamVal = *ParamValPtr;
   if (a <= -1)
   {
      return ParamVal;
   }
   
   // Move into array level 0
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "} does not contain a perl reference (to an array).\n";
      exit(Errno);
   }
   AV *ArrayLevel0 = (AV*) SvRV(ParamVal);
   SV **ScalarLevel1Ptr = av_fetch(ArrayLevel0,a,FALSE);
   if (ScalarLevel1Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "] does not exist.\n";
      exit(Errno);
   }
   
   SV *ScalarLevel1 = *ScalarLevel1Ptr;
   if (b <= -1)
   {
      return ScalarLevel1;
   }
   
   // move into array level 1
   if ((SvTYPE(ScalarLevel1) != SVt_RV) || (SvTYPE(SvRV(ScalarLevel1)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "] does not contain a perl reference (to an array).\n";
      exit(Errno);
   }
   
   AV *ArrayLevel1 = (AV*) SvRV(ScalarLevel1);
   SV **ScalarLevel2Ptr = av_fetch(ArrayLevel1,b,FALSE);
   if (ScalarLevel2Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "] does not exist.\n";
      exit(Errno);
   }
   
   SV *ScalarLevel2 = *ScalarLevel2Ptr;
   if (c <= -1)
   {
      return ScalarLevel2;
   }
   
   // move into array level 2
   cout << "0 at level 2\n";
   if ((SvTYPE(ScalarLevel2) != SVt_RV) || (SvTYPE(SvRV(ScalarLevel2)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "][" << b << "]["
           << c << "] does not contain a perl referece (to an array).\n";
      exit(Errno);
   }
   cout << "1 at level 2\n";
   AV *ArrayLevel2 = (AV*) SvRV(ScalarLevel2);
   SV **ScalarLevel3Ptr = av_fetch(ArrayLevel2,c,FALSE);
   cout << "2 at level 2\n";
   if (ScalarLevel3Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "] does not exist.\n";
      exit(Errno);
   }
   
   SV *ScalarLevel3 = *ScalarLevel3Ptr;
   if (d <= -1)
   {
      return ScalarLevel3;
   }
   
   // move into array level 3
   if ((SvTYPE(ScalarLevel3) != SVt_RV) || (SvTYPE(SvRV(ScalarLevel3)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "][" << b << "]["
           << c << "] does not contain a reference (to an array).\n";
      exit(Errno);
   }
   
   AV *ArrayLevel3 = (AV*) SvRV(ScalarLevel3);
   SV **ScalarLevel4Ptr = av_fetch(ArrayLevel3,d,FALSE);
   if (ScalarLevel4Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "]["
           << d << "] does not exist.\n";
      exit(Errno);
   }
   
   SV *ScalarLevel4 = *ScalarLevel4Ptr;
   if (e <= -1)
   {
      return ScalarLevel4;
   }
   
   // move into array level 4
   if ((SvTYPE(ScalarLevel4) != SVt_RV) || (SvTYPE(SvRV(ScalarLevel4)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "][" << b << "][" << c << "]["
           << d << "] does not contain a reference (to an array).\n";
      exit(Errno);
   }
   
   AV *ArrayLevel4 = (AV*) SvRV(ScalarLevel4);
   SV **ScalarLevel5Ptr = av_fetch(ArrayLevel4,e,FALSE);
   if (ScalarLevel5Ptr == NULL)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "][" << d << "]["
           << e << "] does not exist.\n";
      exit(Errno);
   }
   
   SV *ScalarLevel5 = *ScalarLevel5Ptr;
   if (!SvNOK(ScalarLevel5))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "][" << b << "][" << c << "]["
           << d << "][" << e << "] does not contain a double.\n";
      exit(Errno);
   }
   else
   {
      return ScalarLevel4;
   }
}

PerlInput::HashStruct PerlInput::getHash(HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   int Errno = -6;
   HashStruct NewHash;
   
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d,e);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVHV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      if (e > -1) cerr << "[" << e << "]";
      cerr << " does not contain a hash.\n";
      exit(Errno);
   }

   NewHash.Ptr = (HV*) SvRV(ParamVal);
   if (e > -1)
      sprintf(NewHash.Name,"%s{%s}[%u][%u][%u][%u][%u]",Hash.Name,ParamName,a,b,c,d,e);
   else if (d > -1)
      sprintf(NewHash.Name,"%s{%s}[%u][%u][%u][%u]",Hash.Name,ParamName,a,b,c,d);
   else if (c > -1)
      sprintf(NewHash.Name,"%s{%s}[%u][%u][%u]",Hash.Name,ParamName,a,b,c);
   else if (b > -1)
      sprintf(NewHash.Name,"%s{%s}[%u][%u]",Hash.Name,ParamName,a,b);
   else if (a > -1)
      sprintf(NewHash.Name,"%s{%s}[%u]",Hash.Name,ParamName,a);
   else
      sprintf(NewHash.Name,"%s{%s}",Hash.Name,ParamName);
   return NewHash;
}

int PerlInput::HashOK(char *HashName)
{
   int exists = 1;

   HV *HashPtr = get_hv(HashName,FALSE);
   if (HashPtr == NULL)
   {
      exists = 0;
   }

   return exists;
}

int PerlInput::ParameterOK(HashStruct Hash,char *ParamName)
{
   int exists = 1;
   
   SV **ParamValPtr = hv_fetch(Hash.Ptr,ParamName,strlen(ParamName),FALSE);
   if (ParamValPtr == NULL)
   {
      exists = 0;
   }

   return exists;
}

unsigned PerlInput::getArrayLength(HashStruct Hash,char *ParamName,
                                   int a,int b,int c,int d)
{
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d);
   int Errno=-4;
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " does not contain an array.\n";
      exit(Errno);
   }
   
   return av_len((AV*) SvRV(ParamVal))+1;
}


double PerlInput::getDouble(HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d,e);
   int Errno=-5;
   if (!SvNOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      if (e > -1) cerr << "[" << e << "]";
      cerr << " does not contain a double.\n";
      exit(Errno);
   }
   else
   {
      return SvNV(ParamVal);
   }
}

int PerlInput::getInt(HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d,e);
   int Errno=-5;
   if (!SvIOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      if (e > -1) cerr << "[" << e << "]";
      cerr << " does not contain an int.\n";
      exit(Errno);
   }
   else
   {
      return SvIV(ParamVal);
   }
}

int PerlInput::getUnsigned(HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d,e);
   int Errno=-5;
   if (!SvIOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      if (e > -1) cerr << "[" << e << "]";
      cerr << " does not contain an integer.\n";
      exit(Errno);
   }
   else
   {
      int tmp = SvIV(ParamVal);
      if (tmp < 0)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         if (d > -1) cerr << "[" << d << "]";
         if (e > -1) cerr << "[" << e << "]";
         cerr << " is negative!\n";
         exit(Errno);
      }
      
      return (unsigned) tmp;
   }
}

const char *PerlInput::getString(HashStruct Hash,char *ParamName,int a,int b,int c,int d,int e)
{
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d,e);
   int Errno=-5;
   if (!SvPOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      if (e > -1) cerr << "[" << e << "]";
      cerr << " does not contain a string.\n";
      exit(Errno);
   }
   else
   {
      return SvPV_nolen(ParamVal);
   }
}

void PerlInput::getVector(Vector &Vctr,HashStruct Hash,char *ParamName,
                          int a,int b,int c,int d)
{
   int Errno = -5;
   int len = Vctr.Dim();
   
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((av_len(ArrayPtr)+1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << len << ".\n";
      exit(Errno);
   }
   
   SV **ParamValPtr;
   for (int i=0;i<len;++i)
   {
      ParamValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if (SvTYPE(ParamVal) != SVt_NV)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         if (d > -1) cerr << "[" << d << "]";
         cerr << "[" << i << "] is not of type double\n";
         exit(Errno);
      }
      else
      {
         Vctr[i] = SvNV(ParamVal);
      }
   }
}

void PerlInput::getMatrix(Matrix &Mtrx,HashStruct Hash,char *ParamName,int a,int b,int c)
{
   int Errno = -5;
   int rows=Mtrx.Rows(),cols=Mtrx.Cols();
   
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((av_len(ArrayPtr)+1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << rows << ".\n";
      exit(Errno);
   }
   
   SV **RowValPtr, **ParamValPtr;
   SV *RowVal;
   AV *Row;
   for (int i=0;i<rows;++i)
   {
      RowValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }
      
      Row = (AV*) SvRV(RowVal);
      if ((av_len(Row)+1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row)+1 << " not " << cols << ".\n";
         exit(Errno);
      }
      
      for (int j=0;j<cols;++j)
      {
         ParamValPtr = av_fetch(Row,j,FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;
         
         if (SvTYPE(ParamVal) != SVt_NV)
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1) cerr << "[" << a << "]";
            if (b > -1) cerr << "[" << b << "]";
            if (c > -1) cerr << "[" << c << "]";
            cerr << "[" << i << "][" << j << "] is not of type double\n";
            exit(Errno);
         }
         else
         {
            Mtrx[i][j] = SvNV(ParamVal);
         }
      }
   }
}

void PerlInput::getIntVector(int *IntArry,unsigned len,HashStruct Hash,char *ParamName,
                                  int a,int b,int c,int d)
{
   int Errno = -5;
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((unsigned) (av_len(ArrayPtr)+1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << len << ".\n";
      exit(Errno);
   }
   
   SV **ParamValPtr;
   for (unsigned i=0;i<len;++i)
   {
      ParamValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if (SvTYPE(ParamVal) != SVt_IV)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         if (d > -1) cerr << "[" << d << "]";
         cerr << "[" << i << "] is not of type int\n";
         exit(Errno);
      }
      else
      {
         IntArry[i] = SvIV(ParamVal);
      }
   }
}

void PerlInput::getUnsignedVector(unsigned *UnsignedArry,unsigned len,HashStruct Hash,
                                  char *ParamName,int a,int b,int c,int d)
{
   int Errno = -5;
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c,d);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((unsigned) (av_len(ArrayPtr)+1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      if (d > -1) cerr << "[" << d << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << len << ".\n";
      exit(Errno);
   }
   
   SV **ParamValPtr;
   for (unsigned i=0;i<len;++i)
   {
      ParamValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if (SvTYPE(ParamVal) != SVt_IV)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         if (d > -1) cerr << "[" << d << "]";
         cerr << "[" << i << "] is not of type int\n";
         exit(Errno);
      }
      else
      {
         int tmp = SvIV(ParamVal);
         if (tmp < 0)
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1) cerr << "[" << a << "]";
            if (b > -1) cerr << "[" << b << "]";
            if (c > -1) cerr << "[" << c << "]";
            if (d > -1) cerr << "[" << d << "]";
            cerr << " is negative!\n";
            exit(Errno);
         }
      
         UnsignedArry[i] = (unsigned) tmp;
      }
   }
}

void PerlInput::getIntMatrix(int *IntMtrx,unsigned rows,unsigned cols,HashStruct Hash,
                                  char *ParamName,int a,int b,int c)
{
   int Errno = -5;
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((unsigned) (av_len(ArrayPtr)+1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << rows << ".\n";
      exit(Errno);
   }
   
   SV **RowValPtr, **ParamValPtr;
   SV *RowVal;
   AV *Row;
   for (unsigned i=0;i<rows;++i)
   {
      RowValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }
      
      Row = (AV*) SvRV(RowVal);
      if ((unsigned) (av_len(Row)+1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row)+1 << " not " << cols << ".\n";
         exit(Errno);
      }
      
      for (unsigned j=0;j<cols;++j)
      {
         ParamValPtr = av_fetch(Row,j,FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;
         
         if (SvTYPE(ParamVal) != SVt_IV)
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1) cerr << "[" << a << "]";
            if (b > -1) cerr << "[" << b << "]";
            if (c > -1) cerr << "[" << c << "]";
            cerr << "[" << i << "][" << j << "] is not of type double\n";
            exit(Errno);
         }
         else
         {
            IntMtrx[i*cols + j] = SvIV(ParamVal);
         }
      }
   }
}

void PerlInput::getUnsignedMatrix(unsigned *UnsignedMtrx,unsigned rows,unsigned cols,
                                  HashStruct Hash,char *ParamName,int a,int b,int c)
{
   int Errno = -5;
   SV *ParamVal = getScalar(Hash,ParamName,a,b,c);
   
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }
   
   
   AV *ArrayPtr = (AV*) SvRV(ParamVal);
   
   if ((unsigned) (av_len(ArrayPtr)+1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1) cerr << "[" << a << "]";
      if (b > -1) cerr << "[" << b << "]";
      if (c > -1) cerr << "[" << c << "]";
      cerr << " is an array of length " << av_len(ArrayPtr)+1 << " not " << rows << ".\n";
      exit(Errno);
   }
   
   SV **RowValPtr, **ParamValPtr;
   SV *RowVal;
   AV *Row;
   for (unsigned i=0;i<rows;++i)
   {
      RowValPtr = av_fetch(ArrayPtr,i,FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }
      
      Row = (AV*) SvRV(RowVal);
      if ((unsigned) (av_len(Row)+1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1) cerr << "[" << a << "]";
         if (b > -1) cerr << "[" << b << "]";
         if (c > -1) cerr << "[" << c << "]";
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row)+1 << " not " << cols << ".\n";
         exit(Errno);
      }
      
      for (unsigned j=0;j<cols;++j)
      {
         ParamValPtr = av_fetch(Row,j,FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;
         
         if (SvTYPE(ParamVal) != SVt_IV)
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1) cerr << "[" << a << "]";
            if (b > -1) cerr << "[" << b << "]";
            if (c > -1) cerr << "[" << c << "]";
            cerr << "[" << i << "][" << j << "] is not of type int.\n";
            exit(Errno);
         }
         else
         {
            int tmp = SvIV(ParamVal);
            if (tmp < 0)
            {
               cerr << "Error: Perl hash variable: " << Hash.Name
                    << "{" << ParamName << "}";
               if (a > -1) cerr << "[" << a << "]";
               if (b > -1) cerr << "[" << b << "]";
               if (c > -1) cerr << "[" << c << "]";
               cerr << " is negative!\n";
               exit(Errno);
            }
            
            UnsignedMtrx[i*cols + j] = (unsigned) tmp;
         }
      }
   }
}

