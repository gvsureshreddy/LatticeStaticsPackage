#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "PerlInput.h"
#include <perl.h>

static PerlInterpreter* my_perl = 0;
static bool PerlInitialized = false;

void PerlInput::Initialize()
{
   ReconstructedInput_ << scientific << setprecision(14);
   const char* args[] = {"perl", "-W", "-e", "0"};
   int four = 4;
   char*** nu = 0;
   char** pargs = const_cast<char**>(args);
   if (my_perl != 0)
   {
      cerr << "PerlInput Error: Can only have one instance of PerlInput!\n";
      exit(-3);
   }

   if (!PerlInitialized)
   {
     PERL_SYS_INIT3(&four, &pargs, nu);
     atexit(&PerlInput::AtExit);
     PerlInitialized = true;
   }

   if (my_perl)
   {
     PL_perl_destruct_level = 1;
     perl_destruct(my_perl);
     perl_free(my_perl);
   }

   my_perl = perl_alloc();
   perl_construct(my_perl);
   PL_exit_flags |= PERL_EXIT_DESTRUCT_END;
   perl_parse(my_perl, 0, 4, pargs, (char**) 0);
   perl_run(my_perl);
}

void PerlInput::ClearHash(char const* const hashname)
{
   char tmp[256];
   SV* retval;

   sprintf(tmp, "undef %%%s;", hashname);
   retval = eval_pv(tmp, TRUE);
   // retval is always "undefined"
}

void PerlInput::EvaluateString(char const* const expression)
{
   SV* retval;
   retval = eval_pv(expression, TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in EvaluateString(), exiting...\n";
      exit(-5);
   }
}

void PerlInput::Readfile(char const* const datafile)
{
   char tmp[256];
   SV* retval;

   LastInputFileName_ = datafile;
   sprintf(tmp, "do '%s';", datafile);
   retval = eval_pv(tmp, TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }
}

void PerlInput::Readfile(char const* const datafile, char const* const prefix)
{
   char templt[]
      = "$___dmy = \"\"; open(___INPUT,\"%s\") or die \"can't open file for input from perl.\";"
        "while(<___INPUT>){ if ( /^%s/ ) {$_ =~ s/^%s//;"
        "$___dmy = $___dmy . $_;}}; close(___INPUT);";
   char temp[512];
   sprintf(temp, templt, datafile, prefix, prefix);

   LastInputFileName_ = datafile;

   SV* retval;
   retval = eval_pv(temp, TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }

   retval = eval_pv("eval $___dmy;", TRUE);
   if (!SvOK(retval))
   {
      cerr << "Perl Error in file: " << datafile << ", exiting...\n";
      exit(-2);
   }
}

PerlInput::PerlInput(char const* const datafile)
{
   Initialize();
   Readfile(datafile);
}

PerlInput::~PerlInput()
{
   PL_perl_destruct_level = 1;
   perl_destruct(my_perl);
   perl_free(my_perl);
   my_perl = 0;
}

void PerlInput::AtExit()
{
  PERL_SYS_TERM();
}

PerlInput::HashStruct PerlInput::getHash(char const* const HashName) const
{
   int Errno = -7;

   HashStruct Hash;
   Hash.Ptr = get_hv(HashName, FALSE);
   if (Hash.Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << HashName << " does not exist.\n";
      exit(Errno);
   }

   Hash.Name = HashName;
   return Hash;
}

PerlInput::HashStruct PerlInput::useHash(char const* const HashName) const
{
   HashStruct Hash;
   Hash.Ptr = 0;
   Hash.Name = HashName;

   return Hash;
}

PerlInput::HashStruct PerlInput::useHash(HashStruct const& Hash, char const* const ParamName)
const
{
   HashStruct NewHash;

   NewHash.Ptr = 0;
   stringstream tmp;
   tmp.str("");
   tmp << Hash.Name << "{" << ParamName << "}";
   NewHash.Name = tmp.str();

   return NewHash;
}

SV* const getScalar(PerlInput::HashStruct const& Hash, char const* const ParamName,
                    int const& a = -1, int const& b = -1, int const& c = -1, int const& d = -1,
                    int const& e = -1);

SV* const getScalar(PerlInput::HashStruct const& Hash, char const* const ParamName,
                    int const& a, int const& b, int const& c, int const& d, int const& e)
{
   int Errno = -5;

   SV** ParamValPtr = hv_fetch(Hash.Ptr, ParamName, strlen(ParamName), FALSE);
   if (ParamValPtr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "} does not exist.\n";
      exit(Errno);
   }

   SV* ParamVal = *ParamValPtr;
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
   AV* ArrayLevel0 = (AV*) SvRV(ParamVal);
   SV** ScalarLevel1Ptr = av_fetch(ArrayLevel0, a, FALSE);
   if (ScalarLevel1Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "] does not exist.\n";
      exit(Errno);
   }

   SV* ScalarLevel1 = *ScalarLevel1Ptr;
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

   AV* ArrayLevel1 = (AV*) SvRV(ScalarLevel1);
   SV** ScalarLevel2Ptr = av_fetch(ArrayLevel1, b, FALSE);
   if (ScalarLevel2Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "] does not exist.\n";
      exit(Errno);
   }

   SV* ScalarLevel2 = *ScalarLevel2Ptr;
   if (c <= -1)
   {
      return ScalarLevel2;
   }

   // move into array level 2
   if ((SvTYPE(ScalarLevel2) != SVt_RV) || (SvTYPE(SvRV(ScalarLevel2)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a
           << "][" << b << "]["
           << c << "] does not contain a perl referece (to an array).\n";
      exit(Errno);
   }
   AV* ArrayLevel2 = (AV*) SvRV(ScalarLevel2);
   SV** ScalarLevel3Ptr = av_fetch(ArrayLevel2, c, FALSE);
   if (ScalarLevel3Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "] does not exist.\n";
      exit(Errno);
   }

   SV* ScalarLevel3 = *ScalarLevel3Ptr;
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

   AV* ArrayLevel3 = (AV*) SvRV(ScalarLevel3);
   SV** ScalarLevel4Ptr = av_fetch(ArrayLevel3, d, FALSE);
   if (ScalarLevel4Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "]["
           << d << "] does not exist.\n";
      exit(Errno);
   }

   SV* ScalarLevel4 = *ScalarLevel4Ptr;
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

   AV* ArrayLevel4 = (AV*) SvRV(ScalarLevel4);
   SV** ScalarLevel5Ptr = av_fetch(ArrayLevel4, e, FALSE);
   if (ScalarLevel5Ptr == 0)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}[" << a << "]["
           << b << "][" << c << "][" << d << "]["
           << e << "] does not exist.\n";
      exit(Errno);
   }

   SV* ScalarLevel5 = *ScalarLevel5Ptr;
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
      return ScalarLevel5;
   }
}

PerlInput::HashStruct PerlInput::getHash(HashStruct const& Hash, char const* const ParamName,
                                         int const& a, int const& b, int const& c, int const& d,
                                         int const& e) const
{
   int Errno = -6;
   HashStruct NewHash;

   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d, e);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVHV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " does not contain a hash.\n";
      exit(Errno);
   }

   NewHash.Ptr = (HV*) SvRV(ParamVal);
   stringstream tmp;
   if (e > -1)
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}["
          << a << "][" << b << "][" << c << "][" << d << "][" << e << "]";
      NewHash.Name = tmp.str();
   }
   else if (d > -1)
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}["
          << a << "][" << b << "][" << c << "][" << d << "]";
      NewHash.Name = tmp.str();
   }
   else if (c > -1)
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}["
          << a << "][" << b << "][" << c << "]";
      NewHash.Name = tmp.str();
   }
   else if (b > -1)
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}["
          << a << "][" << b << "]";
      NewHash.Name = tmp.str();
   }
   else if (a > -1)
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}["
          << a << "]";
      NewHash.Name = tmp.str();
   }
   else
   {
      tmp.str("");
      tmp << Hash.Name << "{" << ParamName << "}";
      NewHash.Name = tmp.str();
   }
   return NewHash;
}

int PerlInput::HashOK(char const* const HashName) const
{
   int exists = 1;

   HV* HashPtr = get_hv(HashName, FALSE);
   if (HashPtr == 0)
   {
      exists = 0;
   }

   return exists;
}

int PerlInput::ParameterOK(HashStruct const& Hash, char const* const ParamName) const
{
   int exists = 1;

   SV** ParamValPtr = hv_fetch(Hash.Ptr, ParamName, strlen(ParamName), FALSE);
   if (ParamValPtr == 0)
   {
      exists = 0;
   }

   return exists;
}

int PerlInput::getArrayLength(HashStruct const& Hash, char const* const ParamName,
                              int const& a, int const& b, int const& c, int const& d) const
{
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d);
   int Errno = -4;
   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " does not contain an array.\n";
      exit(Errno);
   }

   return av_len((AV*) SvRV(ParamVal)) + 1;
}


double PerlInput::getDouble(HashStruct const& Hash, char const* const ParamName,
                            int const& a, int const& b, int const& c, int const& d, int const& e)
const
{
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d, e);
   int Errno = -5;
   if (!SvNOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " does not contain a double.\n";
      exit(Errno);
   }
   else
   {
      ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
      if (a > -1)
      {
         ReconstructedInput_ << "[" << a << "]";
      }
      if (b > -1)
      {
         ReconstructedInput_ << "[" << b << "]";
      }
      if (c > -1)
      {
         ReconstructedInput_ << "[" << c << "]";
      }
      if (d > -1)
      {
         ReconstructedInput_ << "[" << d << "]";
      }
      if (e > -1)
      {
         ReconstructedInput_ << "[" << e << "]";
      }
      ReconstructedInput_ << " = " << SvNV(ParamVal) << ";\n";
      return SvNV(ParamVal);
   }
}

double PerlInput::useDouble(double const& DefaultValue,
                            HashStruct const& Hash, char const* const ParamName, int const& a,
                            int const& b, int const& c, int const& d, int const& e) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   if (e > -1)
   {
      ReconstructedInput_ << "[" << e << "]";
   }
   ReconstructedInput_ << " = " << DefaultValue << ";                # Default Value\n";

   return DefaultValue;
}

void PerlInput::writeDouble(ostream& out, double const& Value,
                            HashStruct const& Hash, char const* const ParamName, int const& a,
                            int const& b, int const& c, int const& d, int const& e) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   if (e > -1)
   {
      out << "[" << e << "]";
   }
   out << " = " << Value << ";\n";
}

int PerlInput::getInt(HashStruct const& Hash, char const* const ParamName,
                      int const& a, int const& b, int const& c, int const& d, int const& e) const
{
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d, e);
   int Errno = -5;
   if (!SvIOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " does not contain an int.\n";
      exit(Errno);
   }
   else
   {
      ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
      if (a > -1)
      {
         ReconstructedInput_ << "[" << a << "]";
      }
      if (b > -1)
      {
         ReconstructedInput_ << "[" << b << "]";
      }
      if (c > -1)
      {
         ReconstructedInput_ << "[" << c << "]";
      }
      if (d > -1)
      {
         ReconstructedInput_ << "[" << d << "]";
      }
      if (e > -1)
      {
         ReconstructedInput_ << "[" << e << "]";
      }
      ReconstructedInput_ << " = " << SvIV(ParamVal) << ";\n";
      return SvIV(ParamVal);
   }
}

int PerlInput::useInt(int const& DefaultValue,
                      HashStruct const& Hash, char const* const ParamName, int const& a,
                      int const& b, int const& c, int const& d, int const& e) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   if (e > -1)
   {
      ReconstructedInput_ << "[" << e << "]";
   }
   ReconstructedInput_ << " = " << DefaultValue << ";                # Default Value\n";

   return DefaultValue;
}

void PerlInput::writeInt(ostream& out, int const& Value,
                         HashStruct const& Hash, char const* const ParamName, int const& a,
                         int const& b, int const& c, int const& d, int const& e) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   if (e > -1)
   {
      out << "[" << e << "]";
   }
   out << " = " << Value << ";\n";
}

int PerlInput::getPosInt(HashStruct const& Hash, char const* const ParamName,
                         int const& a, int const& b, int const& c, int const& d, int const& e)
const
{
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d, e);
   int Errno = -5;
   if (!SvIOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
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
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         if (e > -1)
         {
            cerr << "[" << e << "]";
         }
         cerr << " is negative!\n";
         exit(Errno);
      }

      ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
      if (a > -1)
      {
         ReconstructedInput_ << "[" << a << "]";
      }
      if (b > -1)
      {
         ReconstructedInput_ << "[" << b << "]";
      }
      if (c > -1)
      {
         ReconstructedInput_ << "[" << c << "]";
      }
      if (d > -1)
      {
         ReconstructedInput_ << "[" << d << "]";
      }
      if (e > -1)
      {
         ReconstructedInput_ << "[" << e << "]";
      }
      ReconstructedInput_ << " = " << tmp << ";\n";
      return tmp;
   }
}

int PerlInput::usePosInt(int const& DefaultValue, HashStruct const& Hash,
                         char const* const ParamName, int const& a, int const& b, int const& c,
                         int const& d, int const& e) const
{
   if (DefaultValue < 0)
   {
      cerr << "Error: Perl hash DEFAULT variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " is negative!\n";
      exit(34);
   }

   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   if (e > -1)
   {
      ReconstructedInput_ << "[" << e << "]";
   }
   ReconstructedInput_ << " = " << DefaultValue << ";                # Default Value\n";

   return DefaultValue;
}

void PerlInput::writePosInt(ostream& out, int const& Value, HashStruct const& Hash,
                            char const* const ParamName, int const& a, int const& b, int const& c,
                            int const& d, int const& e) const
{
   if (Value < 0)
   {
      cerr << "Error: writePosInt variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " is negative!\n";
      exit(34);
   }

   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   if (e > -1)
   {
      out << "[" << e << "]";
   }
   out << " = " << Value << ";\n";
}

char const* const PerlInput::getString(HashStruct const& Hash, char const* const ParamName,
                                       int const& a, int const& b, int const& c, int const& d,
                                       int const& e) const
{
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d, e);
   int Errno = -5;
   if (!SvPOK(ParamVal))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      if (e > -1)
      {
         cerr << "[" << e << "]";
      }
      cerr << " does not contain a string.\n";
      exit(Errno);
   }
   else
   {
      ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
      if (a > -1)
      {
         ReconstructedInput_ << "[" << a << "]";
      }
      if (b > -1)
      {
         ReconstructedInput_ << "[" << b << "]";
      }
      if (c > -1)
      {
         ReconstructedInput_ << "[" << c << "]";
      }
      if (d > -1)
      {
         ReconstructedInput_ << "[" << d << "]";
      }
      if (e > -1)
      {
         ReconstructedInput_ << "[" << e << "]";
      }
      ReconstructedInput_ << " = \"" << SvPV_nolen(ParamVal) << "\";\n";
      return SvPV_nolen(ParamVal);
   }
}

char const* const PerlInput::useString(char const* const DefaultValue, HashStruct const& Hash,
                                       char const* const ParamName, int const& a, int const& b,
                                       int const& c, int const& d, int const& e) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   if (e > -1)
   {
      ReconstructedInput_ << "[" << e << "]";
   }
   ReconstructedInput_ << " = \"" << DefaultValue << "\";                # Default Value\n";

   return DefaultValue;
}

void PerlInput::writeString(ostream& out, char const* const Value, HashStruct const& Hash,
                            char const* const ParamName, int const& a, int const& b,
                            int const& c, int const& d, int const& e) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   if (e > -1)
   {
      out << "[" << e << "]";
   }
   out << " = \"" << Value << "\";\n";
}

void PerlInput::getVector(Vector& Vctr, HashStruct const& Hash, char const* const ParamName,
                          int const& a, int const& b, int const& c, int const& d) const
{
   int Errno = -5;
   int len = Vctr.Dim();

   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }

   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << len << ".\n";
      exit(Errno);
   }

   SV** ParamValPtr;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < len; ++i)
   {
      ParamValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if ((SvTYPE(ParamVal) == SVt_NV) || (SvTYPE(ParamVal) == SVt_PVNV))
      {
         ReconstructedInput_ << SvNV(ParamVal) << ((i < len - 1) ? ", " : "];\n");
         Vctr[i] = SvNV(ParamVal);
      }
      else
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         cerr << "[" << i << "] is not of type double\n";
         exit(Errno);
      }
   }
}

void PerlInput::useVector(Vector const& DefaultVect, HashStruct const& Hash,
                          char const* const ParamName, int const& a, int const& b, int const& c,
                          int const& d) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < DefaultVect.Dim() - 1; ++i)
   {
      ReconstructedInput_ << DefaultVect[i] << ", ";
   }
   ReconstructedInput_ << DefaultVect[DefaultVect.Dim() - 1]
                       << "];                # Default Value\n";
}

void PerlInput::writeVector(ostream& out, Vector const& Vect, HashStruct const& Hash,
                            char const* const ParamName, int const& a, int const& b, int const& c,
                            int const& d) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   out << " = [";
   for (int i = 0; i < Vect.Dim() - 1; ++i)
   {
      out << Vect[i] << ", ";
   }
   out << Vect[Vect.Dim() - 1] << "];\n";
}

void PerlInput::getMatrix(Matrix& Mtrx, HashStruct const& Hash, char const* const ParamName,
                          int const& a, int const& b, int const& c) const
{
   int Errno = -5;
   int rows = Mtrx.Rows(), cols = Mtrx.Cols();
   int spaces = 0;

   SV* ParamVal = getScalar(Hash, ParamName, a, b, c);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }


   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << rows << ".\n";
      exit(Errno);
   }

   SV** RowValPtr, ** ParamValPtr;
   SV* RowVal;
   AV* Row;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [["; spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      RowValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }

      Row = (AV*) SvRV(RowVal);
      if ((av_len(Row) + 1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row) + 1 << " not " << cols << ".\n";
         exit(Errno);
      }

      for (int j = 0; j < cols; ++j)
      {
         ParamValPtr = av_fetch(Row, j, FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;

         if ((SvTYPE(ParamVal) == SVt_NV) || (SvTYPE(ParamVal) == SVt_PVNV))
         {
            ReconstructedInput_ << SvNV(ParamVal) << ((j < cols - 1) ? ", " : "]");
            Mtrx[i][j] = SvNV(ParamVal);
         }
         else
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            cerr << "[" << i << "][" << j << "] is not of type double\n";
            exit(Errno);
         }
      }
      if (i < rows - 1)
      {
         ReconstructedInput_ << ",\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "];\n";
      }
   }
}

void PerlInput::useMatrix(Matrix const& DefaultMtrx, HashStruct const& Hash,
                          char const* const ParamName, int const& a, int const& b, int const& c)
const
{
   int spaces = 0;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [["; spaces += 4;
   for (int i = 0; i < DefaultMtrx.Rows(); ++i)
   {
      for (int j = 0; j < DefaultMtrx.Cols() - 1; ++j)
      {
         ReconstructedInput_ << DefaultMtrx[i][j] << ", ";
      }
      ReconstructedInput_ << DefaultMtrx[i][DefaultMtrx.Cols() - 1] << "]";
      if (DefaultMtrx.Rows() - 1 == i)
      {
         ReconstructedInput_ << "];";
      }
      else
      {
         ReconstructedInput_ << ",";
      }
      if (0 == i)
      {
         ReconstructedInput_ << "                # Default Value";
      }
      if (DefaultMtrx.Rows() - 1 > i)
      {
         ReconstructedInput_ << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "\n";
      }
   }
}

void PerlInput::writeMatrix(ostream& out, Matrix const& Mtrx, HashStruct const& Hash,
                            char const* const ParamName, int const& a, int const& b, int const& c)
const
{
   int spaces = 0;
   out << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      out << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      out << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      out << "[" << c << "]"; spaces += 3;
   }
   out << " = [["; spaces += 4;
   for (int i = 0; i < Mtrx.Rows(); ++i)
   {
      for (int j = 0; j < Mtrx.Cols() - 1; ++j)
      {
         out << Mtrx[i][j] << ", ";
      }
      out << Mtrx[i][Mtrx.Cols() - 1] << "]";
      if (Mtrx.Rows() - 1 == i)
      {
         out << "];";
      }
      else
      {
         out << ",";
      }
      if (Mtrx.Rows() - 1 > i)
      {
         out << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            out << " ";
         }
         out << "[";
      }
      else
      {
         out << "\n";
      }
   }
}

void PerlInput::getIntVector(int* const IntArry, int const& len, HashStruct const& Hash,
                             char const* const ParamName, int const& a, int const& b,
                             int const& c, int const& d) const
{
   int Errno = -5;
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }

   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << len << ".\n";
      exit(Errno);
   }

   SV** ParamValPtr;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < len; ++i)
   {
      ParamValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if ((SvTYPE(ParamVal) == SVt_IV) || (SvTYPE(ParamVal) == SVt_PVIV))
      {
         ReconstructedInput_ << SvIV(ParamVal) << ((i < len - 1) ? ", " : "];\n");
         IntArry[i] = SvIV(ParamVal);
      }
      else
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         cerr << "[" << i << "] is not of type int\n";
         exit(Errno);
      }
   }
}

void PerlInput::useIntVector(int const* const DefaultIntArry, int const& len,
                             HashStruct const& Hash, char const* const ParamName, int const& a,
                             int const& b, int const& c, int const& d) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < len - 1; ++i)
   {
      ReconstructedInput_ << DefaultIntArry[i] << ", ";
   }
   ReconstructedInput_ << DefaultIntArry[len - 1]
                       << "];                # Default Value\n";
}

void PerlInput::writeIntVector(ostream& out, int const* const IntArry, int const& len,
                               HashStruct const& Hash, char const* const ParamName, int const& a,
                               int const& b, int const& c, int const& d) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   out << " = [";
   for (int i = 0; i < len - 1; ++i)
   {
      out << IntArry[i] << ", ";
   }
   out << IntArry[len - 1] << "];\n";
}

void PerlInput::getPosIntVector(int* const PosIntArry, int const& len, HashStruct const& Hash,
                                char const* const ParamName, int const& a, int const& b,
                                int const& c, int const& d) const
{
   int Errno = -5;
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c, d);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }

   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != len)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << len << ".\n";
      exit(Errno);
   }

   SV** ParamValPtr;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < len; ++i)
   {
      ParamValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      ParamVal = *ParamValPtr;
      if ((SvTYPE(ParamVal) == SVt_IV) || (SvTYPE(ParamVal) == SVt_PVIV))
      {
         int tmp = SvIV(ParamVal);
         if (tmp < 0)
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            if (d > -1)
            {
               cerr << "[" << d << "]";
            }
            cerr << " is negative!\n";
            exit(Errno);
         }

         ReconstructedInput_ << tmp << ((i < len - 1) ? ", " : "];\n");
         PosIntArry[i] = tmp;
      }
      else
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         cerr << "[" << i << "] is not of type int\n";
         exit(Errno);
      }
   }
}

void PerlInput::usePosIntVector(int const* const DefaultPosIntArry, int const& len,
                                HashStruct const& Hash, char const* const ParamName, int const& a,
                                int const& b, int const& c, int const& d) const
{
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]";
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]";
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]";
   }
   if (d > -1)
   {
      ReconstructedInput_ << "[" << d << "]";
   }
   ReconstructedInput_ << " = [";
   for (int i = 0; i < len - 1; ++i)
   {
      if (DefaultPosIntArry[i] < 0)
      {
         cerr << "Error: Perl hash DEFAULT variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         cerr << "[" << i << "]";
         cerr << " is negative!\n";
         exit(34);
      }
      ReconstructedInput_ << DefaultPosIntArry[i] << ", ";
   }
   if (DefaultPosIntArry[len - 1] < 0)
   {
      cerr << "Error: Perl hash DEFAULT variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << "[" << len - 1 << "]";
      cerr << " is negative!\n";
      exit(34);
   }
   ReconstructedInput_ << DefaultPosIntArry[len - 1]
                       << "];                # Default Value\n";
}

void PerlInput::writePosIntVector(ostream& out, int const* const PosIntArry,
                                  int const& len, HashStruct const& Hash,
                                  char const* const ParamName, int const& a, int const& b,
                                  int const& c, int const& d) const
{
   out << "$" << Hash.Name << "{" << ParamName << "}";
   if (a > -1)
   {
      out << "[" << a << "]";
   }
   if (b > -1)
   {
      out << "[" << b << "]";
   }
   if (c > -1)
   {
      out << "[" << c << "]";
   }
   if (d > -1)
   {
      out << "[" << d << "]";
   }
   out << " = [";
   for (int i = 0; i < len - 1; ++i)
   {
      if (PosIntArry[i] < 0)
      {
         cerr << "Error: writePosIntVector() variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         if (d > -1)
         {
            cerr << "[" << d << "]";
         }
         cerr << "[" << i << "]";
         cerr << " is negative!\n";
         exit(34);
      }
      out << PosIntArry[i] << ", ";
   }
   if (PosIntArry[len - 1] < 0)
   {
      cerr << "Error: writePosIntVector variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      if (d > -1)
      {
         cerr << "[" << d << "]";
      }
      cerr << "[" << len - 1 << "]";
      cerr << " is negative!\n";
      exit(34);
   }
   out << PosIntArry[len - 1] << "];\n";
}

void PerlInput::getIntMatrix(int* const IntMtrx, int const& rows, int const& cols,
                             HashStruct const& Hash, char const* const ParamName, int const& a,
                             int const& b, int const& c) const
{
   int Errno = -5;
   int spaces = 0;
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }


   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << rows << ".\n";
      exit(Errno);
   }

   SV** RowValPtr, ** ParamValPtr;
   SV* RowVal;
   AV* Row;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [[";
   spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      RowValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }

      Row = (AV*) SvRV(RowVal);
      if ((av_len(Row) + 1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row) + 1 << " not " << cols << ".\n";
         exit(Errno);
      }

      for (int j = 0; j < cols; ++j)
      {
         ParamValPtr = av_fetch(Row, j, FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;

         if ((SvTYPE(ParamVal) == SVt_IV) || (SvTYPE(ParamVal) == SVt_PVIV))
         {
            ReconstructedInput_ << SvIV(ParamVal) << ((j < cols - 1) ? ", " : "]");
            IntMtrx[i * cols + j] = SvIV(ParamVal);
         }
         else
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            cerr << "[" << i << "][" << j << "] is not of type double\n";
            exit(Errno);
         }
      }
      if (i < rows - 1)
      {
         ReconstructedInput_ << ",\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "];\n";
      }
   }
}

void PerlInput::useIntMatrix(int const* const DefaultIntMtrx, int const& rows, int const& cols,
                             HashStruct const& Hash, char const* const ParamName, int const& a,
                             int const& b, int const& c) const
{
   int spaces = 0;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [["; spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      for (int j = 0; j < cols - 1; ++j)
      {
         ReconstructedInput_ << DefaultIntMtrx[i * cols + j] << ", ";
      }
      ReconstructedInput_ << DefaultIntMtrx[i * cols + cols - 1] << "]";
      if (rows - 1 == i)
      {
         ReconstructedInput_ << "];";
      }
      else
      {
         ReconstructedInput_ << ",";
      }
      if (0 == i)
      {
         ReconstructedInput_ << "                # Default Value";
      }
      if (rows - 1 > i)
      {
         ReconstructedInput_ << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "\n";
      }
   }
}

void PerlInput::writeIntMatrix(ostream& out, int const* const IntMtrx, int const& rows,
                               int const& cols, HashStruct const& Hash,
                               char const* const ParamName, int const& a, int const& b,
                               int const& c) const
{
   int spaces = 0;
   out << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      out << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      out << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      out << "[" << c << "]"; spaces += 3;
   }
   out << " = [["; spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      for (int j = 0; j < cols - 1; ++j)
      {
         out << IntMtrx[i * cols + j] << ", ";
      }
      out << IntMtrx[i * cols + cols - 1] << "]";
      if (rows - 1 == i)
      {
         out << "];";
      }
      else
      {
         out << ",";
      }
      if (rows - 1 > i)
      {
         out << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            out << " ";
         }
         out << "[";
      }
      else
      {
         out << "\n";
      }
   }
}

void PerlInput::getPosIntMatrix(int* const PosIntMtrx, int const& rows, int const& cols,
                                HashStruct const& Hash, char const* const ParamName,
                                int const& a, int const& b, int const& c) const
{
   int Errno = -5;
   int spaces = 0;
   SV* ParamVal = getScalar(Hash, ParamName, a, b, c);

   if ((SvTYPE(ParamVal) != SVt_RV) || (SvTYPE(SvRV(ParamVal)) != SVt_PVAV))
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " does not contain a reference to an array.\n";
      exit(Errno);
   }


   AV* ArrayPtr = (AV*) SvRV(ParamVal);

   if ((av_len(ArrayPtr) + 1) != rows)
   {
      cerr << "Error: Perl hash variable: " << Hash.Name
           << "{" << ParamName << "}";
      if (a > -1)
      {
         cerr << "[" << a << "]";
      }
      if (b > -1)
      {
         cerr << "[" << b << "]";
      }
      if (c > -1)
      {
         cerr << "[" << c << "]";
      }
      cerr << " is an array of length " << av_len(ArrayPtr) + 1 << " not " << rows << ".\n";
      exit(Errno);
   }

   SV** RowValPtr, ** ParamValPtr;
   SV* RowVal;
   AV* Row;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [[";
   spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      RowValPtr = av_fetch(ArrayPtr, i, FALSE); // Should always be ok because of len check.
      RowVal = *RowValPtr;
      if ((SvTYPE(RowVal) != SVt_RV) || (SvTYPE(SvRV(RowVal)) != SVt_PVAV))
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "] is not an array.\n";
         exit(Errno);
      }

      Row = (AV*) SvRV(RowVal);
      if ((av_len(Row) + 1) != cols)
      {
         cerr << "Error: Perl hash variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "]";
         cerr << " is an array of length " << av_len(Row) + 1 << " not " << cols << ".\n";
         exit(Errno);
      }

      for (int j = 0; j < cols; ++j)
      {
         ParamValPtr = av_fetch(Row, j, FALSE); // Should be ok because of len check
         ParamVal = *ParamValPtr;

         if ((SvTYPE(ParamVal) == SVt_IV) || (SvTYPE(ParamVal) == SVt_PVIV))
         {
            int tmp = SvIV(ParamVal);
            if (tmp < 0)
            {
               cerr << "Error: Perl hash variable: " << Hash.Name
                    << "{" << ParamName << "}";
               if (a > -1)
               {
                  cerr << "[" << a << "]";
               }
               if (b > -1)
               {
                  cerr << "[" << b << "]";
               }
               if (c > -1)
               {
                  cerr << "[" << c << "]";
               }
               cerr << " is negative!\n";
               exit(Errno);
            }

            ReconstructedInput_ << tmp << ((j < cols - 1) ? ", " : "]");
            PosIntMtrx[i * cols + j] = tmp;
         }
         else
         {
            cerr << "Error: Perl hash variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            cerr << "[" << i << "][" << j << "] is not of type int.\n";
            exit(Errno);
         }
      }
      if (i < rows - 1)
      {
         ReconstructedInput_ << ",\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "];\n";
      }
   }
}

void PerlInput::usePosIntMatrix(int const* const DefaultPosIntMtrx, int const& rows,
                                int const& cols, HashStruct const& Hash,
                                char const* const ParamName, int const& a, int const& b,
                                int const& c) const
{
   int spaces = 0;
   ReconstructedInput_ << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      ReconstructedInput_ << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      ReconstructedInput_ << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      ReconstructedInput_ << "[" << c << "]"; spaces += 3;
   }
   ReconstructedInput_ << " = [["; spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      for (int j = 0; j < cols - 1; ++j)
      {
         if (DefaultPosIntMtrx[i * cols + j] < 0)
         {
            cerr << "Error: Perl hash DEFAULT variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            cerr << "[" << i << "][" << j << "]";
            cerr << " is negative!\n";
            exit(34);
         }
         ReconstructedInput_ << DefaultPosIntMtrx[i * cols + j] << ", ";
      }
      if (DefaultPosIntMtrx[i * cols + cols - 1] < 0)
      {
         cerr << "Error: Perl hash DEFAULT variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "][" << cols - 1 << "]";
         cerr << " is negative!\n";
         exit(34);
      }
      ReconstructedInput_ << DefaultPosIntMtrx[i * cols + cols - 1] << "]";
      if (rows - 1 == i)
      {
         ReconstructedInput_ << "];";
      }
      else
      {
         ReconstructedInput_ << ",";
      }
      if (0 == i)
      {
         ReconstructedInput_ << "                # Default Value";
      }
      if (rows - 1 > i)
      {
         ReconstructedInput_ << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            ReconstructedInput_ << " ";
         }
         ReconstructedInput_ << "[";
      }
      else
      {
         ReconstructedInput_ << "\n";
      }
   }
}

void PerlInput::writePosIntMatrix(ostream& out, int const* const PosIntMtrx, int const& rows,
                                  int const& cols, HashStruct const& Hash,
                                  char const* const ParamName, int const& a, int const& b,
                                  int const& c) const
{
   int spaces = 0;
   out << "$" << Hash.Name << "{" << ParamName << "}";
   spaces += 3 + Hash.Name.length() + strlen(ParamName);
   if (a > -1)
   {
      out << "[" << a << "]"; spaces += 3;
   }
   if (b > -1)
   {
      out << "[" << b << "]"; spaces += 3;
   }
   if (c > -1)
   {
      out << "[" << c << "]"; spaces += 3;
   }
   out << " = [["; spaces += 4;
   for (int i = 0; i < rows; ++i)
   {
      for (int j = 0; j < cols - 1; ++j)
      {
         if (PosIntMtrx[i * cols + j] < 0)
         {
            cerr << "Error: writePosIntMatrix() variable: " << Hash.Name
                 << "{" << ParamName << "}";
            if (a > -1)
            {
               cerr << "[" << a << "]";
            }
            if (b > -1)
            {
               cerr << "[" << b << "]";
            }
            if (c > -1)
            {
               cerr << "[" << c << "]";
            }
            cerr << "[" << i << "][" << j << "]";
            cerr << " is negative!\n";
            exit(34);
         }
         out << PosIntMtrx[i * cols + j] << ", ";
      }
      if (PosIntMtrx[i * cols + cols - 1] < 0)
      {
         cerr << "Error: writePosIntMatrix() variable: " << Hash.Name
              << "{" << ParamName << "}";
         if (a > -1)
         {
            cerr << "[" << a << "]";
         }
         if (b > -1)
         {
            cerr << "[" << b << "]";
         }
         if (c > -1)
         {
            cerr << "[" << c << "]";
         }
         cerr << "[" << i << "][" << cols - 1 << "]";
         cerr << " is negative!\n";
         exit(34);
      }
      out << PosIntMtrx[i * cols + cols - 1] << "]";
      if (rows - 1 == i)
      {
         out << "];";
      }
      else
      {
         out << ",";
      }
      if (rows - 1 > i)
      {
         out << "\n";
         for (int z = 0; z < spaces; ++z)
         {
            out << " ";
         }
         out << "[";
      }
      else
      {
         out << "\n";
      }
   }
}
