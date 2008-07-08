#ifndef __PerlInput
#define __PerlInput

#include <Matrix.h>
#include <Vector.h>
#include <EXTERN.h>
#include <perl.h>

using namespace std;

#define HASH_NAME_LEN 256

class PerlInput
{
private:
   void Initialize();
   
public:
   
   struct HashStruct
   {
      HV *Ptr;
      char Name[HASH_NAME_LEN];
   };
   
   PerlInput() {Initialize();}
   PerlInput(char const* const datafile);
   virtual ~PerlInput();
   
   void Readfile(char const* const datafile);
   void Readfile(char const* const datafile,char const* const prefix);
   void ClearHash(char const* const hashname);
   
   int HashOK(char const* const HashName) const;
   int ParameterOK(HashStruct const& Hash,char const* const ParamName) const;
   int ParameterOK(char const* const HashName,char const* const ParamName) const
   {return ParameterOK(getHash(HashName),ParamName);}
   
   HashStruct getHash(char const* const HashName) const;
   HashStruct getHash(HashStruct const& Hash,char const* const ParamName,
                      int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                      int const& e=-1) const;
   HashStruct getHash(char const* const HashName,char const* const ParamName,
                      int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                      int const& e=-1) const
   {return getHash(getHash(HashName),ParamName,a,b,c,d,e);}
   
   int getArrayLength(HashStruct const& Hash,char const* const ParamName,
                      int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const;
   int getArrayLength(char const* const HashName,char const* const ParamName,
                      int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const
   {return getArrayLength(getHash(HashName),ParamName,a,b,c,d);}
   
   double getDouble(HashStruct const& Hash,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const;
   double getDouble(char const* const HashName,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const
   {return getDouble(getHash(HashName),ParamName,a,b,c,d,e);}
   
   int getInt(HashStruct const& Hash,char const* const ParamName,
              int const& a=-1,int const&  b=-1,int const& c=-1,int const& d=-1,
              int const& e=-1) const;
   int getInt(char const* const HashName,char const* const ParamName,
              int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,int const& e=-1)
      const
   {return getInt(getHash(HashName),ParamName,a,b,c,d,e);}
   
   int getPosInt(HashStruct const& Hash,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const;
   int getPosInt(char const* const HashName,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const
   {return getPosInt(getHash(HashName),ParamName,a,b,c,d,e);}
   
   char const* const getString(HashStruct const& Hash,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const;
   char const* const getString(char const* const HashName,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const
   {return getString(getHash(HashName),ParamName,a,b,c,d,e);}
   
   void getVector(Vector& Vect,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const;
   void getVector(Vector& Vect,char const* const HashName,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const
   {getVector(Vect,getHash(HashName),ParamName,a,b,c,d);}
   
   void getMatrix(Matrix& Mtrx,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1) const;
   void getMatrix(Matrix& Mtrx,char const* const HashName,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1) const
   {getMatrix(Mtrx,getHash(HashName),ParamName,a,b,c);}
   
   
   void getIntVector(int* const IntArry,int const& len,HashStruct const& Hash,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const;
   void getIntVector(int* const IntArry,int const& len,char const* const HashName,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const
   {getIntVector(IntArry,len,getHash(HashName),ParamName,a,b,c,d);}
   
   void getPosIntVector(int* const PosIntArry,int const& len,HashStruct const& Hash,
                        char const* const ParamName,int const& a=-1,int const& b=-1,
                        int const& c=-1,int const& d=-1) const;
   void getPosIntVector(int* const PosIntArry,int const& len,char const* const HashName,
                        char const* const ParamName,int const& a=-1,int const& b=-1,
                        int const& c=-1,int const& d=-1) const
   {getPosIntVector(PosIntArry,len,getHash(HashName),ParamName,a,b,c,d);}
   
   void getIntMatrix(int* const IntMtrx,int const& rows,int const& cols,HashStruct const& Hash,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1) const;
   void getIntMatrix(int* const IntMtrx,int const& rows,int const& cols,
                     char const* const HashName,char const* const ParamName,
                     int const& a=-1,int const& b=-1,int const& c=-1) const
   {getIntMatrix(IntMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}
   
   void getPosIntMatrix(int* const PosIntMtrx,int const& rows,int const& cols,
                        HashStruct const& Hash,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const;
   void getPosIntMatrix(int* const PosIntMtrx,int const& rows,int const& cols,
                        char const* const HashName,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const
   {getPosIntMatrix(PosIntMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}
};

#endif
