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
   PerlInput(char *datafile);
   virtual ~PerlInput();

   void Readfile(const char *datafile);
   void Readfile(const char *datafile,const char *prefix);

   int HashOK(char *HashName);
   int ParameterOK(HashStruct Hash,char *ParamName);
   int ParameterOK(char *HashName,char *ParamName)
   {return ParameterOK(getHash(HashName),ParamName);}

   HashStruct getHash(char *HashName);
   HashStruct getHash(HashStruct Hash,char *ParamName,
               int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);
   HashStruct getHash(char *HashName,char *ParamName,
               int a=-1,int b=-1,int c=-1,int d=-1,int e=-1)
   {return getHash(getHash(HashName),ParamName,a,b,c,d,e);}
   
   unsigned getArrayLength(HashStruct Hash,char *ParamName,
                           int a=-1,int b=-1,int c=-1,int d=-1);
   unsigned getArrayLength(char *HashName,char *ParamName,
                           int a=-1,int b=-1,int c=-1,int d=-1)
   {return getArrayLength(getHash(HashName),ParamName,a,b,c,d);}

   double getDouble(HashStruct Hash,char *ParamName,
                    int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);
   double getDouble(char *HashName,char *ParamName,
                    int a=-1,int b=-1,int c=-1,int d=-1,int e=-1)
   {return getDouble(getHash(HashName),ParamName,a,b,c,d,e);}
   
   int getInt(HashStruct Hash,char *ParamName,
              int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);
   int getInt(char *HashName,char *ParamName,
              int a=-1,int b=-1,int c=-1,int d=-1,int e=-1)
   {return getInt(getHash(HashName),ParamName,a,b,c,d,e);}

   int getUnsigned(HashStruct Hash,char *ParamName,
                   int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);
   int getUnsigned(char *HashName,char *ParamName,
                   int a=-1,int b=-1,int c=-1,int d=-1,int e=-1)
   {return getUnsigned(getHash(HashName),ParamName,a,b,c,d,e);}

   const char *getString(HashStruct Hash,char *ParamName,
                         int a=-1,int b=-1,int c=-1,int d=-1,int e=-1);
   const char *getString(char *HashName,char *ParamName,
                         int a=-1,int b=-1,int c=-1,int d=-1,int e=-1)
   {return getString(getHash(HashName),ParamName,a,b,c,d,e);}
   
   void getVector(Vector &Vect,HashStruct Hash,char *ParamName,
                  int a=-1,int b=-1,int c=-1,int d=-1);
   void getVector(Vector &Vect,char *HashName,char *ParamName,
                  int a=-1,int b=-1,int c=-1,int d=-1)
   {getVector(Vect,getHash(HashName),ParamName,a,b,c,d);}
  
   void getMatrix(Matrix &Mtrx,HashStruct Hash,char *ParamName,
                  int a=-1,int b=-1,int c=-1);
   void getMatrix(Matrix &Mtrx,char *HashName,char *ParamName,
                  int a=-1,int b=-1,int c=-1)
   {getMatrix(Mtrx,getHash(HashName),ParamName,a,b,c);}

   
   void getIntVector(int *IntArry,unsigned len,HashStruct Hash,char *ParamName,
                     int a=-1,int b=-1,int c=-1,int d=-1);
   void getIntVector(int *IntArry,unsigned len,char *HashName,char *ParamName,
                     int a=-1,int b=-1,int c=-1,int d=-1)
   {getIntVector(IntArry,len,getHash(HashName),ParamName,a,b,c,d);}
  
   void getUnsignedVector(unsigned *UnsignedArry,unsigned len,HashStruct Hash,char *ParamName,
                          int a=-1,int b=-1,int c=-1,int d=-1);
   void getUnsignedVector(unsigned *UnsignedArry,unsigned len,char *HashName,char *ParamName,
                          int a=-1,int b=-1,int c=-1,int d=-1)
   {getUnsignedVector(UnsignedArry,len,getHash(HashName),ParamName,a,b,c,d);}
   
   void getIntMatrix(int *IntMtrx,unsigned rows,unsigned cols,HashStruct Hash,char *ParamName,
                     int a=-1,int b=-1,int c=-1);
   void getIntMatrix(int *IntMtrx,unsigned rows,unsigned cols,char *HashName,char *ParamName,
                     int a=-1,int b=-1,int c=-1)
   {getIntMatrix(IntMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}

   void getUnsignedMatrix(unsigned *UnsignedMtrx,unsigned rows,unsigned cols,
                          HashStruct Hash,char *ParamName,int a=-1,int b=-1,int c=-1);   
   void getUnsignedMatrix(unsigned *UnsignedMtrx,unsigned rows,unsigned cols,
                          char *HashName,char *ParamName,int a=-1,int b=-1,int c=-1)
   {getUnsignedMatrix(UnsignedMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}
};

#endif
