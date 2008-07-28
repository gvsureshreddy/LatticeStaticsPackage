#ifndef RSE__PerlInput
#define RSE__PerlInput

#include <string>
#include <sstream>

#include <Matrix.h>
#include <Vector.h>
#include <EXTERN.h>
#include <perl.h>

using namespace std;

class PerlInput
{
private:
   string LastInputFileName_;
   mutable ostringstream ReconstructedInput_;

   void Initialize();

public:
   
   struct HashStruct
   {
      HV* Ptr;
      string Name;
   };
   
   PerlInput() {Initialize();}
   PerlInput(char const* const datafile);
   virtual ~PerlInput();
   
   void Readfile(char const* const datafile);
   void Readfile(char const* const datafile,char const* const prefix);
   char const* const LastInputFileName() const {return LastInputFileName_.c_str();}
   void EvaluateString(char const* const expression);
   void ClearHash(char const* const hashname);

   string ReconstructedInput() const {return ReconstructedInput_.str();}
   void EndofInputSection() const {ReconstructedInput_ << "\n";}
   
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

   HashStruct useHash(char const* const HashName) const;
   HashStruct useHash(HashStruct const& Hash,char const* const ParamName) const;
   HashStruct useHash(char const* const HashName,char const* const ParamName) const
   {return useHash(useHash(HashName),ParamName);}
   
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

   double useDouble(double const& DefaultValue,
                    HashStruct const& Hash,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const;
   double useDouble(double const& DefaultValue,
                    char const* const HashName,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const
   {return useDouble(DefaultValue,useHash(HashName),ParamName,a,b,c,d,e);}

   void writeDouble(ostream& out,double const& Value,
                    HashStruct const& Hash,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const;
   void writeDouble(ostream& out,double const& Value,
                    char const* const HashName,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const
   {writeDouble(out,Value,useHash(HashName),ParamName,a,b,c,d,e);}

   int getInt(HashStruct const& Hash,char const* const ParamName,
              int const& a=-1,int const&  b=-1,int const& c=-1,int const& d=-1,
              int const& e=-1) const;
   int getInt(char const* const HashName,char const* const ParamName,
              int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,int const& e=-1)
      const
   {return getInt(getHash(HashName),ParamName,a,b,c,d,e);}

   int useInt(int const& DefaultValue,
              HashStruct const& Hash,char const* const ParamName,
              int const& a=-1,int const&  b=-1,int const& c=-1,int const& d=-1,
              int const& e=-1) const;
   int useInt(int const& DefaultValue,
              char const* const HashName,char const* const ParamName,
              int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,int const& e=-1)
      const
   {return useInt(DefaultValue,useHash(HashName),ParamName,a,b,c,d,e);}

   void writeInt(ostream& out,int const& Value,
                 HashStruct const& Hash,char const* const ParamName,
                 int const& a=-1,int const&  b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const;
   void writeInt(ostream& out,int const& Value,
                 char const* const HashName,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const
   {writeInt(out,Value,useHash(HashName),ParamName,a,b,c,d,e);}

   int getPosInt(HashStruct const& Hash,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const;
   int getPosInt(char const* const HashName,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const
   {return getPosInt(getHash(HashName),ParamName,a,b,c,d,e);}

   int usePosInt(int const& DefaultValue,
                 HashStruct const& Hash,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const;
   int usePosInt(int const& DefaultValue,
                 char const* const HashName,char const* const ParamName,
                 int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                 int const& e=-1) const
   {return usePosInt(DefaultValue,useHash(HashName),ParamName,a,b,c,d,e);}

   void writePosInt(ostream& out,int const& Value,
                    HashStruct const& Hash,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const;
   void writePosInt(ostream& out,int const& Value,
                    char const* const HashName,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const
   {writePosInt(out,Value,useHash(HashName),ParamName,a,b,c,d,e);}

   char const* const getString(HashStruct const& Hash,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const;
   char const* const getString(char const* const HashName,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const
   {return getString(getHash(HashName),ParamName,a,b,c,d,e);}

   char const* const useString(char const* const DefaultValue,
                               HashStruct const& Hash,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const;
   char const* const useString(char const* const DefaultValue,
                               char const* const HashName,char const* const ParamName,
                               int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                               int const& e=-1) const
   {return useString(DefaultValue,useHash(HashName),ParamName,a,b,c,d,e);}

   void writeString(ostream& out,char const* const Value,
                    HashStruct const& Hash,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const;
   void writeString(ostream& out,char const* const Value,
                    char const* const HashName,char const* const ParamName,
                    int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1,
                    int const& e=-1) const
   {writeString(out,Value,useHash(HashName),ParamName,a,b,c,d,e);}

   void getVector(Vector& Vect,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const;
   void getVector(Vector& Vect,char const* const HashName,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const
   {getVector(Vect,getHash(HashName),ParamName,a,b,c,d);}

   void useVector(Vector const& DefaultVect,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const;
   void useVector(Vector const& DefaultVect,char const* const HashName,
                  char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1,
                  int const& d=-1) const
   {useVector(DefaultVect,useHash(HashName),ParamName,a,b,c,d);}

   void writeVector(ostream& out,Vector const& Vect,HashStruct const& Hash,
                    char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1,
                    int const& d=-1) const;
   void writeVector(ostream& out,Vector const& Vect,char const* const HashName,
                    char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1,
                    int const& d=-1) const
   {writeVector(out,Vect,useHash(HashName),ParamName,a,b,c,d);}

   void getMatrix(Matrix& Mtrx,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1) const;
   void getMatrix(Matrix& Mtrx,char const* const HashName,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1) const
   {getMatrix(Mtrx,getHash(HashName),ParamName,a,b,c);}

   void useMatrix(Matrix const& DefaultMtrx,HashStruct const& Hash,char const* const ParamName,
                  int const& a=-1,int const& b=-1,int const& c=-1) const;
   void useMatrix(Matrix const& DefaultMtrx,char const* const HashName,
                  char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1)
      const
   {useMatrix(DefaultMtrx,useHash(HashName),ParamName,a,b,c);}

   void writeMatrix(ostream& out,Matrix const& Mtrx,HashStruct const& Hash,
                    char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1)
      const;
   void writeMatrix(ostream& out,Matrix const& Mtrx,char const* const HashName,
                    char const* const ParamName,int const& a=-1,int const& b=-1,int const& c=-1)
      const
   {writeMatrix(out,Mtrx,useHash(HashName),ParamName,a,b,c);}

   void getIntVector(int* const IntArry,int const& len,HashStruct const& Hash,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const;
   void getIntVector(int* const IntArry,int const& len,char const* const HashName,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const
   {getIntVector(IntArry,len,getHash(HashName),ParamName,a,b,c,d);}

   void useIntVector(int const* const DefaultIntArry,int const& len,HashStruct const& Hash,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const;
   void useIntVector(int const* const DefaultIntArry,int const& len,char const* const HashName,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1,int const& d=-1) const
   {useIntVector(DefaultIntArry,len,useHash(HashName),ParamName,a,b,c,d);}

   void writeIntVector(ostream& out,int const* const IntArry,int const& len,
                       HashStruct const& Hash,char const* const ParamName,int const& a=-1,
                       int const& b=-1,int const& c=-1,int const& d=-1) const;
   void writeIntVector(ostream& out,int const* const IntArry,int const& len,
                       char const* const HashName,char const* const ParamName,int const& a=-1,
                       int const& b=-1,int const& c=-1,int const& d=-1) const
   {writeIntVector(out,IntArry,len,useHash(HashName),ParamName,a,b,c,d);}

   void getPosIntVector(int* const PosIntArry,int const& len,HashStruct const& Hash,
                        char const* const ParamName,int const& a=-1,int const& b=-1,
                        int const& c=-1,int const& d=-1) const;
   void getPosIntVector(int* const PosIntArry,int const& len,char const* const HashName,
                        char const* const ParamName,int const& a=-1,int const& b=-1,
                        int const& c=-1,int const& d=-1) const
   {getPosIntVector(PosIntArry,len,getHash(HashName),ParamName,a,b,c,d);}

   void usePosIntVector(int const* const DefaultPosIntArry,int const& len,
                        HashStruct const& Hash,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const;
   void usePosIntVector(int const* const DefaultPosIntArry,int const& len,
                        char const* const HashName,char const* const ParamName,int const& a=-1,
                        int const& b=-1,int const& c=-1,int const& d=-1) const
   {usePosIntVector(DefaultPosIntArry,len,useHash(HashName),ParamName,a,b,c,d);}

   void writePosIntVector(ostream& out,int const* const PosIntArry,int const& len,
                          HashStruct const& Hash,char const* const ParamName,
                          int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1)
      const;
   void writePosIntVector(ostream& out,int const* const PosIntArry,int const& len,
                          char const* const HashName,char const* const ParamName,
                          int const& a=-1,int const& b=-1,int const& c=-1,int const& d=-1) const
   {writePosIntVector(out,PosIntArry,len,useHash(HashName),ParamName,a,b,c,d);}

   void getIntMatrix(int* const IntMtrx,int const& rows,int const& cols,HashStruct const& Hash,
                     char const* const ParamName,int const& a=-1,int const& b=-1,
                     int const& c=-1) const;
   void getIntMatrix(int* const IntMtrx,int const& rows,int const& cols,
                     char const* const HashName,char const* const ParamName,
                     int const& a=-1,int const& b=-1,int const& c=-1) const
   {getIntMatrix(IntMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}

   void useIntMatrix(int const* const DefaultIntMtrx,int const& rows,int const& cols,
                     HashStruct const& Hash,char const* const ParamName,int const& a=-1,
                     int const& b=-1,int const& c=-1) const;
   void useIntMatrix(int const* const DefaultIntMtrx,int const& rows,int const& cols,
                     char const* const HashName,char const* const ParamName,
                     int const& a=-1,int const& b=-1,int const& c=-1) const
   {useIntMatrix(DefaultIntMtrx,rows,cols,useHash(HashName),ParamName,a,b,c);}

   void writeIntMatrix(ostream& out,int const* const IntMtrx,int const& rows,int const& cols,
                       HashStruct const& Hash,char const* const ParamName,int const& a=-1,
                       int const& b=-1,int const& c=-1) const;
   void writeIntMatrix(ostream& out,int const* const IntMtrx,int const& rows,int const& cols,
                       char const* const HashName,char const* const ParamName,
                       int const& a=-1,int const& b=-1,int const& c=-1) const
   {writeIntMatrix(out,IntMtrx,rows,cols,useHash(HashName),ParamName,a,b,c);}

   void getPosIntMatrix(int* const PosIntMtrx,int const& rows,int const& cols,
                        HashStruct const& Hash,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const;
   void getPosIntMatrix(int* const PosIntMtrx,int const& rows,int const& cols,
                        char const* const HashName,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const
   {getPosIntMatrix(PosIntMtrx,rows,cols,getHash(HashName),ParamName,a,b,c);}

   void usePosIntMatrix(int const* const DefaultPosIntMtrx,int const& rows,int const& cols,
                        HashStruct const& Hash,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const;
   void usePosIntMatrix(int const* const DefaultPosIntMtrx,int const& rows,int const& cols,
                        char const* const HashName,char const* const ParamName,
                        int const& a=-1,int const& b=-1,int const& c=-1) const
   {usePosIntMatrix(DefaultPosIntMtrx,rows,cols,useHash(HashName),ParamName,a,b,c);}

   void writePosIntMatrix(ostream& out,int const* const PosIntMtrx,int const& rows,
                          int const& cols,HashStruct const& Hash,char const* const ParamName,
                          int const& a=-1,int const& b=-1,int const& c=-1) const;
   void writePosIntMatrix(ostream& out,int const* const PosIntMtrx,int const& rows,
                          int const& cols,char const* const HashName,
                          char const* const ParamName,int const& a=-1,int const& b=-1,
                          int const& c=-1) const
   {writePosIntMatrix(out,PosIntMtrx,rows,cols,useHash(HashName),ParamName,a,b,c);}
};

#endif
