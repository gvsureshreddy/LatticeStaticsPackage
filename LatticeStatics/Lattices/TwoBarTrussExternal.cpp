#include "TwoBarTrussExternal.h"
#include <fstream>

using namespace std;

TwoBarTrussExternal::~TwoBarTrussExternal()
{
   cout << "TwoBarTrussExternal Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

TwoBarTrussExternal::TwoBarTrussExternal(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input,Echo),
   DOFS_(2),
   DOF_(DOFS_,0.0),
   Lambda_(0.0),
   Width_(Width),
   E1CachedValue_(2),
   E1DLoadCachedValue_(2),
   E2CachedValue_(2,2),
   EmptyV_(2,0.0),
   EmptyM_(2,2,0.0)
{
   LoadParameter_ = Load;
   for (int i=0;i<cachesize;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
}

void TwoBarTrussExternal::UpdateValues() const
{
   fstream in;
   fstream out;
   out.open("TwoBarTrussModelInput",ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << "TwoBarTrussModelInput" << " for write"
           << "\n";
      exit(-1);
   }
   out << setiosflags(ios::scientific) << setprecision(20);
   out << "use Math::Trig;\n\n"
       << "%Main = (FieldWidth => " << 30 << ",\n"
       << "         Precision  => " << 20 << ");\n\n"
       << "%Lattice = (TwoBarTruss => {Caching => No,\n"
       << "                            Theta => 65*atan(1,1)*4/180});\n\n"
       << "%TwoBarTrussModel = (DOF    => [" << setw(30) << DOF_[0] << "," << setw(30) << DOF_[1] << "],\n"
       << "                     Lambda => " << setw(30) << Lambda_ << ");\n\n";

   out.close();

   system("TwoBarTrussModel TwoBarTrussModelInput TwoBarTrussModelOutput >& /dev/null");

   in.open("TwoBarTrussModelOutput",ios::in);
   if (in.fail())
   {
      cerr << "Error: Unable to open file : " << "TwoBarTrussModelOutput" << "for read" << "\n";
      exit(-2);
   }
   in >> E0CachedValue_;
   in >> E1CachedValue_;
   in >> E1DLoadCachedValue_;
   in >> E2CachedValue_;

   for (int i=0;i<cachesize;++i) Cached_[i] = 1;

   in.close();
}

double TwoBarTrussExternal::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues();
      Cached_[0] = 1;
      CallCount_[0]++;
   }
   
   return E0CachedValue_;
}

Vector const& TwoBarTrussExternal::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues();
      Cached_[1] = 1;
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& TwoBarTrussExternal::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues();
      Cached_[2] = 1;
      CallCount_[2]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& TwoBarTrussExternal::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues();
      Cached_[3] = 1;
      CallCount_[3]++;
   }
   
   return E2CachedValue_;
}

void TwoBarTrussExternal::Print(ostream& out,PrintDetail const& flag,
                                PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   double mintestfunct;
   Matrix
      stiff(DOFS_,DOFS_);
   Vector str(DOFS_);
   Vector TestFunctVals(NumTestFunctions());
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   engy = E0();
   str = E1();
   stiff = E2();
   
   NoNegTestFunctions=TestFunctions(TestFunctVals,LHS);
   mintestfunct = TestFunctVals[0];
   for (int i=0;i<DOFS_;++i)
   {
      if (mintestfunct > TestFunctVals[i])
         mintestfunct = TestFunctVals[i];
   }
   
   switch (flag)
   {
      case PrintLong:
         out << "TwoBarTrussExternal:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "TwoBarTrussExternal:" << "\n" << "\n";
         }
         // passthrough to short
      case PrintShort:
         out << "Lambda: " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Stress:" << "\n" << setw(W) << str << "\n\n"
             << "Stiffness:" << setw(W) << stiff
             << "Eigenvalue Info:"  << "\n"<< setw(W) << TestFunctVals << "\n"
             << "Bifurcation Info:" << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda: " << setw(W) << Lambda_ << "\n"
                 << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                 << "Potential Value:" << setw(W) << engy << "\n";

            cout << "Stress:" << "\n" << setw(W) << str << "\n\n"
                 << "Stiffness:" << setw(W) << stiff
                 << "Eigenvalue Info:"  << "\n" << setw(W) << TestFunctVals <<"\n"
                 << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out,TwoBarTrussExternal& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
