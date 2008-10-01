#include "QC.h"

using namespace std;

extern "C" void get_qc_(int& mode,int& nfree,double* u,double& t,double& E,double* Eu,
                        double* Euu,double* Eut);

QC::~QC()
{
   cout << "QC Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

QC::QC(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input),
   Lambda_(0.0),
   Echo_(Echo),
   Width_(Width)
{
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash,"QC");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   DOF_.Resize(DOFS_,0.0);
   E1CachedValue_.Resize(DOFS_);
   E1DLoadCachedValue_.Resize(DOFS_);
   E2CachedValue_.Resize(DOFS_,DOFS_);
   EmptyV_.Resize(DOFS_,0.0);
   EmptyM_.Resize(DOFS_,DOFS_,0.0);

   LoadParameter_ = Load;
   for (int i=0;i<cachesize;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
}

void QC::UpdateValues(UpdateFlag flag) const
{
   if (NoStiffness==flag)
   {
      int mode=0;
      get_qc_(mode,DOFS_,&(DOF_[0]),Lambda_,E0CachedValue_,&(E1CachedValue_[0]),0,0);
      Cached_[0]=1;
      Cached_[1]=1;
   }
   else if (NeedStiffness==flag)
   {
      int mode=1;
      get_qc_(mode,DOFS_,&(DOF_[0]),Lambda_,E0CachedValue_,&(E1CachedValue_[0]),
              &(E2CachedValue_[0][0]),&(E1DLoadCachedValue_[0]));
      Cached_[0]=1;
      Cached_[1]=1;
      Cached_[2]=1;
      Cached_[3]=1;
   }
   else
   {
      cerr << "Error in QC::UpdateValues(), unknown UpdateFlag.\n";
      exit(-45);
   }
}

double QC::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues(NoStiffness);
      CallCount_[0]++;
   }
   
   return E0CachedValue_;
}

Vector const& QC::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& QC::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues(NeedStiffness);
      CallCount_[2]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& QC::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
      CallCount_[3]++;
   }
   
   return E2CachedValue_;
}

void QC::Print(ostream& out,PrintDetail const& flag)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   double mintestfunct;
   Vector TestFunctVals(DOFS_);
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   engy = E0();
   
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
         out << "QC:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "QC:" << "\n" << "\n";
         }
         // passthrough to short
      case PrintShort:
         out << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Bifurcation Info:" << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
                 << "Potential Value:" << setw(W) << engy << "\n";

            cout << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out,QC& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
