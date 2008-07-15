#include "TwoBarTruss.h"

using namespace std;

TwoBarTruss::~TwoBarTruss()
{
   cout << "TwoBarTruss Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

TwoBarTruss::TwoBarTruss(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input),
   DOFS_(2),
   DOF_(DOFS_,0.0),
   Lambda_(0.0),
   Echo_(Echo),
   Width_(Width),
   E1CachedValue_(1,2),
   E1DLoadCachedValue_(1,2),
   E2CachedValue_(2,2)
{
   LoadParameter_ = Load;
   for (int i=0;i<4;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
   
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash,"TwoBarTruss");
   const char* const caching = Input.getString(Hash,"Caching");
   if (!strcmp("Yes",caching))
   {
      Caching_ = 1;
   }
   else
   {
      Caching_ = 0;
   }
   Theta_ = Input.getDouble(Hash,"Theta");
   Input.EndofInputSection();
   COSTheta_ = cos(Theta_);
   SINTheta_ = sin(Theta_);
}

double TwoBarTruss::E0() const
{
   if ((!Caching_) || (!Cached_[0]))
   {
      E0CachedValue_ = 0.125*((DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] + 2.0*COSTheta_*DOF_[0]
                              - 2.0*SINTheta_*DOF_[1])
                             *(DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] + 2.0*COSTheta_*DOF_[0]
                               - 2.0*SINTheta_*DOF_[1])
                             +(DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] - 2.0*COSTheta_*DOF_[0]
                               - 2.0*SINTheta_*DOF_[1])
                             *(DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] + 2.0*COSTheta_*DOF_[0]
                               - 2.0*SINTheta_*DOF_[1]))
         - Lambda_*DOF_[1];
      Cached_[0] = 1;
      CallCount_[0]++;
   }
   
   return E0CachedValue_;
}

Matrix const& TwoBarTruss::E1() const
{
   if ((!Caching_) || (!Cached_[1]))
   {
      E1CachedValue_[0][0] = DOF_[0]*(DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] - 2.0*SINTheta_*DOF_[1]
                                     + 2.0*COSTheta_*COSTheta_);
      E1CachedValue_[0][1] = DOF_[1]*DOF_[1]*DOF_[1] + DOF_[0]*DOF_[0]*DOF_[1]
         - 3.0*SINTheta_*DOF_[1]*DOF_[1]
         - SINTheta_*DOF_[0]*DOF_[0] + 2.0*SINTheta_*SINTheta_*DOF_[1] - Lambda_;

      Cached_[1] = 1;
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Matrix const& TwoBarTruss::E1DLoad() const
{
   if ((!Caching_) || (!Cached_[2]))
   {
      E1DLoadCachedValue_[0][0] = 0.0;
      E1DLoadCachedValue_[0][1] = -1.0;

      Cached_[2] = 1;
      CallCount_[2]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& TwoBarTruss::E2() const
{
   if ((!Caching_) || (!Cached_[3]))
   {
      E2CachedValue_[0][0] = 3.0*DOF_[0]*DOF_[0] + DOF_[1]*DOF_[1] - 2.0*SINTheta_*DOF_[1]
         + 2.0*COSTheta_*COSTheta_;
      E2CachedValue_[0][1] = E2CachedValue_[1][0]
         = 2.0*DOF_[0]*DOF_[1] - 2.0*SINTheta_*DOF_[0];
      E2CachedValue_[1][1] = 3.0*DOF_[1]*DOF_[1] + DOF_[0]*DOF_[0] - 6.0*SINTheta_*DOF_[1]
         + 2.0*SINTheta_*SINTheta_;

      Cached_[3] = 1;
      CallCount_[3]++;
   }
   
   return E2CachedValue_;
}

void TwoBarTruss::Print(ostream& out,PrintDetail const& flag)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   Matrix
      str(1,DOFS_),
      stiff(DOFS_,DOFS_);
   Vector TestFunctVals(DOFS_);
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   engy = E0();
   str = E1();
   stiff = E2();
   
   NoNegTestFunctions=TestFunctions(TestFunctVals,LHS);
   
   switch (flag)
   {
      case PrintLong:
         out << "TwoBarTruss:" << "\n" << "\n";
         out << "Theta:" << setw(W) << Theta_ << "\n";

         if (Echo_)
         {
            cout << "TwoBarTruss:" << "\n" << "\n";
            cout << "Theta:" << setw(W) << Theta_ << "\n";
         }
         // passthrough to short
      case PrintShort:
         out << "Lambda: " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Stress:" << setw(W) << str
             << "Stiffness:" << setw(W) << stiff
             << "Eigenvalue Info:"  << "\n"<< setw(W) << TestFunctVals << "\n"
             << "Bifurcation Info:" << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda: " << setw(W) << Lambda_ << "\n"
                 << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                 << "Potential Value:" << setw(W) << engy << "\n";

            cout << "Stress:" << setw(W) << str
                 << "Stiffness:" << setw(W) << stiff
                 << "Eigenvalue Info:"  << "\n" << setw(W) << TestFunctVals <<"\n"
                 << "Bifurcation Info:" << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out,TwoBarTruss& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
