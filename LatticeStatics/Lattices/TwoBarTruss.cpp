#include "TwoBarTruss.h"
#include <cstdlib>
#include <cstring>

using namespace std;

TwoBarTruss::~TwoBarTruss()
{
   cout << "TwoBarTruss Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n"
        << "\tE3 calls - " << CallCount_[4] << "\n"
        << "\tE4 calls - " << CallCount_[5] << "\n";
}

TwoBarTruss::TwoBarTruss(PerlInput const& Input, int const& Echo, int const& Width) :
   Lattice(Input, Echo),
   DOFS_(2),
   DOF_(DOFS_, 0.0),
   Lambda_(0.0),
   Gamma_(0.0),
   Width_(Width),
   E1CachedValue_(2),
   E1DLoadCachedValue_(2),
   E2CachedValue_(2, 2),
   E3CachedValue_(4, 2),
   E4CachedValue_(4, 4),
   EmptyV_(2, 0.0),
   EmptyM_(2, 2, 0.0)
{
   LoadParameter_ = Load;
   for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }

   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash, "TwoBarTruss");
   const char* const caching = Input.getString(Hash, "Caching");
   if (!strcmp("Yes", caching))
   {
      Caching_ = 1;
   }
   else
   {
      Caching_ = 0;
   }
   Theta_ = Input.getDouble(Hash, "Theta");
   COSTheta_ = cos(Theta_);
   SINTheta_ = sin(Theta_);

   Gamma_ = Input.getDouble(Hash, "Gamma");

   if (Input.ParameterOK(Hash, "NumExtraTFs"))
   {
      NumExtraTFs_ = Input.getPosInt(Hash, "NumExtraTFs");
   }
   else
   {
      NumExtraTFs_ = 0;
      Input.usePosInt(0, Hash, "NumExtraTFs"); // Default Value
   }
   if (NumExtraTFs_ > 0)
   {
      if (Input.ParameterOK(Hash, "ExtraTFs"))
      {
         if (Input.getArrayLength(Hash, "ExtraTFs") == NumExtraTFs_)
         {
            ExtraTestFunctions_.Resize(NumExtraTFs_);
            Input.getVector(ExtraTestFunctions_, Hash, "ExtraTFs");
         }
         else
         {
            cerr << "Error: ArrayLength of " << Hash.Name
                 << "{ExtraTFs} is not equal to Lattice{NumExtraTFs}.\n";
            exit(-2);
         }
      }
      else
      {
         cerr << "Error: ExtraTFs not defined but Lattice{NumExtraTFs} = "
              << NumExtraTFs_ << ".\n";
         exit(-3);
      }
   }
   else
   {
      ExtraTestFunctions_.Resize(0);
   }

   Input.EndofInputSection();
}

double TwoBarTruss::E0() const
{
   if ((!Caching_) || (!Cached_[0]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     - 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);

      E0CachedValue_ = 0.5 * (Gamma_ * eps1_ * eps1_ + eps2_ * eps2_) - Lambda_ * DOF_[1];
      Cached_[0] = 1;
      CallCount_[0]++;
   }

   return E0CachedValue_;
}

Vector const& TwoBarTruss::E1() const
{
   if ((!Caching_) || (!Cached_[1]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     - 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps1u_ = DOF_[0] + COSTheta_;
      eps1v_ = DOF_[1] - SINTheta_;
      eps2u_ = DOF_[0] - COSTheta_;
      eps2v_ = DOF_[1] - SINTheta_;

      E1CachedValue_[0] = Gamma_ * eps1_ * eps1u_ + eps2_ * eps2u_;
      E1CachedValue_[1] = Gamma_ * eps1_ * eps1v_ + eps2_ * eps2v_ - Lambda_;

      Cached_[1] = 1;
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& TwoBarTruss::E1DLoad() const
{
   if ((!Caching_) || (!Cached_[2]))
   {
      E1DLoadCachedValue_[0] = 0.0;
      E1DLoadCachedValue_[1] = -1.0;

      Cached_[2] = 1;
      CallCount_[2]++;
   }

   return E1DLoadCachedValue_;
}

Matrix const& TwoBarTruss::E2() const
{
   if ((!Caching_) || (!Cached_[3]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     - 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps1u_ = DOF_[0] + COSTheta_;
      eps1v_ = DOF_[1] - SINTheta_;
      eps2u_ = DOF_[0] - COSTheta_;
      eps2v_ = DOF_[1] - SINTheta_;
      eps1uu_ = 1.0;
      eps1vv_ = 1.0;
      eps1uv_ = 0.0;
      eps2uu_ = 1.0;
      eps2vv_ = 1.0;
      eps2uv_ = 0.0;

      E2CachedValue_[0][0] = Gamma_ * (eps1u_ * eps1u_ + eps1_ * eps1uu_)
                             + (eps2u_ * eps2u_ + eps2_ * eps2uu_);
      E2CachedValue_[0][1] = Gamma_ * (eps1v_ * eps1u_ + eps1_ * eps1uv_)
                             + (eps2v_ * eps2u_ + eps2_ * eps2uv_);
      E2CachedValue_[1][0] = E2CachedValue_[0][1];
      E2CachedValue_[1][1] = Gamma_ * (eps1v_ * eps1v_ + eps1_ * eps1vv_)
                             + (eps2v_ * eps2v_ + eps2_ * eps2vv_);

      Cached_[3] = 1;
      CallCount_[3]++;
   }

   return E2CachedValue_;
}

Matrix const& TwoBarTruss::E3() const
{
   if ((!Caching_) || (!Cached_[4]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     - 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps1u_ = DOF_[0] + COSTheta_;
      eps1v_ = DOF_[1] - SINTheta_;
      eps2u_ = DOF_[0] - COSTheta_;
      eps2v_ = DOF_[1] - SINTheta_;
      eps1uu_ = 1.0;
      eps1vv_ = 1.0;
      eps1uv_ = 0.0;
      eps2uu_ = 1.0;
      eps2vv_ = 1.0;
      eps2uv_ = 0.0;

      // uuu
      E3CachedValue_[0][0] = Gamma_ * (3.0 * eps1uu_ * eps1u_)
                             + (3.0 * eps2uu_ * eps2u_);
      // uuv
      E3CachedValue_[0][1] = Gamma_ * (2.0 * eps1uv_ * eps1u_ + eps1v_ * eps1uu_)
                             + (2.0 * eps2uv_ * eps2u_ + eps2v_ * eps2uu_);
      // uvu
      E3CachedValue_[1][0] = E3CachedValue_[0][1];
      // uvv
      E3CachedValue_[1][1] = Gamma_ * (2.0 * eps1uv_ * eps1v_ + eps1u_ * eps1vv_)
                             + (2.0 * eps2uv_ * eps2v_ + eps2u_ * eps2vv_);
      // vuu
      E3CachedValue_[2][0] = E3CachedValue_[0][1];
      // vuv
      E3CachedValue_[2][1] = E3CachedValue_[1][1];
      // vvu
      E3CachedValue_[3][0] = E3CachedValue_[1][1];
      // vvv
      E3CachedValue_[3][1] = Gamma_ * (3.0 * eps1vv_ * eps1v_)
                             + (3.0 * eps2vv_ * eps2v_);

      Cached_[4] = 1;
      CallCount_[4]++;
   }

   return E3CachedValue_;
}

Matrix const& TwoBarTruss::E4() const
{
   if ((!Caching_) || (!Cached_[5]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     - 2.0 * COSTheta_ * DOF_[0] - 2.0 * SINTheta_ * DOF_[1]);
      eps1u_ = DOF_[0] + COSTheta_;
      eps1v_ = DOF_[1] - SINTheta_;
      eps2u_ = DOF_[0] - COSTheta_;
      eps2v_ = DOF_[1] - SINTheta_;
      eps1uu_ = 1.0;
      eps1vv_ = 1.0;
      eps1uv_ = 0.0;
      eps2uu_ = 1.0;
      eps2vv_ = 1.0;
      eps2uv_ = 0.0;

      // uuuu
      E4CachedValue_[0][0] = Gamma_ * (3.0 * eps1uu_ * eps1uu_)
                             + (3.0 * eps2uu_ * eps2uu_);
      // uuuv
      E4CachedValue_[0][1] = Gamma_ * (3.0 * eps1uu_ * eps1uv_)
                             + (3.0 * eps2uu_ * eps2uv_);
      // uuvu
      E4CachedValue_[0][2] = E4CachedValue_[0][1];
      // uuvv
      E4CachedValue_[0][3] = Gamma_ * (2.0 * eps1uv_ * eps1uv_ + eps1vv_ * eps1uu_)
                             + (2.0 * eps2uv_ * eps2uv_ + eps2vv_ * eps2uu_);
      // uvuu
      E4CachedValue_[1][0] = E4CachedValue_[0][1];
      // uvuv
      E4CachedValue_[1][1] = E4CachedValue_[0][3];
      // uvvu
      E4CachedValue_[1][2] = E4CachedValue_[0][3];
      // uvvv
      E4CachedValue_[1][3] = Gamma_ * (3.0 * eps1vv_ * eps1uv_)
                             + (3.0 * eps2vv_ * eps2uv_);
      // vuuu
      E4CachedValue_[2][0] = E4CachedValue_[0][2];
      // vuuv
      E4CachedValue_[2][1] = E4CachedValue_[1][2];
      // vuvu
      E4CachedValue_[2][2] = E4CachedValue_[1][2];
      // vuvv
      E4CachedValue_[2][3] = E4CachedValue_[1][3];
      // vvuu
      E4CachedValue_[3][0] = E4CachedValue_[0][3];
      // vvuv
      E4CachedValue_[3][1] = E4CachedValue_[1][3];
      // vvvu
      E4CachedValue_[3][2] = E4CachedValue_[2][3];
      // vvvv
      E4CachedValue_[3][3] = Gamma_ * (3.0 * eps1vv_ * eps1vv_)
                             + (3.0 * eps2vv_ * eps2vv_);

      Cached_[5] = 1;
      CallCount_[5]++;
   }

   return E4CachedValue_;
}

void TwoBarTruss::ExtraTestFunctions(Vector& TF) const
{
   for (int i = 0; i < NumExtraTFs_; ++i)
   {
      TF[i] = (ExtraTestFunctions_[i] - Lambda());
   }
}

void TwoBarTruss::Print(ostream& out, PrintDetail const& flag,
                        PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions = 0;
   double engy;
   double mintestfunct;
   Matrix
   stiff(DOFS_, DOFS_);
   Vector str(DOFS_);
   Vector TestFunctVals(NumTestFunctions());

   W = out.width();

   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }

   engy = E0();
   str = E1();
   stiff = E2();

   TestFunctions(TestFunctVals, LHS);
   mintestfunct = TestFunctVals[0];
   for (int i = 0; i < NumTestFunctions(); ++i)
   {
      if ((TestFunctVals[i] < 0.0) && (i < DOFS_))
      {
         ++NoNegTestFunctions;
      }
      if (mintestfunct > TestFunctVals[i])
      {
         mintestfunct = TestFunctVals[i];
      }
   }

   switch (flag)
   {
      case PrintLong:
         out << "TwoBarTruss:" << "\n" << "\n";
         out << "Theta:" << setw(W) << Theta_ << "\n";
         out << "Gamma:" << setw(W) << Gamma_ << "\n";

         if (Echo_)
         {
            cout << "TwoBarTruss:" << "\n" << "\n";
            cout << "Theta:" << setw(W) << Theta_ << "\n";
            cout << "Gamma:" << setw(W) << Gamma_ << "\n";
         }
      // passthrough to short
      case PrintShort:
         out << "Lambda: " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Stress:" << "\n" << setw(W) << str << "\n\n"
             << "Stiffness:" << setw(W) << stiff
             << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
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
                 << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
                 << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out, TwoBarTruss& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}
