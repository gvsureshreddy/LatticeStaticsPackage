#include "FourBarTruss.h"

// Author: Karthikreddy Ginnavaram

using namespace std;

FourBarTruss::~FourBarTruss()
{
   cout << "FourBarTruss Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n"
        << "\tE3 calls - " << CallCount_[4] << "\n"
        << "\tE4 calls - " << CallCount_[5] << "\n";
}

FourBarTruss::FourBarTruss(PerlInput const& Input, int const& Echo, int const& Width) :
   Lattice(Input, Echo),
   DOFS_(3),
   DOF_(DOFS_, 0.0),
   Lambda_(0.0),
   Gamma_(0.0),
   Width_(Width),
   E1CachedValue_(3),
   E1DLoadCachedValue_(3),
   E2CachedValue_(3, 3),
   E3CachedValue_(9, 3),
   E4CachedValue_(9, 9),
   EmptyV_(3, 0.0),
   EmptyM_(3, 3, 0.0)
{
   LoadParameter_ = Load;
   for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash, "FourBarTruss");
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
   Psi_ = Input.getDouble(Hash, "Psi");
   COSPsi_ = cos(Psi_);
   SINPsi_ = sin(Psi_);
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

double FourBarTruss::E0() const
{
   if ((!Caching_) || (!Cached_[0]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps3_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps4_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);


      E0CachedValue_ = 0.5 * (Gamma_ * eps1_ * eps1_
                              + Gamma_ * eps2_ * eps2_ + eps3_ * eps3_
                              + eps4_ * eps4_) - Lambda_ * DOF_[2];
      Cached_[0] = 1;
      CallCount_[0]++;
   }

   return E0CachedValue_;
}


Vector const& FourBarTruss::E1() const
{
   if ((!Caching_) || (!Cached_[1]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps3_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps4_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);

      eps1u_ = DOF_[0] + COSTheta_ * SINPsi_;
      eps1v_ = DOF_[1] - SINTheta_ * SINPsi_;
      eps1w_ = DOF_[2] - COSPsi_;
      eps2u_ = DOF_[0] + SINTheta_ * SINPsi_;
      eps2v_ = DOF_[1] + COSTheta_ * SINPsi_;
      eps2w_ = DOF_[2] - COSPsi_;
      eps3u_ = DOF_[0] - COSTheta_ * SINPsi_;
      eps3v_ = DOF_[1] + SINTheta_ * SINPsi_;
      eps3w_ = DOF_[2] - COSPsi_;
      eps4u_ = DOF_[0] - SINTheta_ * SINPsi_;
      eps4v_ = DOF_[1] - COSTheta_ * SINPsi_;
      eps4w_ = DOF_[2] - COSPsi_;

      E1CachedValue_[0] = Gamma_ * eps1_ * eps1u_
                          + Gamma_ * eps2_ * eps2u_ + eps3_ * eps3u_
                          + eps4_ * eps4u_;
      E1CachedValue_[1] = Gamma_ * eps1_ * eps1v_
                          + Gamma_ * eps2_ * eps2v_ + eps3_ * eps3v_
                          + eps4_ * eps4v_;
      E1CachedValue_[2] = Gamma_ * eps1_ * eps1w_
                          + Gamma_ * eps2_ * eps2w_ + eps3_ * eps3w_
                          + eps4_ * eps4w_ - Lambda_;


      Cached_[1] = 1;
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& FourBarTruss::E1DLoad() const
{
   if ((!Caching_) || (!Cached_[2]))
   {
      E1DLoadCachedValue_[0] = 0.0;
      E1DLoadCachedValue_[1] = 0.0;
      E1DLoadCachedValue_[2] = -1.0;
      Cached_[2] = 1;
      CallCount_[2]++;
   }

   return E1DLoadCachedValue_;
}

Matrix const& FourBarTruss::E2() const
{
   if ((!Caching_) || (!Cached_[3]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps3_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps4_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);

      eps1u_ = DOF_[0] + COSTheta_ * SINPsi_;
      eps1v_ = DOF_[1] - SINTheta_ * SINPsi_;
      eps1w_ = DOF_[2] - COSPsi_;
      eps2u_ = DOF_[0] + SINTheta_ * SINPsi_;
      eps2v_ = DOF_[1] + COSTheta_ * SINPsi_;
      eps2w_ = DOF_[2] - COSPsi_;
      eps3u_ = DOF_[0] - COSTheta_ * SINPsi_;
      eps3v_ = DOF_[1] + SINTheta_ * SINPsi_;
      eps3w_ = DOF_[2] - COSPsi_;
      eps4u_ = DOF_[0] - SINTheta_ * SINPsi_;
      eps4v_ = DOF_[1] - COSTheta_ * SINPsi_;
      eps4w_ = DOF_[2] - COSPsi_;

      eps1uu_ = 1.0;
      eps1vv_ = 1.0;
      eps1ww_ = 1.0;
      eps1uv_ = 0.0;
      eps1vw_ = 0.0;
      eps1uw_ = 0.0;
      eps1vu_ = 0.0;
      eps1wv_ = 0.0;
      eps1wu_ = 0.0;

      eps2uu_ = 1.0;
      eps2vv_ = 1.0;
      eps2ww_ = 1.0;
      eps2uv_ = 0.0;
      eps2uw_ = 0.0;
      eps2vw_ = 0.0;
      eps2vu_ = 0.0;
      eps2wv_ = 0.0;
      eps2wu_ = 0.0;
      eps3uu_ = 1.0;
      eps3vv_ = 1.0;
      eps3ww_ = 1.0;
      eps3uv_ = 0.0;
      eps3uw_ = 0.0;
      eps3vw_ = 0.0;
      eps3vu_ = 0.0;
      eps3wv_ = 0.0;
      eps3wu_ = 0.0;
      eps4uu_ = 1.0;
      eps4vv_ = 1.0;
      eps4ww_ = 1.0;
      eps4uv_ = 0.0;
      eps4uw_ = 0.0;
      eps4vw_ = 0.0;
      eps4vu_ = 0.0;
      eps4wv_ = 0.0;
      eps4wu_ = 0.0;



      E2CachedValue_[0][0] = Gamma_ * (eps1u_ * eps1u_
                                       + eps1_ * eps1uu_ + eps2u_ * eps2u_
                                       + eps2_ * eps2uu_)
                             + (eps3u_ * eps3u_ + eps3_ * eps3uu_)
                             + (eps3u_ * eps3u_ + eps3_ * eps3uu_);
      E2CachedValue_[0][1] = Gamma_ * (eps1_ * eps1uv_
                                       + eps1v_ * eps1u_ + eps2v_ * eps2u_
                                       + eps2_ * eps2uv_)
                             + (eps3v_ * eps3u_ + eps3_ * eps3uv_)
                             + (eps4v_ * eps4u_ + eps4_ * eps4uv_);

      E2CachedValue_[0][2] = Gamma_ * (eps1_ * eps1uw_
                                       + eps1w_ * eps1u_ + eps2w_ * eps2u_
                                       + eps2_ * eps2uw_)
                             + (eps3w_ * eps3u_ + eps3_ * eps3uw_)
                             + (eps4w_ * eps4u_ + eps4_ * eps4uw_);


      E2CachedValue_[1][0] = E2CachedValue_[0][1];
      E2CachedValue_[1][1] = Gamma_ * (eps1v_ * eps1v_
                                       + eps1_ * eps1vv_ + eps2v_ * eps2v_
                                       + eps2_ * eps2vv_)
                             + (eps3v_ * eps3v_ + eps3_ * eps3vv_)
                             + (eps4v_ * eps4v_ + eps4_ * eps4vv_);
      E2CachedValue_[1][2] = Gamma_ * (eps1w_ * eps1v_
                                       + eps1_ * eps1vw_ + eps2w_ * eps2v_
                                       + eps2_ * eps2vw_)
                             + (eps3_ * eps3vw_ +
                                eps3w_ * eps3v_)
                             + (eps4w_ * eps4v_ + eps4_ * eps4vw_);
      E2CachedValue_[2][0] = E2CachedValue_[0][2];
      E2CachedValue_[2][1] = E2CachedValue_[1][2];
      E2CachedValue_[2][2] = Gamma_ * (eps1w_ * eps1w_
                                       + eps1_ * eps1ww_ + eps2w_ * eps2w_
                                       + eps2_ * eps2ww_)
                             + (eps3w_ * eps3w_ + eps3_ * eps3ww_)
                             + (eps4w_ * eps4w_ + eps4_ * eps4ww_);






      Cached_[3] = 1;
      CallCount_[3]++;
   }

   return E2CachedValue_;
}

Matrix const& FourBarTruss::E3() const
{
   eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                  + DOF_[2] * DOF_[2] + 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                  - 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
   eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                  + DOF_[2] * DOF_[2] + 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                  + 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
   eps3_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                  + DOF_[2] * DOF_[2] - 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                  + 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
   eps4_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                  + DOF_[2] * DOF_[2] - 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                  - 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);

   eps1u_ = DOF_[0] + COSTheta_ * SINPsi_;
   eps1v_ = DOF_[1] - SINTheta_ * SINPsi_;
   eps1w_ = DOF_[2] - COSPsi_;
   eps2u_ = DOF_[0] + SINTheta_ * SINPsi_;
   eps2v_ = DOF_[1] + COSTheta_ * SINPsi_;
   eps2w_ = DOF_[2] - COSPsi_;
   eps3u_ = DOF_[0] - COSTheta_ * SINPsi_;
   eps3v_ = DOF_[1] + SINTheta_ * SINPsi_;
   eps3w_ = DOF_[2] - COSPsi_;
   eps4u_ = DOF_[0] - SINTheta_ * SINPsi_;
   eps4v_ = DOF_[1] - COSTheta_ * SINPsi_;
   eps4w_ = DOF_[2] - COSPsi_;


   eps1uu_ = 1.0;
   eps1vv_ = 1.0;
   eps1ww_ = 1.0;
   eps1uv_ = 0.0;
   eps1vw_ = 0.0;
   eps1uw_ = 0.0;
   eps1vu_ = 0.0;
   eps1wv_ = 0.0;
   eps1wu_ = 0.0;
   eps2uu_ = 1.0;
   eps2vv_ = 1.0;
   eps2ww_ = 1.0;   // InitializeLat
   eps2uv_ = 0.0;
   eps2uw_ = 0.0;
   eps2vw_ = 0.0;
   eps2vu_ = 0.0;
   eps2wv_ = 0.0;
   eps2wu_ = 0.0;
   eps3uu_ = 1.0;
   eps3vv_ = 1.0;
   eps3ww_ = 1.0;
   eps3uv_ = 0.0;
   eps3uw_ = 0.0;
   eps3vw_ = 0.0;
   eps3vu_ = 0.0;
   eps3wv_ = 0.0;
   eps3wu_ = 0.0;
   eps4uu_ = 1.0;
   eps4vv_ = 1.0;
   eps4ww_ = 1.0;
   eps4uv_ = 0.0;
   eps4uw_ = 0.0;
   eps4vw_ = 0.0;
   eps4vu_ = 0.0;
   eps4wv_ = 0.0;
   eps4wu_ = 0.0;

   if ((!Caching_) || (!Cached_[4]))
   {
      // uuu
      E3CachedValue_[0][0] = Gamma_ * (3.0 * eps1uu_ * eps1u_
                                       + 3.0 * eps2uu_ * eps2u_)
                             + (3.0 * eps3uu_ * eps3u_)
                             + (3.0 * eps4uu_ * eps4u_);
      // uuv
      E3CachedValue_[0][1] = Gamma_ * (2.0 * eps1uv_ * eps1u_
                                       + eps1v_ * eps1uu_
                                       + 2.0 * eps2uv_ * eps2u_
                                       + eps2uu_ * eps2v_)
                             + 2.0 * eps3uv_ * eps3u_
                             + eps3v_ * eps3uu_
                             + 2.0 * eps4uv_ * eps4u_
                             + eps4v_ * eps4uu_;
      // uuw
      E3CachedValue_[0][2] = Gamma_ * (2.0 * eps1uw_ * eps1u_
                                       + eps1w_ * eps1uu_
                                       + 2.0 * eps2uw_ * eps2u_
                                       + eps2uu_ * eps2w_)
                             + eps3uu_ * eps3w_
                             + 2.0 * eps3u_ * eps3uw_
                             + 2.0 * eps4uw_ * eps4u_
                             + eps4uu_ * eps4w_;


      // uvu
      E3CachedValue_[1][0] = E3CachedValue_[0][1];
      // uvv
      E3CachedValue_[1][1] = Gamma_ * (2.0 * eps1uv_ * eps1v_
                                       + eps1u_ * eps1vv_
                                       + eps2vv_ * eps2u_
                                       + 2.0 * eps2v_ * eps2uv_)
                             + eps3vv_ * eps3u_
                             + 2.0 * eps3uv_ * eps3v_
                             + eps4vv_ * eps4u_
                             + 2.0 * eps4uv_ * eps4v_;
      // uvw
      E3CachedValue_[1][2] = Gamma_ * (eps1uv_ * eps1w_
                                       + eps1u_ * eps1vw_
                                       + eps1v_ * eps1uw_
                                       + eps2vw_ * eps2u_
                                       + eps2uw_ * eps2v_
                                       + eps2uv_ * eps2w_)
                             + eps3uw_ * eps3v_
                             + eps3vw_ * eps3u_
                             + eps3uv_ * eps3w_
                             + eps4uw_ * eps4v_
                             + eps4vw_ * eps4u_
                             + eps4uv_ * eps4w_;

      // uwu
      E3CachedValue_[2][0] = E3CachedValue_[0][2];
      // uwv
      E3CachedValue_[2][1] = E3CachedValue_[1][2];
      // uww
      E3CachedValue_[2][2] = Gamma_ * (2.0 * eps1uw_ * eps1w_
                                       + eps1ww_ * eps1u_
                                       + 2.0 * eps2w_ * eps2uw_
                                       + eps2u_ * eps2ww_)
                             + eps3ww_ * eps3u_
                             + 2.0 * eps3uw_ * eps3w_
                             + 2.0 * eps4uw_ * eps4w_
                             + eps4ww_ * eps4u_;

      // vuu
      E3CachedValue_[3][0] = E3CachedValue_[0][1];
      // vuv
      E3CachedValue_[3][1] = E3CachedValue_[1][1];
      // vuw
      E3CachedValue_[3][2] = E3CachedValue_[1][2];
      // vvu
      E3CachedValue_[4][0] = E3CachedValue_[1][1];
      // vvv
      E3CachedValue_[4][1] = Gamma_ * (3.0 * eps1vv_ * eps1v_
                                       + 3.0 * eps2vv_ * eps2v_)
                             + (3.0 * eps3vv_ * eps3v_)
                             + (3.0 * eps4vv_ * eps4v_);

      // vvw
      E3CachedValue_[4][2] = Gamma_ * (2.0 * eps1vw_ * eps1v_
                                       + eps1vv_ * eps1w_
                                       + 2.0 * eps2v_ * eps2vw_
                                       + eps2w_ * eps2vv_)
                             + 2.0 * eps3vw_ * eps3v_
                             + eps3vv_ * eps3w_
                             + 2.0 * eps4vw_ * eps4v_
                             + eps4vv_ * eps4w_;

      // vwu
      E3CachedValue_[5][0] = E3CachedValue_[1][2];
      // vwv
      E3CachedValue_[5][1] = E3CachedValue_[4][2];
      // vww
      E3CachedValue_[5][2] = Gamma_ * (eps1ww_ * eps1v_
                                       + 2.0 * eps1vw_ * eps1w_
                                       + eps2v_ * eps2ww_
                                       + 2.0 * eps2w_ * eps2vw_)
                             + 2.0 * eps3vw_ * eps3w_
                             + eps3ww_ * eps3v_
                             + 2.0 * eps4vw_ * eps4w_
                             + eps4ww_ * eps4v_;

      // wuu
      E3CachedValue_[6][0] = E3CachedValue_[0][2];
      // wuv
      E3CachedValue_[6][1] = E3CachedValue_[1][2];
      // wuw
      E3CachedValue_[6][2] = E3CachedValue_[2][2];

      // wvu
      E3CachedValue_[7][0] = E3CachedValue_[1][2];
      // wvv
      E3CachedValue_[7][1] = E3CachedValue_[4][2];
      // wvw
      E3CachedValue_[7][2] = E3CachedValue_[5][2];

      // wwu
      E3CachedValue_[8][0] = E3CachedValue_[2][2];
      // wwv
      E3CachedValue_[8][1] = E3CachedValue_[5][2];
      // www
      E3CachedValue_[8][2] = Gamma_ * (3.0 * eps1ww_ * eps1w_
                                       + 3.0 * eps2ww_ * eps2w_)
                             + (3.0 * eps3ww_ * eps3w_)
                             + (3.0 * eps4ww_ * eps4w_);
   }

   return E3CachedValue_;
}

Matrix const& FourBarTruss::E4() const
{
   if ((!Caching_) || (!Cached_[5]))
   {
      eps1_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps2_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] + 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps3_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * COSTheta_ * SINPsi_ * DOF_[0]
                     + 2.0 * SINTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps4_ = 0.5 * (DOF_[0] * DOF_[0] + DOF_[1] * DOF_[1]
                     + DOF_[2] * DOF_[2] - 2.0 * SINTheta_ * SINPsi_ * DOF_[0]
                     - 2.0 * COSTheta_ * SINPsi_ * DOF_[1] - 2 * COSPsi_ * DOF_[2]);
      eps1u_ = DOF_[0] + COSTheta_ * SINPsi_;
      eps1v_ = DOF_[1] - SINTheta_ * SINPsi_;
      eps1w_ = DOF_[2] - COSPsi_;
      eps2u_ = DOF_[0] + SINTheta_ * SINPsi_;
      eps2v_ = DOF_[1] + COSTheta_ * SINPsi_;
      eps2w_ = DOF_[2] - COSPsi_;
      eps3u_ = DOF_[0] - COSTheta_ * SINPsi_;
      eps3v_ = DOF_[1] + SINTheta_ * SINPsi_;
      eps3w_ = DOF_[2] - COSPsi_;
      eps4u_ = DOF_[0] - SINTheta_ * SINPsi_;
      eps4v_ = DOF_[1] - COSTheta_ * SINPsi_;
      eps4w_ = DOF_[2] - COSPsi_;
      eps1uu_ = 1.0;
      eps1vv_ = 1.0;
      eps1ww_ = 1.0;
      eps1uv_ = 0.0;
      eps1vw_ = 0.0;
      eps1uw_ = 0.0;
      eps1vu_ = 0.0;
      eps1wv_ = 0.0;
      eps1wu_ = 0.0;
      eps2uu_ = 1.0;
      eps2vv_ = 1.0;
      eps2ww_ = 1.0;
      eps2uv_ = 0.0;
      eps2uw_ = 0.0;
      eps2vw_ = 0.0;
      eps2vu_ = 0.0;
      eps2wv_ = 0.0;
      eps2wu_ = 0.0;
      eps3uu_ = 1.0;
      eps3vv_ = 1.0;
      eps3ww_ = 1.0;
      eps3uv_ = 0.0;
      eps3uw_ = 0.0;
      eps3vw_ = 0.0;
      eps3vu_ = 0.0;
      eps3wv_ = 0.0;
      eps3wu_ = 0.0;
      eps4uu_ = 1.0;
      eps4vv_ = 1.0;
      eps4ww_ = 1.0;
      eps4uv_ = 0.0;
      eps4uw_ = 0.0;
      eps4vw_ = 0.0;
      eps4vu_ = 0.0;
      eps4wv_ = 0.0;
      eps4wu_ = 0.0;

      // uuuu
      E4CachedValue_[0][0] = Gamma_ * (3.0 * eps1uu_ * eps1uu_
                                       + 3.0 * eps2uu_ * eps2uu_)
                             + 3.0 * eps3uu_ * eps3uu_
                             + 3.0 * eps4uu_ * eps4uu_;
      // uuuv
      E4CachedValue_[0][1] = Gamma_ * (3.0 * eps1uu_ * eps1uv_
                                       + 3.0 * eps2uu_ * eps2uv_)
                             + 3.0 * eps3uv_ * eps3uu_
                             + 3.0 * eps4uv_ * eps4uu_;
      // uuuw
      E4CachedValue_[0][2] = Gamma_ * (3.0 * eps1uu_ * eps1uw_
                                       + 3.0 * eps2uu_ * eps2uw_)
                             + 3.0 * eps3uu_ * eps3uw_
                             + 3.0 * eps4uu_ * eps4uw_;
      // uuvu
      E4CachedValue_[0][3] = E4CachedValue_[0][1];
      // uuvv
      E4CachedValue_[0][4] = Gamma_ * (2.0 * eps1uv_ * eps1uv_
                                       + eps1vv_ * eps1uu_
                                       + 2.0 * eps2uv_ * eps2uv_
                                       + eps2vv_ * eps2uu_)
                             + 2.0 * eps3uv_ * eps3uv_
                             + eps3vv_ * eps3uu_
                             + 2.0 * eps4uv_ * eps4uv_
                             + eps4vv_ * eps4uu_;
      // uuvw
      E4CachedValue_[0][5] = Gamma_ * (2.0 * eps1uv_ * eps1uv_
                                       + eps1vw_ * eps1uu_
                                       + 2.0 * eps2uu_ * eps2vv_
                                       + eps2vw_ * eps2uu_)
                             + eps3vw_ * eps3uu_
                             + 2.0 * eps3uv_ * eps3uw_
                             + 2.0 * eps4uw_ * eps4uv_
                             + eps4uu_ * eps4vw_;
      // uuwu
      E4CachedValue_[0][6] = E4CachedValue_[0][2];
      // uuwv
      E4CachedValue_[0][7] = E4CachedValue_[0][5];
      // uuww
      E4CachedValue_[0][8] = Gamma_ * (2.0 * eps1uw_ * eps1uw_
                                       + eps1ww_ * eps1uu_
                                       + 2.0 * eps2uw_ * eps2uw_
                                       + eps2ww_ * eps2uu_)
                             + eps3ww_ * eps3uu_
                             + 2 * eps3uw_ * eps3uw_
                             + 2.0 * eps4uw_ * eps4uw_
                             + eps4ww_ * eps4uu_;

      // uvuu
      E4CachedValue_[1][0] = E4CachedValue_[0][1];
      // uvuv
      E4CachedValue_[1][1] = E4CachedValue_[0][4];
      // uvuw
      E4CachedValue_[1][2] = E4CachedValue_[0][5];
      // uvvu
      E4CachedValue_[1][3] = E4CachedValue_[0][4];
      // uvvv
      E4CachedValue_[1][4] = Gamma_ * (3.0 * eps1vv_ * eps1uv_
                                       + 3.0 * eps2vv_ * eps2uv_)
                             + 3.0 * eps3vv_ * eps3uv_
                             + 3.0 * eps4vv_ * eps4uv_;
      // uvvw
      E4CachedValue_[1][5] = Gamma_ * (2.0 * eps1uv_ * eps1vw_
                                       + eps1vv_ * eps1uw_
                                       + 2.0 * eps2uv_ * eps2vw_
                                       + eps2vv_ * eps2uw_)
                             + 2.0 * eps3vw_ * eps3uv_
                             + eps3vv_ * eps3uw_
                             + 2.0 * eps4vw_ * eps4uv_
                             + eps4vv_ * eps4uw_;
      // uvwu

      E4CachedValue_[1][6] = Gamma_ * (2.0 * eps1uw_ * eps1uv_
                                       + eps1uu_ * eps1wv_
                                       + 2.0 * eps2uv_ * eps2uw_
                                       + eps2uu_ * eps2vw_)
                             + 2.0 * eps3vw_ * eps3uu_
                             + eps3vw_ * eps3uu_
                             + 2.0 * eps4vw_ * eps4uu_
                             + eps4vw_ * eps4uu_;
      // uvwv
      E4CachedValue_[1][7] = E4CachedValue_[1][5];
      // uvww
      E4CachedValue_[1][8] = Gamma_ * (eps1ww_ * eps1uv_
                                       + 2.0 * eps1uw_ * eps1vw_
                                       + eps2uv_ * eps2ww_
                                       + 2.0 * eps2vw_ * eps2uw_)
                             + 2.0 * eps3vw_ * eps3uw_
                             + eps3ww_ * eps3uv_
                             + 2.0 * eps4vw_ * eps4uw_
                             + eps4ww_ * eps4uv_;

      // uwuu
      E4CachedValue_[2][0] = E4CachedValue_[0][2];
      // uwuv
      E4CachedValue_[2][1] = E4CachedValue_[1][2];
      // uwuw
      E4CachedValue_[2][2] = E4CachedValue_[0][8];
      // uwvu
      E4CachedValue_[2][3] = E4CachedValue_[0][5];
      // uwvv
      E4CachedValue_[2][4] = E4CachedValue_[1][5];
      // uwvw
      E4CachedValue_[2][5] = E4CachedValue_[1][8];
      // uwwu
      E4CachedValue_[2][6] = E4CachedValue_[0][8];
      // uwwv
      E4CachedValue_[2][7] = E4CachedValue_[1][8];
      // uwww
      E4CachedValue_[2][8] = Gamma_ * (3.0 * eps1ww_ * eps1uw_
                                       + 3.0 * eps2ww_ * eps2uw_)
                             + 3.0 * eps3ww_ * eps3uw_
                             + 3.0 * eps4ww_ * eps4uw_;

      // vuuu
      E4CachedValue_[3][0] = E4CachedValue_[0][1];
      // vuuv
      E4CachedValue_[3][1] = E4CachedValue_[0][4];
      // vuuw
      E4CachedValue_[3][2] = E4CachedValue_[0][5];
      // vuvu
      E4CachedValue_[3][3] = E4CachedValue_[0][4];
      // vuvv
      E4CachedValue_[3][4] = E4CachedValue_[1][4];
      // vuvw
      E4CachedValue_[3][5] = E4CachedValue_[1][5];
      // vuwu
      E4CachedValue_[3][6] = E4CachedValue_[0][5];
      // vuwv
      E4CachedValue_[3][7] = E4CachedValue_[1][5];
      // vuww
      E4CachedValue_[3][8] = E4CachedValue_[1][8];

      // vvuu
      E4CachedValue_[4][0] = E4CachedValue_[0][4];
      // vvuv
      E4CachedValue_[4][1] = E4CachedValue_[1][4];
      // vvuw
      E4CachedValue_[4][2] = E4CachedValue_[1][5];
      // vvvu
      E4CachedValue_[4][3] = E4CachedValue_[1][4];
      // vvvv
      E4CachedValue_[4][4] = Gamma_ * (3.0 * eps1vv_ * eps1vv_
                                       + 3.0 * eps2vv_ * eps2vv_)
                             + 3.0 * eps3vv_ * eps3vv_
                             + 3.0 * eps4vv_ * eps4vv_;

      // vvvw
      E4CachedValue_[4][5] = Gamma_ * (3.0 * eps1vv_ * eps1vw_
                                       + 3.0 * eps2vv_ * eps2vw_)
                             + 3.0 * eps3vv_ * eps3vw_
                             + 3.0 * eps4vv_ * eps4vw_;
      // vvwu
      E4CachedValue_[4][6] = E4CachedValue_[1][5];
      // vvwv
      E4CachedValue_[4][7] = E4CachedValue_[4][5];
      // vvww
      E4CachedValue_[4][8] = Gamma_ * (2.0 * eps1vw_ * eps1vw_
                                       + eps1ww_ * eps1vv_
                                       + 2.0 * eps2vw_ * eps2vw_
                                       + eps2ww_ * eps2vv_)
                             + eps3ww_ * eps3vv_
                             + 2 * eps3vw_ * eps3vw_
                             + 2.0 * eps4vw_ * eps4vw_
                             + eps4ww_ * eps4vv_;

      // vwuu
      E4CachedValue_[5][0] = E4CachedValue_[0][7];
      // vwuv
      E4CachedValue_[5][1] = E4CachedValue_[1][5];
      // vwuw
      E4CachedValue_[5][2] = E4CachedValue_[1][8];
      // vwvu
      E4CachedValue_[5][3] = E4CachedValue_[1][5];
      // vwvv
      E4CachedValue_[5][4] = E4CachedValue_[4][5];
      // vwvw
      E4CachedValue_[5][5] = E4CachedValue_[4][8];
      // vwwu
      E4CachedValue_[5][6] = E4CachedValue_[0][5];
      // vwwv
      E4CachedValue_[5][7] = E4CachedValue_[4][8];
      // vwww
      E4CachedValue_[5][8] = Gamma_ * (3.0 * eps1ww_ * eps1vw_
                                       + 3.0 * eps2ww_ * eps2vw_)
                             + 3.0 * eps3ww_ * eps3vw_
                             + 3.0 * eps4ww_ * eps4vw_;

      // wuuu
      E4CachedValue_[6][0] = E4CachedValue_[0][2];
      // wuuv
      E4CachedValue_[6][1] = E4CachedValue_[1][2];
      // wuuw
      E4CachedValue_[6][2] = E4CachedValue_[0][8];
      // wuvu
      E4CachedValue_[6][3] = E4CachedValue_[0][5];
      // wuvv
      E4CachedValue_[6][4] = E4CachedValue_[1][5];
      // wuvw
      E4CachedValue_[6][5] = E4CachedValue_[1][8];
      // wuwu
      E4CachedValue_[6][6] = E4CachedValue_[0][8];
      // wuwv
      E4CachedValue_[6][7] = E4CachedValue_[1][8];
      // wuww
      E4CachedValue_[6][8] = E4CachedValue_[2][8];

      // wvuu
      E4CachedValue_[7][0] = E4CachedValue_[0][5];
      // wvuv
      E4CachedValue_[7][1] = E4CachedValue_[1][5];
      // wvuw
      E4CachedValue_[7][2] = E4CachedValue_[1][8];
      // wvvu
      E4CachedValue_[7][3] = E4CachedValue_[1][5];
      // wvvv
      E4CachedValue_[7][4] = E4CachedValue_[4][7];
      // wvvw
      E4CachedValue_[7][5] = E4CachedValue_[4][8];
      // wvwu
      E4CachedValue_[7][6] = E4CachedValue_[1][8];
      // wvwv
      E4CachedValue_[7][7] = E4CachedValue_[4][8];
      // wvww
      E4CachedValue_[7][8] = E4CachedValue_[5][8];

      // wwuu
      E4CachedValue_[8][0] = E4CachedValue_[0][8];
      // wwuv
      E4CachedValue_[8][1] = E4CachedValue_[1][8];
      // wwuw
      E4CachedValue_[8][2] = E4CachedValue_[2][8];
      // wwvu
      E4CachedValue_[8][3] = E4CachedValue_[1][8];
      // wwvv
      E4CachedValue_[8][4] = E4CachedValue_[5][7];
      // wwvw
      E4CachedValue_[8][5] = E4CachedValue_[5][8];
      // wwwu
      E4CachedValue_[8][6] = E4CachedValue_[2][8];
      // wwwv
      E4CachedValue_[8][7] = E4CachedValue_[3][8];
      // wwww
      E4CachedValue_[8][8] = Gamma_ * (3.0 * eps1ww_ * eps1ww_
                                       + 3.0 * eps2ww_ * eps2ww_)
                             + 3.0 * eps3ww_ * eps3ww_
                             + 3.0 * eps4ww_ * eps4ww_;


      Cached_[5] = 1;
      CallCount_[5]++;
   }

   return E4CachedValue_;
}

void FourBarTruss::ExtraTestFunctions(Vector& TF) const
{
   for (int i = 0; i < NumExtraTFs_; ++i)
   {
      TF[i] = (ExtraTestFunctions_[i] - Lambda());
   }
}

void FourBarTruss::Print(ostream& out, PrintDetail const& flag,
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
         out << "FourBarTruss:" << "\n" << "\n";
         out << "Theta:" << setw(W) << Theta_ << "\n";
         out << "Psi:" << setw(W) << Psi_ << "\n";
         out << "Gamma:" << setw(W) << Gamma_ << "\n";

         if (Echo_)
         {
            cout << "FourBarTruss:" << "\n" << "\n";
            cout << "Theta:" << setw(W) << Theta_ << "\n";
            cout << "Psi:" << setw(W) << Psi_ << "\n";
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

ostream& operator<<(ostream& out, FourBarTruss& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}
