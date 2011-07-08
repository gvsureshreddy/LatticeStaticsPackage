#ifndef RSE__FourBarTruss
#define RSE__FourBarTruss

// Author: Karthikreddy Ginnavaram

#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class FourBarTruss : public Lattice
{
private:
   int DOFS_;

   // DOF[i] = [u v w]
   Vector DOF_;
   enum LDeriv {L0, DL};
   double Lambda_;  // applied load
   double Gamma_;   // imperfection in modulus
   double Theta_;
   double Psi_;
   double COSTheta_;
   double SINTheta_;
   double COSPsi_;
   double SINPsi_;

   int Width_;

   int Caching_;
   static const int cachesize = 6;
   mutable int Cached_[cachesize];
   mutable double E0CachedValue_;
   mutable Vector E1CachedValue_;
   mutable Vector E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable Matrix E3CachedValue_;
   mutable Matrix E4CachedValue_;
   mutable Vector ExtraTestFunctions_;
   mutable int CallCount_[cachesize];

public:
   // Functions provided by FourBarTruss
   FourBarTruss(PerlInput const& Input, int const& Echo = 1, int const& Width = 20);
   ~FourBarTruss();

   // Virtual Functions required by Lattice
   Vector const& DOF() const
   {
      return DOF_;
   }

   void SetDOF(Vector const& dof)
   {
      DOF_ = dof; if (Caching_)
      {
         for (int i = 0; i < cachesize; ++i)
         {
            Cached_[i] = 0;
         }
      }
   }

   double Lambda() const
   {
      return Lambda_;
   }

   void SetLambda(double const& lambda)
   {
      Lambda_ = lambda; if (Caching_)
      {
         for (int i = 0; i < cachesize; ++i)
         {
            Cached_[i] = 0;
         }
      }
   }

   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Vector const& E1DLoad() const;
   virtual Matrix const& E2() const;
   virtual Matrix const& E3() const;
   virtual Matrix const& E4() const;
   virtual void ExtraTestFunctions(Vector& TF) const;
   virtual char const* const Type() const
   {
      return "FourBarTruss";
   }
   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   friend ostream& operator<<(ostream& out, FourBarTruss& A);

   // ignore these
   double Entropy() const
   {
      return 0.0;
   }

   double HeatCapacity() const
   {
      return 0.0;
   }

   Vector const& StressDT() const
   {
      return EmptyV_;
   }

   Matrix const& StiffnessDT() const
   {
      return EmptyM_;
   }

   double Temp() const
   {
      return 0.0;
   }

   void SetTemp(double const& Ntemp)
   {
   }

   Vector const& StressDL() const
   {
      return E1DLoad();
   }

   Matrix const& StiffnessDL() const
   {
      return EmptyM_;
   }

   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
   {
   }

   virtual void SetGridSize(int const& Grid)
   {
   }

private:
   // temp storage space
   mutable double eps1_;
   mutable double eps2_;
   mutable double eps3_;
   mutable double eps4_;
   mutable double eps1u_;
   mutable double eps2u_;
   mutable double eps3u_;
   mutable double eps4u_;
   mutable double eps1v_;
   mutable double eps2v_;
   mutable double eps3v_;
   mutable double eps4v_;
   mutable double eps1w_;
   mutable double eps2w_;
   mutable double eps3w_;
   mutable double eps4w_;
   mutable double eps1uu_;
   mutable double eps1vv_;
   mutable double eps1ww_;
   mutable double eps1uv_;
   mutable double eps1uw_;
   mutable double eps1vw_;
   mutable double eps1vu_;
   mutable double eps1wu_;
   mutable double eps1wv_;
   mutable double eps2uu_;
   mutable double eps2vv_;
   mutable double eps2ww_;
   mutable double eps2uv_;
   mutable double eps2uw_;
   mutable double eps2vw_;
   mutable double eps2vu_;
   mutable double eps2wu_;
   mutable double eps2wv_;
   mutable double eps3uu_;
   mutable double eps3vv_;
   mutable double eps3ww_;
   mutable double eps3uv_;
   mutable double eps3uw_;
   mutable double eps3vw_;
   mutable double eps3vu_;
   mutable double eps3wu_;
   mutable double eps3wv_;
   mutable double eps4uu_;
   mutable double eps4vv_;
   mutable double eps4ww_;
   mutable double eps4uv_;
   mutable double eps4uw_;
   mutable double eps4vw_;
   mutable double eps4vu_;
   mutable double eps4wu_;
   mutable double eps4wv_;

   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;
};

#endif
