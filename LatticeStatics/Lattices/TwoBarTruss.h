#ifndef __TwoBarTruss
#define __TwoBarTruss

#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class TwoBarTruss : public Lattice
{
private:
   int DOFS_;
   
   // DOF[i] = [u v]
   Vector DOF_;
   enum LDeriv {L0,DL};
   double Lambda_;
   double Theta_;
   double COSTheta_;
   double SINTheta_;

   int Echo_;
   int Width_;

   int Caching_;
   mutable int Cached_[4];
   mutable double E0CachedValue_;
   mutable Matrix E1CachedValue_;
   mutable Matrix E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable int CallCount_[4];
   
public:
   // Functions provided by TwoBarTruss
   TwoBarTruss(PerlInput const& Input,int const& Echo=1,int const& Width=20);
   ~TwoBarTruss();
   
   // Virtual Functions required by Lattice
   Vector const& DOF() const {return DOF_;}
   void SetDOF(Vector const& dof)
   {DOF_ = dof; if (Caching_) for (int i=0;i<4;++i) Cached_[i]=0;}

   double Lambda() const {return Lambda_;}
   void SetLambda(double const& lambda)
   {Lambda_ = lambda; if (Caching_) for (int i=0;i<4;++i) Cached_[i]=0;}
   
   virtual double E0() const;
   virtual Matrix const& E1() const;
   virtual Matrix const& E1DLoad() const;
   virtual Matrix const& E2() const;
   virtual void Print(ostream& out,PrintDetail const& flag);

   friend ostream& operator<<(ostream& out,TwoBarTruss& A);
   
   // ignore these
   virtual void CriticalPointInfo(Vector const& DrDt,int const& NumZeroEigenVals,
                                  double const& Tolerance,int const& Width,ostream& out) {}
   double Entropy() const {return 0.0;}
   double HeatCapacity() const {return 0.0;}
   Matrix const& StressDT() const {return Empty_;}
   Matrix const& StiffnessDT() const {return Empty_;}
   double Temp() const {return 0.0;}
   void SetTemp(double const& Ntemp) {}
   Matrix const& StressDL() const {return Empty_;}
   Matrix const& StiffnessDL() const {return Empty_;}
   virtual Matrix const& E3() const {return Empty_;}
   virtual Matrix const& E4() const {return Empty_;}
   virtual void SetParameters(double const* const Vals,int const& ResetRef = 1) {}
   virtual void SetGridSize(int const& Grid) {}
   
private:
   // place holder
   Matrix Empty_;
   
};

#endif
