#ifndef RSE__LJ
#define RSE__LJ

#include "PairPotentials.h"

using namespace std;

class LJ: public PairPotentials
{
protected:
   mutable int EpsChk_[DTmax];
   mutable int SigmaChk_[DTmax];
   mutable int Gchk_[DYmax][DTmax];
   mutable int Hchk_[DYmax][DTmax];
   
   mutable double EpsVal_[DTmax];
   mutable double SigmaVal_[DTmax];
   mutable double Gval_[DYmax][DTmax];
   mutable double Hval_[DYmax][DTmax];
   
   double Eps0_,Eps1_,Sigma0_,Sigma1_;
   
public:
   
   LJ() {};
   LJ(double const& Eps0,double const& Eps1,double const& Sigma0,double const& Sigma1);
   ~LJ() {};
   friend ostream& operator<<(ostream& out,LJ const& A);
   double PairPotential(double const& NTemp,double const& r2,YDeriv const& dy=Y0,
                        TDeriv const& dt=T0) const;
   virtual int GetNoParameters() const {return 4;}
   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const {return "LJ";}
   
   double const& Eps0() const {return Eps0_;}
   double const& Eps1() const {return Eps1_;}
   double const& Sigma0() const {return Sigma0_;}
   double const& Sigma1() const {return Sigma1_;}
   
   void SetEps0(double const& Eps0) {Eps0_=Eps0;}
   void SetEps1(double const& Eps1) {Eps1_=Eps1;}
   void SetSigma0(double const& Sigma0) {Sigma0_=Sigma0;}
   void SetSigma1(double const& Sigma1) {Sigma1_=Sigma1;}
private:
   double Eps(double const& NTemp,TDeriv const& dt=T0) const;
   inline double e(double const& NTemp,TDeriv const& dt=T0) const
   {return (EpsChk_[dt])? EpsVal_[dt] : EpsVal_[dt]=Eps(NTemp,dt);}
   
   double Sigma(double const& NTemp,TDeriv const& dt=T0) const;
   inline double s(double const& NTemp,TDeriv const& dt=T0) const
   {return (SigmaChk_[dt])? SigmaVal_[dt] : SigmaVal_[dt]=Sigma(NTemp,dt);}
   
   double G(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const;
   inline double g(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const
   {return (Gchk_[dy][dt])? Gval_[dy][dt] : Gval_[dy][dt]=G(NTemp,r2,dy,dt);}
   
   double H(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const;
   inline double h(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const
   {return (Hchk_[dy][dt])? Hval_[dy][dt] : Hval_[dy][dt]=H(NTemp,r2,dy,dt);}
};

#endif
