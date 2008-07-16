#ifndef RSE__Dobson
#define RSE__Dobson

#include "PairPotentials.h"

#define D3_jFACT -8000.0
#define D3_RCUTL 1.5
#define D3_RCUTU 1.7

using namespace std;

class Dobson: public PairPotentials
{
private:
   double Eps0_,Eps1_,Sigma0_,Sigma1_,rcut_;
   
public:
   
   Dobson() {};
   Dobson(double const& Eps0,double const& Eps1,double const& Sigma0,double const& Sigma1,
          double const& rcut);
   ~Dobson() {};
   friend ostream& operator<<(ostream& out,Dobson const& A);
   double PairPotential(double const& NTemp,double const& r2,YDeriv const& dy=Y0,
                        TDeriv const& dt=T0) const;
   virtual int GetNoParameters() const {return 5;}
   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const {return "Dobson";}
   
   double const& Eps0() const {return Eps0_;}
   double const& Eps1() const {return Eps1_;}
   double const& Sigma0() const {return Sigma0_;}
   double const& Sigma1() const {return Sigma1_;}
   double const& rcut() const {return rcut_;}
   
   void SetEps0(double const& Eps0) {Eps0_=Eps0;}
   void SetEps1(double const& Eps1) {Eps1_=Eps1;}
   void SetSigma0(double const& Sigma0) {Sigma0_=Sigma0;}
   void SetSigma1(double const& Sigma1) {Sigma1_=Sigma1;}
   void Setrcut(double const& rcut) {rcut_ = rcut;}
   
private:
   double Eps(double const& NTemp,TDeriv const& dt=T0) const;
   double Sigma(double const& NTemp,TDeriv const& dt=T0) const;
   double j(double const& NTemp,double const& r2,YDeriv const& dy=Y0,TDeriv const& dt=T0) const;
   double A(double const& NTemp,TDeriv const& dt=T0) const;
   double B(double const& NTemp,TDeriv const& dt=T0) const;
   
};

#endif
