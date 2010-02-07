#ifndef RSE__LJLinearCutoff
#define RSE__LJLinearCutoff

#include "LJ.h"

using namespace std;

class LJLinearCutoff: public LJ
{
private:
   double Cutoff_;
   
public:
   
   LJLinearCutoff() {};
   LJLinearCutoff(double const& Eps0,double const& Eps1,double const& Sigma0,double const& Sigma1,
            double const& Cutoff);
   ~LJLinearCutoff() {};
   friend ostream& operator<<(ostream& out,LJLinearCutoff const& A);
   double PairPotential(double const& NTemp,double const& r2,YDeriv const& dy=Y0,
                        TDeriv const& dt=T0) const;
   virtual int GetNoParameters() const {return (1+LJ::GetNoParameters());}
   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const {return "LJLinearCutoff";}
   
   double const& Cutoff() const {return Cutoff_;}
   
   void SetCutoff(double const& Cutoff) {Cutoff_=Cutoff;}
private:
   
   double CutoffFunction(double const& NTemp,double const& r2,YDeriv const& dy=Y0,
                         TDeriv const& dt=T0) const;
   
};

#endif
