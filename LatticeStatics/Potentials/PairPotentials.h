#ifndef RSE__PairPotentials
#define RSE__PairPotentials

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;


class PairPotentials
{
public:
   enum YDeriv {Y0, DY, D2Y, D3Y, D4Y, D5Y, D6Y, DYmax};
   enum TDeriv {T0, DT, D2T, DTmax};

   virtual ~PairPotentials()
   {
   }

   virtual double PairPotential(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                                TDeriv const& dt = T0) const = 0;
   virtual int GetNoParameters() const = 0;
   virtual void SetParameters(double const* const Vals) = 0;
   virtual char const* const Type() const = 0;
   virtual void Print(ostream& out) const = 0;
   friend ostream& operator<<(ostream& out, PairPotentials const& PP)
   {
      PP.Print(out); return out;
   }
};

#endif

