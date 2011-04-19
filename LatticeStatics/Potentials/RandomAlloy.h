#ifndef RSE__RandomAlloy
#define RSE__RandomAlloy

#include "PerlInput.h"
#include "PairPotentials.h"

using namespace std;

#define BINARYPAIRBONDS 3

class RandomAlloy : public PairPotentials
{
private:
   double CompOfSubLatA_; // Fraction of type 1 atoms (0.0 - 100% type 0; 1.0 - 100% type 1)
   double CompOfSubLatB_; // Fraction of type 1 atoms (0.0 - 100% type 0; 1.0 - 100% type 1)
   PairPotentials* BinaryAlloyPots[BINARYPAIRBONDS];

public:
   RandomAlloy()
   {
   }
   RandomAlloy(double const& CompA, double const& CompB, PerlInput const& Input,
               char const* const HashName);
   ~RandomAlloy()
   {
      for (int i = 0; i < BINARYPAIRBONDS; ++i)
      {
         delete BinaryAlloyPots[i];
      }
   }
   friend ostream& operator<<(ostream& out, RandomAlloy const& A);
   double PairPotential(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                        TDeriv const& dt = T0) const;
   // Don't bother with Get/Set-Parameters
   virtual int GetNoParameters() const
   {
      return 0;
   }
   virtual void SetParameters(double const* const Vals)
   {
   }
   // ....
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const
   {
      return "RandomAlloy";
   }

   void SetCompOfSubLatA(double const& CompA)
   {
      CompOfSubLatA_ = CompA;
   }
   void SetCompOfSubLatB(double const& CompB)
   {
      CompOfSubLatB_ = CompB;
   }
private:
};

#endif

