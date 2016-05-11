#include "RandomAlloy.h"
#include "KnownPairPotentials.h"
#include <cstdlib>

RandomAlloy::RandomAlloy(double const& CompA, double const& CompB, PerlInput const& Input,
                         char const* const HashName) :
   CompOfSubLatA_(CompA), CompOfSubLatB_(CompB)
{
   if ((CompOfSubLatA_ < 0.0) || (CompOfSubLatA_ > 1.0))
   {
      cerr << "RandomAlloy: CompOfSubLatA_ out of range!" << endl;
      exit(-1);
   }
   if ((CompOfSubLatB_ < 0.0) || (CompOfSubLatB_ > 1.0))
   {
      cerr << "RandomAlloy: CompOfSubLatB_ out of range!" << endl;
      exit(-1);
   }

   BinaryAlloyPots[0] = InitializePairPotential(HashName, Input, 0, 0);
   BinaryAlloyPots[1] = InitializePairPotential(HashName, Input, 1, 1);
   BinaryAlloyPots[2] = InitializePairPotential(HashName, Input, 0, 1);
}

double RandomAlloy::PairPotential(double const& NTemp, double const& r2, YDeriv const& dy,
                                  TDeriv const& dt) const
{
   return ((1.0 - CompOfSubLatA_) * (1.0 - CompOfSubLatB_)
           * BinaryAlloyPots[0]->PairPotential(NTemp, r2, dy, dt))
          + (CompOfSubLatA_ * CompOfSubLatB_ * BinaryAlloyPots[1]->PairPotential(NTemp, r2, dy, dt))
          + (((1.0 - CompOfSubLatA_) * CompOfSubLatB_ + (1.0 - CompOfSubLatB_) * CompOfSubLatA_)
             * BinaryAlloyPots[2]->PairPotential(NTemp, r2, dy, dt));
}

void RandomAlloy::Print(ostream& out) const
{
   int W = out.width();

   out.width(0);

   out << "CompOfSubLatA=" << setw(W) << CompOfSubLatA_
       << "; CompOfSubLatB=" << setw(W) << CompOfSubLatB_ << "\n";
   for (int i = 0; i < BINARYPAIRBONDS; ++i)
   {
      out << "\tBond" << i << " -- " << BinaryAlloyPots[i]->Type() << " -- " << setw(W);
      BinaryAlloyPots[i]->Print(out);
      out << "\n";
   }
}

ostream& operator<<(ostream& out, RandomAlloy const& A)
{
   A.Print(out);
   return out;
}
