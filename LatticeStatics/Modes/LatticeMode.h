#ifndef RSE__LatticeMode
#define RSE__LatticeMode

#include "Lattice.h"
#define LINELENGTH 600

class LatticeMode
{
public:
   virtual ~LatticeMode() {}
   
   virtual double ModeEnergy() const = 0;
   virtual Vector const& DrDt(Vector const& Diff) const = 0;
   
   virtual Vector const& ModeForce() const = 0;
   virtual Matrix const& ModeStiffness() const = 0;
   virtual Vector const& ModeDOF() const = 0;
   virtual void SetModeDOF(Vector const& dof) = 0;
   virtual void UpdateModeDOF(Vector const& dr) = 0;
   
   virtual char const* const ModeName() const = 0;
   friend ostream& operator<<(ostream& out,LatticeMode const& M)
   {out << M.ModeName(); return out;}
};

#endif
