#ifndef __LatticeMode
#define __LatticeMode

#include "Lattice.h"

class LatticeMode
{
public:
   virtual ~LatticeMode() {}
   
   virtual double ModeEnergy() = 0;
   virtual Vector DrDt(const Vector &Diff) = 0;
   
   virtual Vector ModeForce() = 0;
   virtual Matrix ModeStiffness() = 0;
   virtual Vector ModeDOF() = 0;
   virtual void SetModeDOF(const Vector &dof) = 0;
   virtual void UpdateModeDOF(const Vector &dr) = 0;
   
   virtual char *ModeName() = 0;
   friend ostream &operator<<(ostream &out, LatticeMode *M)
   {out << M->ModeName(); return out;}
};

#endif
