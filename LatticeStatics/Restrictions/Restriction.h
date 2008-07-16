#ifndef RSE__Restriction
#define RSE__Restriction

#include "Lattice.h"
#define LINELENGTH 600

class Restriction
{
public:
   virtual ~Restriction() {}
   
   virtual double Energy() const = 0;
   virtual Vector const& DrDt(Vector const& Diff) const = 0;
   
   virtual Vector const& Force() const = 0;
   virtual Matrix const& Stiffness() const = 0;
   virtual Vector const& DOF() const = 0;
   virtual void SetDOF(Vector const& dof) = 0;
   virtual void UpdateDOF(Vector const& dr) = 0;
   
   virtual char const* const Name() const = 0;
   friend ostream& operator<<(ostream& out,Restriction const& R)
   {out << R.Name(); return out;}
};

#endif
