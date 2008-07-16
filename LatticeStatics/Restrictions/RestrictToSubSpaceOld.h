#ifndef RSE__RestrictToSubSpaceOld
#define RSE__RestrictToSubSpaceOld

#include "PerlInput.h"
#include "Restriction.h"
#include "Lattice.h"

#include <cmath>

using namespace std;

class RestrictToSubSpaceOld : public Restriction
{
private:
   Lattice *Lattice_;
   
   int DOFS_;
   int DOFindlen_[DOFMAX];
   int DOFindex_[DOFMAX][DOFMAX];
   Vector DOFMult_[DOFMAX];
   Vector BaselineDOF_;
   Vector DOF_;
   
   void UpdateLatticeState();
public:
   RestrictToSubSpaceOld(Lattice* const M,PerlInput const& Input);
   
   ~RestrictToSubSpaceOld() {}
   
   // Functions required by LatticeMode
   virtual double Energy() const {return Lattice_->E0();}
   virtual Vector const& DrDt(Vector const& Diff) const;
   
   //----------------------------------------------------------------
   virtual Vector const& Force() const;
   virtual Matrix const& Stiffness() const;
   virtual Vector const& DOF() const {return DOF_;}
   virtual void SetDOF(Vector const& dof);
   virtual void UpdateDOF(Vector const& dr);
   //----------------------------------------------------------------
   virtual char const* const Name() const {return "RestrictToSubSpaceOld";}

private:
   // "static" member variables
   // UpdateLatticeState
   int size_static;
   Vector DOF_static;
   // DrDt
   mutable Vector ddt_static;
   // ModeForce
   mutable Vector force_static;
   mutable Matrix stress_static;
   // ModeStiffness
   mutable Matrix K_static;
   mutable Matrix Stiff_static;
   mutable Matrix stressdt_static;
};

#endif
