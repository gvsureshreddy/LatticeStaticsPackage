#ifndef RSE__MultiMode
#define RSE__MultiMode

#include "PerlInput.h"
#include "LatticeMode.h"
#include "Lattice.h"

#include <cmath>

using namespace std;

class MultiMode : public LatticeMode
{
private:
   Lattice *Lattice_;
   
   int DOFS_;
   int DOFindlen_[DOFMAX];
   int DOFindex_[DOFMAX][DOFMAX];
   Vector DOFMult_[DOFMAX];
   Vector BaselineDOF_;
   Vector ModeDOF_;
   
   void UpdateLatticeState();
public:
   MultiMode(Lattice* const M,PerlInput const& Input);
   
   ~MultiMode() {}
   
   // Functions required by LatticeMode
   virtual double ModeEnergy() const {return Lattice_->E0();}
   virtual Vector const& DrDt(Vector const& Diff) const;
   
   //----------------------------------------------------------------
   virtual Vector const& ModeForce() const;
   virtual Matrix const& ModeStiffness() const;
   virtual Vector const& ModeDOF() const {return ModeDOF_;}
   virtual void SetModeDOF(Vector const& dof);
   virtual void UpdateModeDOF(Vector const& dr);
   //----------------------------------------------------------------
   virtual char const* const ModeName() const {return "MultiMode";}

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
