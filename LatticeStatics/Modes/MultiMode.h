#ifndef __MultiMode
#define __MultiMode

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
   MultiMode(Lattice *M,PerlInput &Input);
   
   ~MultiMode() {}
   
   // Functions required by LatticeMode
   virtual double ModeEnergy() {return Lattice_->E0();}
   virtual Vector DrDt(const Vector &Diff);
   
   //----------------------------------------------------------------
   virtual Vector ModeForce();
   virtual Matrix ModeStiffness();
   virtual Vector ModeDOF() {return ModeDOF_;}
   virtual void SetModeDOF(const Vector &dof);
   virtual void UpdateModeDOF(const Vector &dr);
   //----------------------------------------------------------------
   virtual char *ModeName() {return "MultiMode";}

private:
   // "static" member variables
   // UpdateLatticeState
   int size_static;
   Vector DOF_static;
   // ModeForce
   Vector force_static;
   Matrix stress_static;
   // ModeStiffness
   Matrix K_static;
   Matrix Stiff_static;
   Matrix stressdt_static;
};

#endif
