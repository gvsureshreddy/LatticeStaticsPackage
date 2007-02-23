#ifndef __MultiMode
#define __MultiMode
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
   Vector ModeDOF_;

   int ScnDefParam_;   

public:
   MultiMode(Lattice *M,const char *datafile,const char *prefix);

   ~MultiMode() {}

   // Functions required by LatticeMode
   virtual double ModeEnergy() {return Lattice_->Energy();}
   virtual Vector DrDt(const Vector &Diff);

   //----------------------------------------------------------------
   virtual Vector ModeForce();
   virtual Matrix ModeStiffness();
   virtual Vector ModeDOF() {return ModeDOF_;}
   virtual void SetModeDOF(const Vector &dof);
   virtual void UpdateModeDOF(const Vector &dr);

   //----------------------------------------------------------------
   virtual Vector ArcLenForce(double DS,const Vector &Diff,double Aspect);
   virtual Vector ArcLenDef() {return ModeDOF();}
   virtual void ArcLenSet(const Vector &val) {SetModeDOF(val);}
   virtual void ArcLenUpdate(const Vector &newval) {UpdateModeDOF(newval);}
   virtual double ArcLenAngle(Vector Old,Vector New,double Aspect);
   virtual Matrix ArcLenStiffness(const Vector &Diff,double Aspect);

   //----------------------------------------------------------------
   virtual double ScanningDefParameter();
   virtual void ScanningDefParamSet(const double val);
   virtual void ScanningDefParamUpdate(const double newval);
   virtual double ScanningLoadParameter();
   virtual void ScanningLoadParamSet(const double val);
   virtual void ScanningLoadParamUpdate(const double newval);
   virtual double ScanningStressParameter();
   
   virtual Vector ScanningForce();
   virtual Vector ScanningDef();
   virtual void ScanningSet(const Vector &val);
   virtual void ScanningUpdate(const Vector &newval);
   virtual Matrix ScanningStiffness();

   //----------------------------------------------------------------
   virtual char *ModeName() {return "MultiMode";}

};

#endif
