#ifndef __MultiMode
#define __MultiMode
#include "LatticeMode.h"
#include "Lattice.h"

#include <math.h>

class MultiMode : public LatticeMode
{
private:
   Lattice *Lattice_;

   int DOFS;
   int DOFindlen[DOFMAX];
   int DOFindex[DOFMAX][DOFMAX];
   Vector DOFMult[DOFMAX];

   int ScnDefParam;   

public:
   MultiMode(Lattice *M,const char *datafile,const char *prefix);

   ~MultiMode() {}

   // Functions required by LatticeMode
   virtual double ModeEnergy() {return Lattice_->Energy();}
   virtual Vector ArcLenRHS(double DS,const Vector &Diff,double Aspect);
   virtual Vector ArcLenDef();
   virtual Vector DrDt(const Vector &Diff);
   virtual void ArcLenSet(const Vector &val);
   virtual void ArcLenUpdate(const Vector &newval);
   virtual double ArcLenAngle(Vector Old,Vector New,double Aspect);
   virtual Matrix ArcLenStiffness(const Vector &Diff,double Aspect);

   virtual double ScanningDefParameter();
   virtual void ScanningDefParamSet(const double val);
   virtual void ScanningDefParamUpdate(const double newval);
   virtual double ScanningLoadParameter();
   virtual void ScanningLoadParamSet(const double val);
   virtual void ScanningLoadParamUpdate(const double newval);
   virtual double ScanningStressParameter();
   
   virtual Vector ScanningRHS();
   virtual Vector ScanningDef();
   virtual void ScanningSet(const Vector &val);
   virtual void ScanningUpdate(const Vector &newval);
   virtual Matrix ScanningStiffness();

   virtual char *ModeName() {return "MultiMode";}

};

#endif
