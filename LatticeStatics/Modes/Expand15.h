#ifndef __Expand15
#define __Expand15
#include "LatticeMode.h"
#include "GenericLat.h"

#include <math.h>

class Expand15 : public LatticeMode
{
private:
   GenericLat *Lattice_;

public:
   Expand15(Lattice *M);

   ~Expand15() {}

   // Functions required by LatticeMode
   virtual double ModeEnergy() {return Lattice_->Energy();}
   virtual Vector ArcLenRHS(double DS,const Vector &Diff,double Aspect);
   virtual Vector ArcLenDef();
   virtual Vector DrDt(const Vector &Diff);
   virtual void ArcLenUpdate(const Vector &newval);
   virtual double ArcLenAngle(Vector Old,Vector New,double Aspect);
   virtual Matrix ArcLenStiffness(const Vector &Diff,double Aspect);

   virtual double ScanningDefParameter();
   virtual void ScanningDefParamUpdate(const double newval);
   virtual double ScanningLoadParameter();
   virtual void ScanningLoadParamUpdate(const double newval);
   virtual double ScanningStressParameter();
   
   virtual Vector ScanningRHS();
   virtual Vector ScanningDef();
   virtual void ScanningUpdate(const Vector &newval);
   virtual Matrix ScanningStiffness();

   virtual char *ModeName() {return "Expand15";}

};

#endif
