#ifndef __UniDefTemp3D1MMono
#define __UniDefTemp3D1MMono

#include "LatticeMode.h"
#include "UniDefTempLat.h"
#include <math.h>

class UniDefTemp3D1MMono : public LatticeMode
{
private:
   UniDefTempLat *Lattice_;

public:
   UniDefTemp3D1MMono(UniDefTempLat *M);
   UniDefTemp3D1MMono(Lattice *M);

   ~UniDefTemp3D1MMono() {}

   // Functions required by LatticeMode
   virtual double ModeEnergy() {return Lattice_->Energy();}
   virtual Vector ArcLenRHS(double DS,const Vector &Diff,double Aspect);
   virtual Vector ArcLenDef();
   virtual Vector DrDt(const Vector &Diff) {}
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

   virtual char *ModeName() {return "UniDefTemp3D1MMono";}

};

#endif
