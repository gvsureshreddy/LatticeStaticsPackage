#ifndef __Shuffle3Rhombo
#define __Shuffle3Rhombo

#include "LatticeMode.h"
#include "GenericLat.h"
#include <math.h>

class Shuffle3Rhombo : public LatticeMode
{
private:
   GenericLat *Lattice_;

public:
   Shuffle3Rhombo(Lattice *M);

   ~Shuffle3Rhombo() {}

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

   virtual char *ModeName() {return "Shuffle3Rhombo";}

};

#endif