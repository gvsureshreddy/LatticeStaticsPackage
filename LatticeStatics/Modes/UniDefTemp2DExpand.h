#ifndef __UniDefTemp2DExpand
#define __UniDefTemp2DExpand

#include "LatticeMode.h"
#include "UniDefTempLat.h"
#include <math.h>

class UniDefTemp2DExpand : public LatticeMode
{
private:
   UniDefTempLat *Lattice_;

public:
   UniDefTemp2DExpand(UniDefTempLat *M);
   UniDefTemp2DExpand(Lattice *M);

   ~UniDefTemp2DExpand() {}

   // Functions required by LatticeMode
   virtual Vector ArcLenRHS(double DS,const Vector &Diff,double Aspect);
   virtual Vector ArcLenDef();
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

   virtual char *ModeName() {return "UniDefTemp2DExpand";}

};

#endif
