#ifndef __UniDefTemp3DTetrag
#define __UniDefTemp3DTetrag

#include "LatticeMode.h"
#include "UniDefTempLat.h"
#include <math.h>

class UniDefTemp3DTetrag : public LatticeMode
{
private:
   UniDefTempLat *Lattice_;

public:
   UniDefTemp3DTetrag(UniDefTempLat *M);
   UniDefTemp3DTetrag(Lattice *M);

   ~UniDefTemp3DTetrag() {}

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

   virtual char *ModeName() {return "UniDefTemp3DTetrag";}

};

#endif
