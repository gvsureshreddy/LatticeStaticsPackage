#ifndef __UniDefTemp3DOrtho
#define __UniDefTemp3DOrtho

#include "LatticeMode.h"
#include "UniDefTempLat.h"
#include <math.h>

class UniDefTemp3DOrtho : public LatticeMode
{
private:
   UniDefTempLat *Lattice_;

public:
   UniDefTemp3DOrtho(UniDefTempLat *M);
   UniDefTemp3DOrtho(Lattice *M);

   ~UniDefTemp3DOrtho() {}

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

   virtual char *ModeName() {return "UniDefTemp3DOrtho";}

};

#endif