#ifndef __NaClShuffle
#define __NaClShuffle

#include "LatticeMode.h"
#include "GenericLat.h"
#include <math.h>

class NaClShuffle : public LatticeMode
{
private:
   GenericLat *Lattice_;

public:
   NaClShuffle(Lattice *M);

   ~NaClShuffle() {}

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

   virtual char *ModeName() {return "NaClShuffle";}

};

#endif
