#ifndef __OrthoShuffle
#define __OrthoShuffle
#include "LatticeMode.h"
#include "NiTiShuffleTPPLat.h"

#include <math.h>

class OrthoShuffle : public LatticeMode
{
private:
   NiTiShuffleTPPLat *Lattice_;

public:
   OrthoShuffle(Lattice *M);

   ~OrthoShuffle() {}

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

   virtual char *ModeName() {return "OrthoShuffle";}

};

#endif