#ifndef __LatticeMode
#define __LatticeMode

#include "Lattice.h"

#define DIM 2

class LatticeMode
{
public:
   virtual ~LatticeMode() {}

   virtual Vector ArcLenRHS(double DS,const Vector &Diff,double Aspect) = 0;
   virtual Vector ArcLenDef() = 0;
   virtual void ArcLenUpdate(const Vector &newval) = 0;
   virtual double ArcLenAngle(Vector Old,Vector New,double Aspect) = 0;
   virtual Matrix ArcLenStiffness(const Vector &Diff,double Aspect) = 0;

   virtual double ScanningDefParameter() = 0;
   virtual void ScanningDefParamUpdate(const double newval) = 0;
   virtual double ScanningLoadParameter() = 0;
   virtual void ScanningLoadParamUpdate(const double newval) = 0;
   virtual double ScanningStressParameter() = 0;
   
   virtual Vector ScanningRHS() = 0;
   virtual Vector ScanningDef() = 0;
   virtual void ScanningUpdate(const Vector &newval) = 0;
   virtual Matrix ScanningStiffness() = 0;

   virtual char *ModeName() = 0;
   friend ostream &operator<<(ostream &out, LatticeMode *M)
   {out << M->ModeName(); return out;}
};

#endif
