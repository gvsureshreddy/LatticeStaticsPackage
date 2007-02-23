#ifndef __LatticeMode
#define __LatticeMode

#include "Lattice.h"

class LatticeMode
{
public:
   virtual ~LatticeMode() {}

   virtual double ModeEnergy() = 0;
   virtual Vector DrDt(const Vector &Diff) = 0;

   virtual Vector ModeForce() = 0;
   virtual Matrix ModeStiffness() = 0;
   virtual Vector ModeDOF() = 0;
   virtual void SetModeDOF(const Vector &dof) = 0;
   virtual void UpdateModeDOF(const Vector &dr) = 0;

   virtual Vector ArcLenForce(double DS,const Vector &Diff,double Aspect) = 0;
   virtual Vector ArcLenDef() = 0;
   virtual void ArcLenSet(const Vector &val) = 0;
   virtual void ArcLenUpdate(const Vector &newval) = 0;
   virtual double ArcLenAngle(Vector Old,Vector New,double Aspect) = 0;
   virtual Matrix ArcLenStiffness(const Vector &Diff,double Aspect) = 0;

   virtual double ScanningDefParameter() = 0;
   virtual void ScanningDefParamSet(const double val) = 0;
   virtual void ScanningDefParamUpdate(const double newval) = 0;
   virtual double ScanningLoadParameter() = 0;
   virtual void ScanningLoadParamSet(const double val) = 0;
   virtual void ScanningLoadParamUpdate(const double newval) = 0;
   virtual double ScanningStressParameter() = 0;
   
   virtual Vector ScanningForce() = 0;
   virtual Vector ScanningDef() = 0;
   virtual void ScanningSet(const Vector &val) = 0;
   virtual void ScanningUpdate(const Vector &newval) = 0;
   virtual Matrix ScanningStiffness() = 0;

   virtual char *ModeName() = 0;
   friend ostream &operator<<(ostream &out, LatticeMode *M)
   {out << M->ModeName(); return out;}
};

#endif
