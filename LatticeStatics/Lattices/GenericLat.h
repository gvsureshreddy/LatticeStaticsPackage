#ifndef __GenericLat
#define __GenericLat

#include "Lattice.h"

#define DOFMAX 20

class GenericLat : public Lattice
{
public:
   virtual ~GenericLat() {}

   virtual Vector DOF() = 0;
   virtual void SetDOF(const Vector &dof) = 0;
   virtual Matrix StressDT() = 0;
   virtual Matrix StiffnessDT() = 0;
   virtual double Temp() = 0;
   virtual void SetTemp(const double &temp) = 0;
   virtual void CriticalPointInfo(const Vector &DrDt,double Tolerance,
				  char *datafile,int Width,ostream &out);
};

#endif
   
