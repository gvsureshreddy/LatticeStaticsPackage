#ifndef __Lattice
#define __Lattice

#include <Matrix.h>
#include <Vector.h>
#include <iostream.h>
#include <iomanip.h>


class Lattice
{
public:
   virtual ~Lattice() {}

   virtual double Energy() = 0;
   virtual Matrix Stress() = 0;
   virtual Matrix Stiffness() = 0;
   virtual Matrix Moduli() = 0;
   virtual int StiffnessNulity(double *Min=NULL) = 0;

   enum PrintDetail {PrintLong,PrintShort};
   virtual void Print(ostream &out,PrintDetail flag) = 0;
   friend ostream &operator<<(ostream &out,Lattice *L)
   {L->Print(out,PrintShort); return out;}
};

#endif
