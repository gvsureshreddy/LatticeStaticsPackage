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
   double E0() {return Energy();}
   virtual Matrix Stress() = 0;
   Matrix E1() {return Stress();}
   virtual Matrix Stiffness() = 0;
   Matrix E2() {return Stiffness();}
   virtual Matrix Moduli() = 0;
   virtual Matrix E3() = 0;
   virtual Matrix E4() = 0;
   virtual int StiffnessNulity(double *Min=NULL);
   virtual void CriticalPointInfo(const Vector &DrDt,double Tolerance,
				  char *datafile,int Width,ostream &out) = 0;

   enum PrintDetail {PrintLong,PrintShort};
   virtual void Print(ostream &out,PrintDetail flag) = 0;
   friend ostream &operator<<(ostream &out,Lattice *L)
   {L->Print(out,PrintShort); return out;}
};

#endif
