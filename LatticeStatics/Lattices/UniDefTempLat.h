#ifndef __UniDefTempLat
#define __UniDefTempLat

#include "Lattice.h"

class UniDefTempLat : public Lattice
{
public:
   virtual ~UniDefTempLat() {}

   virtual Matrix DefGrad() = 0;
   virtual void SetDefGrad(const Matrix &defgrad) = 0;
   enum TDeriv {T0,DT};
   virtual Matrix StressDT() = 0;
   virtual double Temp() = 0;
   virtual void SetTemp(const double &temp) = 0;
};

#endif
   
