#ifndef RSE__SymLagrangeCB
#define RSE__SymLagrangeCB

#include <Matrix.h>
#include <Vector.h>
#include "PerlInput.h"
#include "CBKinematics.h"

using namespace std;

class SymLagrangeCB: public CBKinematics
{
private:
   virtual void Reset();
   
public:
   SymLagrangeCB(int const& InternalAtoms,Matrix& RefLattice,Vector* const AtomPositions);
   SymLagrangeCB(PerlInput const& Input,PerlInput::HashStruct const* const ParentHash=0);
   virtual ~SymLagrangeCB() {};
   
#include "UwithoutTransMapping.def"
   
   virtual Vector FractionalPosVec(int const& p) const;
   virtual double DX(double const* const X,int const& p,int const& q,int const& i) const;
   virtual double Dx(double const* const X,int const& p,int const& q,int const& i) const;
   
   virtual double DyDF(double const* const Dx,double const* const DX,int const& r,int const& s)
      const;
   virtual double D2yDFF(double const* const DX,int const& r,int const& s,int const& t,
                         int const& u) const;
   virtual double DyDS(double const* const Dx,int const& p,int const& q,int const& i,
                       int const& j) const;
   virtual double D2yDSS(int const& p,int const& q,int const& i,int const& j,int const& k,
                         int const& l) const;
   virtual double D2yDFS(double const* const Dx,double const* const DX,int const& p,
                         int const& q,int const& i,int const& j,int const& k,int const& l)
      const;
   virtual double D3yDFFS(double const* const DX,int const& p,int const& q,int const& i,
                          int const& j,int const& k,int const& l,int const& m,int const& n)
      const;
   virtual double D3yDSSF(int const& p,int const& q,int const& i,int const& j,int const& k,
                          int const& l,int const& m,int const& n) const;
   virtual double D4yDFFSS(int const& p,int const& q,int const& i,int const& j,int const& k,
                           int const& l,int const& m,int const& n,int const& a,int const& b)
      const;
   
   virtual char const* const IDString() const {return "SymLagrangeCB";}
};

#endif
