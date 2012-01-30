#ifndef RSE__CBKinematics
#define RSE__CBKinematics

#include <Matrix.h>
#include <Vector.h>
#include "PerlInput.h"

// Max number of atoms in unit cell. might need to be changed...
#define CBK_MAX_ATOMS 1000

using namespace std;

class CBKinematics
{
private:
   virtual void Reset() = 0;

public:
   static int const DIM3;

   Vector DOF_;
   int InternalAtoms_;
   int NumberofSpecies_;
   int AtomSpecies_[CBK_MAX_ATOMS];
   Matrix RefLattice_;
   Vector* InternalPOS_;

   Matrix F_;
   Matrix S_;

   CBKinematics(int const& InternalAtoms, Matrix& RefLattice, Vector* const AtomPositions);
   CBKinematics(PerlInput const& Input, PerlInput::HashStruct const* const ParentHash = 0);
   virtual ~CBKinematics()
   {
      delete[] InternalPOS_;
   }

   virtual void InfluenceRegion(double* const InfluenceRegion);

   virtual void SetReferenceToCurrent();
   virtual void SetReferenceDOFs();
   int const& InternalAtoms() const
   {
      return InternalAtoms_;
   }

   int const& NumberofSpecies() const
   {
      return NumberofSpecies_;
   }

   int const& AtomSpecies(int const& i) const
   {
      return AtomSpecies_[i];
   }

   Matrix const& RefLattice() const
   {
      return RefLattice_;
   }

   double DeltaVolume() const
   {
      return F_.Det();
   }

   double RefVolume() const
   {
      return RefLattice_.Det();
   }

   Vector const& DOF() const
   {
      return DOF_;
   }

   void SetDOF(Vector const& dof)
   {
      DOF_ = dof; Reset();
   }

   Vector const& AtomPositions(int const& i) const
   {
      return InternalPOS_[i];
   }

   virtual inline int DOFS() const
   {
      return Fsize() + Ssize();
   }

   virtual int Fsize() const = 0;
   virtual int Ssize() const = 0;
   virtual int NoTrans() const = 0;
   virtual int INDF(int const& i, int const& j) const = 0;
   virtual int INDS(int const& i, int const& j) const = 0;
   virtual int INDFF(int const& k, int const& l, int const& m, int const& n) const = 0;
   virtual int INDSS(int const& k, int const& l, int const& m, int const& n) const = 0;
   virtual int INDFS(int const& i, int const& j, int const& m, int const& n) const = 0;
   virtual int INDSF(int const& m, int const& n, int const& i, int const& j) const = 0;

   virtual Vector CurrentLatticeVec(int const& p) const;
   virtual Vector FractionalPosVec(int const& p) const
   {
      return Vector(DIM3, 0.0);
   }

   // computes supercell information.
   // Don't forget to delete the memory allocated for SuperCellInternalPOS!!!!
   void SuperCellInfo(int const SuperCell[3][3], Matrix& SuperCellRefLattice,
                      int& SuperCellInternalAtoms, Vector*& SuperCellInternalPOS,
                      int* const SuperCellAtomSpecies) const;
   virtual double DX(double const* const X, int const& p, int const& q, int const& i) const = 0;
   virtual double Dx(double const* const X, int const& p, int const& q, int const& i) const = 0;

   virtual double DyDF(double const* const Dx, double const* const DX, int const& r,
                       int const& s) const = 0;
   virtual double D2yDFF(double const* const DX, int const& r, int const& s, int const& t,
                         int const& u) const = 0;
   virtual double DyDS(double const* const Dx, int const& p, int const& q, int const& i,
                       int const& j) const = 0;
   virtual double D2yDSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                         int const& l) const = 0;
   virtual double D2yDFS(double const* const Dx, double const* const DX, int const& p,
                         int const& q, int const& i, int const& j, int const& k, int const& l)
   const = 0;
   virtual double D3yDFFS(double const* const DX, int const& p, int const& q, int const& i,
                          int const& j, int const& k, int const& l, int const& m, int const& n)
   const = 0;
   virtual double D3yDSSF(int const& p, int const& q, int const& i, int const& j, int const& k,
                          int const& l, int const& m, int const& n) const = 0;
   virtual double D4yDFFSS(int const& p, int const& q, int const& i, int const& j, int const& k,
                           int const& l, int const& m, int const& n, int const& a, int const& b)
   const = 0;

   inline double Del(int const& i, int const& j) const
   {
      return i == j;
   }

   inline double DELTA(int const& s, int const& p, int const& q) const
   {
      return Del(s, q) - Del(s, p);
   }

   virtual char const* const IDString() const = 0;
   friend ostream& operator<<(ostream& out, CBKinematics const& CBK)
   {
      out << CBK.IDString(); return out;
   }
};

#endif
