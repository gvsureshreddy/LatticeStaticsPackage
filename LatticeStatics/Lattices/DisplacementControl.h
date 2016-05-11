#ifndef RSE__DisplacementControl
#define RSE__DisplacementControl

#include "PerlInput.h"
#include "Lattice.h"
#include <cstdlib>

using namespace std;

class DisplacementControl : public Lattice
{
 private:
  Lattice* Lat_;
  int LatDOFS_;
  int DispDOFS_;
  Vector Proportion_;
  Vector Reference_;
  int* DispIndex_;
  int* FreeIndex_;

  int DOFS_;
  Vector DOF_;
  enum LDeriv {L0, DL};
  double Lambda_;

  int Width_;

  void SetLatDOF();

  static const int cachesize = 4;
  mutable int Cached_[cachesize];
  mutable double E0CachedValue_;
  mutable Vector E1CachedValue_;
  mutable Vector E1DLoadCachedValue_;
  mutable Matrix E2CachedValue_;
  mutable int CallCount_[cachesize];

 public:
  // Functions provided by DisplacementControl
  DisplacementControl(PerlInput const& Input, int const& Echo = 1,
                      int const& Width = 20);
  ~DisplacementControl();

  // Virtual Functions required by Lattice
  Vector const& DOF() const
  {
    return DOF_;
  }

  void SetDOF(Vector const& dof)
  {
    DOF_ = dof;
    SetLatDOF();
    for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
    }
  }

  double Lambda() const
  {
    return Lambda_;
  }

  void SetLambda(double const& lambda)
  {
    Lambda_ = lambda;
    SetLatDOF();
    for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
    }
  }

  virtual double E0() const;
  virtual Vector const& E1() const;
  virtual Vector const& E1DLoad() const;
  virtual Matrix const& E2() const;
  virtual char const* const Type() const
  {
    return "DisplacementControl";
  }

  virtual void Print(ostream& out, PrintDetail const& flag,
                     PrintPathSolutionType const& SolType = RegularPt);

  friend ostream& operator<<(ostream& out, DisplacementControl& A);

  // ignore these
  double Entropy() const
  {
    return 0.0;
  }

  double HeatCapacity() const
  {
    return 0.0;
  }

  Vector const& StressDT() const
  {
    return EmptyV_;
  }

  Matrix const& StiffnessDT() const
  {
    return EmptyM_;
  }

  double Temp() const
  {
    return 0.0;
  }

  void SetTemp(double const& Ntemp)
  {
  }

  Vector const& StressDL() const
  {
    return E1DLoad();
  }

  Matrix const& StiffnessDL() const
  {
    return EmptyM_;
  }

  virtual Matrix const& E3() const
  {
    cerr << "DisplacementControl::E3() Not Programmed\n"; exit(-1); return EmptyM_;
  }

  virtual Matrix const& E4() const
  {
    cerr << "DisplacementControl::E4() Not Programmed\n"; exit(-1); return EmptyM_;
  }

  virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
  {
  }

  virtual void SetGridSize(int const& Grid)
  {
  }

 private:
  // place holder
  Vector EmptyV_;
  Matrix EmptyM_;
};

#endif
