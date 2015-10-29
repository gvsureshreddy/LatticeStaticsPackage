#ifndef RSE__MultiLatticeKIM
#define RSE__MultiLatticeKIM

#include <string>
#include <sstream>
#include "PerlInput.h"
#include "Lattice.h"
#include "UnitCellIterator.h"
#include "CBKinematics.h"
#include "LagrangeCB.h"
#include "SymLagrangeWTransCB.h"
#include "PPSumKIM.h"
#include <CMatrix.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

using namespace std;

class MultiLatticeKIM : public Lattice
{
 private:
  int InternalAtoms_;

  static const int cachesize = 3;
  mutable int Cached_[cachesize];
  double InfluenceDist_;

  // DOF[i] = [F00 F01 F02 F10 F11 F12 F20 F21 F22 S00 S01 S02 S11 S12 S13 S21... ...]
  // if using a FwithTransMapping CBkinematics
  // or
  // DOF[i] = [U11 U22 U33 U12 U13 U23 S11 S12 S13 S21... ...]
  // if using a UwithoutTransMapping CBkinematics
  CBKinematics* CBK_;
  int KillRotations_;
  void* pkim_;
  int numberOfParticles_;
  int numberOfSpecies_;
  int* particleSpecies_;
  double cutoff_;
  double* coords_;
  double* forces_;
  mutable double energy_;
  int atom_;
  int numNei_;
  int* neiList_;
  double* rVecList_;

  mutable int StiffnessYes_;
  mutable int BlochwaveProcess_;
  mutable double E0CachedValue_;
  ostringstream descriptor_file_;
  int KillTranslations_;
  Vector KillOneRotation_;
  int Density_;
  double REFLambda_;
  double NormModulus_;
  double PhiRef_;
  enum LDeriv {L0, DL};
  double Lambda_;
  Matrix Loading_;
  Vector* BodyForce_;
  mutable PPSumKIM LatSum_;


  int TFType_;
  Vector TFLoad_;

  // Misc
  int FastPrint_;
  double ConvexityDX_;

  double energy(LDeriv const& dl = L0) const;
  Vector const& stress(LDeriv const& dl = L0) const;
  Matrix const& stiffness(LDeriv const& dl = L0) const;

  void UpdateKIMValues() const;
  void Write_KIM_descriptor_file(const char** SpeciesList,
                                 int numberParticleTypes_);
  static int get_neigh(void* kimmdl, int* mode, int* request, int* atom,
                       int* numnei, int** nei1atom, double** Rij);

  static int process_dEdr(void* kimmdl, double* dEdr, double* r,
                          double** dx, int* i, int* j);
  static int process_d2Edr2(void* kimmdl, double* d2Edr2, double** r,
                            double** dx, int** i, int** j);

 public:
  const static int DIM3;
  // Virtual Functions required by Lattice
  Vector const& DOF() const
  {
    return CBK_->DOF();
  }

  void SetDOF(Vector const& dof)
  {
    for (int i=0; i<cachesize; ++i) Cached_[i] = 0;
    CBK_->SetDOF(dof);
  }

  // Needed to satisfy Lattice.h
  double Entropy() const {return 0.0;}
  double HeatCapacity() const {return 0.0;}
  Vector const& StressDT() const {return VecDummy;}
  Matrix const& StiffnessDT() const {return MatDummy;}
  double Temp() const { return 0.0;}
  void SetTemp(double const& Ntemp) {/*do nothing*/}

  Vector const& StressDL() const
  {
    return stress(DL);
  }

  Matrix const& StiffnessDL() const
  {
    return stiffness(DL);
  }

  virtual Vector const& E1DLoad() const
  {
    return StressDL();
  }

  double Lambda() const
  {
    return Lambda_;
  }

  double ConjugateToLambda() const;
  void SetLambda(double const& lambda)
  {
    for (int i=0; i<cachesize; ++i) Cached_[i] = 0;
    Lambda_ = lambda;
  }

  virtual double E0() const;
  virtual Vector const& E1() const;
  virtual Matrix const& E2() const;
  virtual Matrix const& E3() const;
  virtual Matrix const& E4() const;
  virtual int CriticalPointInfo(
      int* const CPCrossingNum, int const& TFIndex, Vector const& DrDt,
      int const& CPorBif, int const& NumZeroEigenVals, double const& Tolerance,
      int const& Width, PerlInput const& Input, ostream& out);
  virtual void ExtraTestFunctions(Vector& TF) const;

  virtual void LongWavelengthModuli(
      double const& dk, int const& gridsize, char const* const prefix,
      ostream& out) const;
  virtual void SetParameters(double const* const Vals, int const& ResetRef);
  virtual void SetGridSize(int const& Grid) { /*do nothing*/ };

  void RefineEqbm(double const& Tol, int const& MaxItr, ostream* const out);
  void RefineCubicEqbm(
      double const& Tol, int const& MaxItr, ostream* const out);
  virtual char const* const Type() const
  {
    return "MultiLatticeKIM";
  }

  virtual void Print(ostream& out, PrintDetail const& flag,
                     PrintPathSolutionType const& SolType);

  // Functions provided by MultiLatticeKIM
  MultiLatticeKIM(PerlInput const& Input, int const& Echo, int const& Width);
  ~MultiLatticeKIM();
  double InfluenceDist() const
  {
    return InfluenceDist_;
  }

  void SetInfluenceDist(double const& InfluenceDist)
  {
    InfluenceDist_ = InfluenceDist;
  }

  inline double Del(int const& i, int const& j) const
  {
    return i == j;
  }

  Matrix const& CondensedModuli() const;
  //Vector const& ThermalExpansion() const;
  friend ostream& operator<<(ostream& out, MultiLatticeKIM& A);

 private:
  int FindLatticeSpacing(int const& iter, bool cubicEqbm);

  // member variables used to avoid repeated memory allocation/deallocation
  // and thus, improve performance.

  // Dummy place holders
  Matrix MatDummy;
  Vector VecDummy;

  // E0
  mutable double Phi0_static;
  mutable double Tsq_static[3];
  mutable double Rsq_static[3];

  // E1
  mutable Vector ME1_static;
  mutable double T_static[3];
  mutable double R_static[3];

  // stress
  mutable Vector S_static;

  // E2
  mutable Matrix ME2_static;

  // stiffness
  mutable Matrix Phi2_static;

  // CondensedModuli
  mutable Matrix CM_static;

  // Print
  mutable Vector str_static;
  mutable Matrix stiff_static;
  mutable Matrix CondEV_static;
  mutable Matrix CondModuli_static;
  mutable Vector TestFunctVals_static;
  mutable Vector TestFunctVals_Print;
};

#endif
