#ifndef RSE__PPSumKIM
#define RSE__PPSumKIM

#include <Matrix.h>
#include <Vector.h>
#include <vector>
#include "PerlInput.h"
#include "CBKinematics.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"

#define PPSUMKIMdatalen 11
#define PPSUMKIMatomstart 0
#define PPSUMKIMdXstart 2
#define PPSUMKIMdxstart 5
#define PPSUMKIMr2start 8
#define PPSUMKIMphi1start 9
#define PPSUMKIMphi2start 10

/***********************************************************************
Helper function declarations
***********************************************************************/


/**************************************************************************
Neighbor Object Class Definition
**************************************************************************/
class neighObject
{
  public:
    neighObject()
    {
      iterVal = 0;
    }
    ~neighObject();
    void resetIterator();
    void incrementIterator();
    int getIterator();

    std::vector<int> numNeigh_;
    std::vector<int> nListAtom_;

  private:
    int iterVal;
};

/***************************************************************************
CBK_KIM class definition
***************************************************************************/
class CBK_KIM
{
private:
  int Recalc_;
  double* InfluenceDist_;
  int InternalAtoms_;
  CBKinematics* CBK_;
  CBKinematics* CBK_F_;

  void* pkim_;

  Vector* BodyForce_;
  ostringstream descriptor_file_;
  mutable int numberOfParticles_;
  mutable int numberOfSpecies_;
  double cutoff_;
  double energy_;
  mutable int StiffnessYes_;
  std::vector<int> particleMap_;
  std::vector<double> coordinates_;
  std::vector<double> forces_;
  std::vector<int> particleSpecies_;

  neighObject neigh_;

  mutable Vector ME1_static;
  mutable Vector ME1_F_static;

  mutable Matrix ME2_static;
  mutable Matrix ME2_F_static;

  void Initialize(PerlInput const& Input);
  void Write_KIM_descriptor_file(const char** SpeciesList,
                            int numberOfSpecies_);
  void ComputeAndUpdate();

public:
  const static int DIM3;

  CBK_KIM()
  {
  }

  ~CBK_KIM();

  CBK_KIM(CBKinematics* const CBK, CBKinematics* const CBK_F,
        PerlInput const& Input);

  void operator()(CBKinematics* const CBK, CBKinematics* const CBK_F,
        PerlInput const& Input);

  void UpdateCoordinatesAndKIMValues();

  void ComputeAndUpdate(int StiffnessYes);

  // get methods
  double get_cutoff()
  {
    return cutoff_;
  }
  Vector get_ME1_static()
  {
    return ME1_static;
  }
  Vector get_ME1_F_static()
  {
    return ME1_F_static;
  }
  Matrix get_ME2_static()
  {
    return ME2_static;
  }
  Matrix get_ME2_F_static()
  {
    return ME2_F_static;
  }
  Vector* get_BodyForce()
  {
    return (BodyForce_);
  }
  double Energy()
  {
    return energy_;
  }

  static int get_neigh(void* kimmdl, int* mode, int* request, int* part,
                int* numnei, int** nei1atom, double** pRij);

  static int process_dEdr(void* kimmdl, double* dEdr, double* r,
                  double** dx, int* iSpec, int* jSpec);
  static int process_d2Edr2(void* kimmdl, double* d2Edr2, double** r,
                   double** dx, int** iSpecs, int** jSpecs);

};

#endif
