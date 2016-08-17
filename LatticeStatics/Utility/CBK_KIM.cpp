#include "CBK_KIM.h"

#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>

using namespace std;

int const CBK_KIM::DIM3 = 3;

/**********************************************************************
Neighbor Object Definition
**********************************************************************/

neighObject::~neighObject()
{
  numNeigh_.clear();
  nListAtom_.clear();
}
void neighObject::resetIterator()
{
  iterVal = 0;
}
void neighObject::incrementIterator()
{
  iterVal++;
}
int neighObject::getIterator()
{
  return iterVal;
}

/**********************************************************************
CBK_KIM Object Definition
**********************************************************************/
CBK_KIM::~CBK_KIM()
{
  delete[] BodyForce_;
  int status;

  status = KIM_API_model_destroy(pkim_);
  if (status < KIM_STATUS_OK)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_model_destroy", status);
  }
  KIM_API_free(&pkim_, &status);
  if(status < KIM_STATUS_OK)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_free", status);
  }
}

void CBK_KIM::operator()(CBKinematics* const CBK, CBKinematics* const CBK_F,
      PerlInput const& Input)
{
  CBK_ = CBK;
  CBK_F_ = CBK_F;

  Initialize(Input);
}

void CBK_KIM::Initialize(PerlInput const& Input)
{
  // Get Lattice definition
  PerlInput::HashStruct Hash = Input.getHash("Lattice", "MultiLatticeKIM");
  PerlInput::HashStruct CBKHash = Input.getHash(Hash, "CBKinematics");

  InternalAtoms_ = CBK_->InternalAtoms();

  // set up body force
  BodyForce_ = new Vector[InternalAtoms_];

  for (int i = 0; i < InternalAtoms_; ++i)
  {
    BodyForce_[i].Resize(DIM3, 0.0);
  }

  int status;
  char* modelname;
  if (Input.ParameterOK(Hash, "KIMModel"))
  {
   // Reads in the name of the model from Input file
   modelname = (char*) Input.getString(Hash, "KIMModel");
  }
  else
  {
   cerr << "No KIM Model in input file" << "\n";
   exit(-1);
  }

  // @@ NumberOfSpecies and SpeciesList need to be moved to CBKinematics object
  // Gets total number of species. From CBKHash
  numberOfSpecies_ = Input.getInt(CBKHash, "NumberOfSpecies");

  const char** SpeciesList = new const char*[numberOfSpecies_];
  const char** AtomSpeciesList = new const char*[InternalAtoms_];

  for (int i = 0; i < numberOfSpecies_; i++)
  {
   // Gets list of Species. From CBKHash
   SpeciesList[i] = Input.getString(CBKHash, "SpeciesList", i);
  }

  // Create empty KIM object conforming to fields in the KIM descriptor files
  // of the Test and Model
  Write_KIM_descriptor_file(SpeciesList, numberOfSpecies_);
  // Calls function to create a compatible descriptor file. Will need to
  // augment so that it can read in info from a model and automatically
  // decides the appropriate tests.
  if (Input.ParameterOK(Hash, "PrintKIM_DescriptorFile"))
  {
   char const* const
      printKIM = Input.getString(Hash, "PrintKIM_DescriptorFile");
   if ((!strcmp(printKIM, "Yes")) || (!strcmp(printKIM, "yes")))
   {
    cout << descriptor_file_.str();
   }
  }
  char* Test_Descriptor_file = new char[descriptor_file_.str().length() + 1];
  strcpy(Test_Descriptor_file, descriptor_file_.str().c_str());

  status = KIM_API_string_init(&pkim_, Test_Descriptor_file, modelname);
  if (KIM_STATUS_OK > status)
  {
   KIM_API_report_error(__LINE__, (char*) __FILE__,
                 (char*) "Test-Model coupling failure "
                 "(see kim.log file for details).", status);
   exit(1);
  }
  delete [] Test_Descriptor_file;

  KIM_API_set_sim_buffer(pkim_, (void*) this, &status);
  if (KIM_STATUS_OK > status)
  {
   KIM_API_report_error(__LINE__, (char*) __FILE__,
                 (char*) "Test-Model coupling failure "
                 "(see kim.log file for details).", status);
   exit(1);
  }

  // Get the cutoff
  KIM_API_set_data(pkim_, "cutoff", 1, &cutoff_);
  status = KIM_API_model_init(pkim_);
  if (KIM_STATUS_OK > status)
  {
   KIM_API_report_error(__LINE__, (char*) __FILE__,
                 (char*) "KIM_API_model_init", status);
  }

  // figures out which atoms corresponds to each element in the species list
  // for each of the internal atoms
  for (int i = 0; i < InternalAtoms_; i++)
  {
    AtomSpeciesList[i] = Input.getString(CBKHash, "AtomSpeciesKIM", i);
    for (int j = 0; j < numberOfSpecies_; j++)
    {
      if (!strcmp(AtomSpeciesList[i], SpeciesList[j]))
      {
        particleSpecies_.push_back(KIM_API_get_species_code(
            pkim_, (char*) SpeciesList[j], &status));
        if (KIM_STATUS_OK > status)
        {
          KIM_API_report_error(__LINE__, (char*) __FILE__,
                               (char*) "KIM_API_get_partcl_type_code",
                               status);
          exit(1);
        }
        j = numberOfSpecies_ + 1;
      }
    }
  }

  // Done with SpeciesList and AtomSpeciesList
  delete [] SpeciesList;
  delete [] AtomSpeciesList;

  // set kim data and calcualte coords
  UpdateCoordinatesAndKIMValues();

  KIM_API_setm_data(
      pkim_, &status, 3 * 4,
      "numberContributingParticles", 1, &InternalAtoms_, 1,
      "numberOfSpecies", 1, &numberOfSpecies_, 1,
      "energy", 1, &energy_, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_data", status);
  }
  KIM_API_setm_method(
      pkim_, & status, 3 * 4,
      "get_neigh", 1, (func_ptr) &CBK_KIM::get_neigh, 1,
      "process_dEdr", 1, (func_ptr) &CBK_KIM::process_dEdr, 1,
      "process_d2Edr2", 1, (func_ptr) &CBK_KIM::process_d2Edr2, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_method", status);
  }

  // reset any KIM Model Published parameters as requested
  char const KIMparams[] = "KIMModelPublishedParameters";
  if (Input.ParameterOK(Hash, KIMparams))
  {
    int const params = Input.getArrayLength(Hash, KIMparams);

    for (int i = 0; i < params; ++i)
    {
      char const* const paramName = Input.getString(Hash, KIMparams, i, 0);
      char const* const paramType = Input.getString(Hash, KIMparams, i, 1);
      cout << "paramName is: " << paramName << "\n"
           << "paramType is: " << paramType << "\n";
      if (!strcmp("integer", paramType))
      {
        int const val = Input.getInt(Hash, KIMparams, i, 2);

        int * const paramVal = (int*)
            KIM_API_get_data(pkim_, paramName, &status);
        if (KIM_STATUS_OK > status)
        {
          KIM_API_report_error(__LINE__, (char*) __FILE__,
                               (char*) "KIM_API_get_data",
                               status);
        }
        else
        {
          *paramVal = val;

          status = KIM_API_model_reinit(pkim_);
          if (KIM_STATUS_OK > status)
          {
            KIM_API_report_error(__LINE__, (char*) __FILE__,
                                 (char*) "KIM_API_model_reinit",
                                 status);
            exit(-1);
          }
        }
      }
      else if (!strcmp("double", paramType))
      {
        double const val = Input.getDouble(Hash, KIMparams, i, 2);
        double * const paramVal = (double*)
            KIM_API_get_data(pkim_, paramName, &status);
        if (KIM_STATUS_OK > status)
        {
          KIM_API_report_error(__LINE__, (char*) __FILE__,
                               (char*) "KIM_API_get_data",
                               status);
        }
        else
        {
          *paramVal = val;

          status = KIM_API_model_reinit(pkim_);
          if (KIM_STATUS_OK > status)
          {
            KIM_API_report_error(__LINE__, (char*) __FILE__,
                                 (char*) "KIM_API_model_reinit",
                                 status);
            exit(-1);
          }
        }
      }
      else
      {
        cerr << "Error (MultiLatticeKIM()): "
             << "Unknown KIM published parameter type"
             << "\n";
        exit(-2);
      }
    }
  }

  // set the neighbor object
  KIM_API_set_data(pkim_, "neighObject", 1, &neigh_);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_data", status);
  }

  // intialize various data storage space
  ME1_static.Resize(CBK_->DOFS(), 0.0);
  ME1_F_static.Resize(CBK_F_->DOFS(), 0.0);
  ME2_static.Resize(CBK_->DOFS(), CBK_->DOFS(), 0.0);
  ME2_F_static.Resize(CBK_F_->DOFS(), CBK_F_->DOFS(), 0.0);
}

void CBK_KIM::Write_KIM_descriptor_file(const char** SpeciesList,
                                int numberOfSpecies_)
{
  descriptor_file_ << "#\n"
    "# BEGINNING OF KIM DESCRIPTOR FILE\n"
    "#\n"
    "# This file is automatically generated from MultiLatticeKIM.\n"
    "#################################################" << endl;
  descriptor_file_ << "KIM_API_Version  := 1.6.0" << endl;
  descriptor_file_ << "Unit_length    := A" << endl;
  descriptor_file_ << "Unit_energy    := eV" << endl;
  descriptor_file_ << "Unit_charge    := e" << endl;
  descriptor_file_ << "Unit_temperature := K" << endl;
  descriptor_file_ << "Unit_time      := ps" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "SUPPORTED_ATOM/PARTICLES_TYPES:" << endl;
  for (int i = 0; i < numberOfSpecies_; i++)
  {
   descriptor_file_ << SpeciesList[i] << " spec  1" << endl;
  }

  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "CONVENTIONS:" << endl;
  descriptor_file_ << "# Name               Type" << endl;
  descriptor_file_ << "ZeroBasedLists          flag" << endl;
  descriptor_file_ << "Neigh_BothAccess        flag" << endl;
  //    descriptor_file_ << "CLUSTER              flag" << endl;
  descriptor_file_ << "NEIGH_PURE_F           flag" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "MODEL_INPUT:" << endl;
  descriptor_file_ << "# Name               Type      Unit"
    "           Shape" << endl;
  descriptor_file_ << "numberOfParticles        integer    none"
    "           []" << endl;
  descriptor_file_ << "numberOfSpecies         integer    none"
    "           []" << endl;
  descriptor_file_ << "particleSpecies         integer    none"
    "           [numberOfParticles]" << endl;
  descriptor_file_ << "coordinates            double     length"
    "          [numberOfParticles,3]" << endl;
  descriptor_file_ << "get_neigh             method     none"
    "           []" << endl;
  descriptor_file_ << "neighObject            pointer    none"
    "           []" << endl;
  descriptor_file_ << "process_dEdr           method     none"
    "           []" << endl;
  descriptor_file_ << "process_d2Edr2          method     none"
    "           []" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "MODEL_OUTPUT:" << endl;
  descriptor_file_ << "# Name               Type      Unit"
    "           Shape" << endl;
  descriptor_file_ << "destroy              method     none"
    "           []" << endl;
  descriptor_file_ << "compute              method     none"
    "           []" << endl;
  descriptor_file_ << "cutoff               double     length"
    "          []" << endl;
  descriptor_file_ << "energy               double     energy"
    "          []" << endl;
  descriptor_file_ << "forces               double     force"
    "          [numberOfParticles,3]" << endl;
  descriptor_file_ <<
    "#\n"
    "# END OF KIM DESCRIPTOR FILE\n"
    "#\n"
             << endl;
}

void CBK_KIM::UpdateCoordinatesAndKIMValues()
{
  int status;
  double X[3];
  double Influencedist[3], tmp;
  int p, q, i;
  int Top[3], Bottom[3], CurrentInfluenceDist;

  CBK_->InfluenceRegion(Influencedist);
  for (i = 0; i < 3; i++)
  {
    Influencedist[i] *= (cutoff_);
  }

  tmp = 1;
  for (p = 0; p < 3; p++)
  {
    // set influence distance based on cube size
    //
    // also setup to be large enough to encompass Eulerian sphere
    CurrentInfluenceDist = int(ceil(Influencedist[p]));
    tmp *= 2.0 * CurrentInfluenceDist;

    Top[p] = CurrentInfluenceDist;
    Bottom[p] = -CurrentInfluenceDist;
  }

  // copy the particle species for the Internal Atoms
  int *speciesTemp = new int[InternalAtoms_];
  for(int i = 0 ; i < InternalAtoms_; i ++)
    speciesTemp[i] = particleSpecies_[i];

  // clear the containers to be sure we start with fresh lists
  neigh_.numNeigh_.clear();
  neigh_.nListAtom_.clear();
  coordinates_.clear();
  forces_.clear();
  particleMap_.clear();
  particleSpecies_.clear();

  // Get the coordinates and map of the Internal atoms
  Vector coordsTemp = CBK_->CBKtoCoords();
  for (int i = 0 ; i < InternalAtoms_ ; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      coordinates_.push_back(coordsTemp[i*3 + j]);
    }
    particleMap_.push_back(i);
  }

  // set up the corrdinates and the mapping
  p = 0;
  for (q = 0; q < InternalAtoms_; q++)
  {
    for (X[0] = Bottom[0]; X[0] <= Top[0]; X[0]++)
    {
      for (X[1] = Bottom[1]; X[1] <= Top[1]; X[1]++)
      {
        for (X[2] = Bottom[2]; X[2] <= Top[2]; X[2]++)
        {
          if (!(X[0] == 0 && X[1] == 0 && X[2] == 0))
          {
            for (i = 0; i < 3; i++)
            {
              coordinates_.push_back(CBK_->Dx(X, p, q, i) + coordinates_[i]);
            }
            particleMap_.push_back(q);
          }
        }
      }
    }
  }
  numberOfParticles_ = particleMap_.size();
  forces_.resize(numberOfParticles_*3, 0.0);

  // map the particle species
  particleSpecies_.resize(numberOfParticles_);
  for(int i = 0; i < numberOfParticles_; i ++)
    particleSpecies_[i] = speciesTemp[particleMap_[i]];

  delete [] speciesTemp;
  // construct neighbor lists according to NEIGH_PURE_F
  double dx[3];
  int const cutoff2 = (cutoff_)*(cutoff_);

  for(int i = 0; i < numberOfParticles_; i++)
  {
    int a = 0;

    if (i > InternalAtoms_)
    {
      neigh_.numNeigh_.push_back(0);
      continue;
    }

    for(int j = 0; j < numberOfParticles_; j++)
    {
      double r2 = 0;
      for(int k = 0; k < 3; k++)
      {
        dx[k] = coordinates_[j*3 + k] - coordinates_[i*3 + k];
        r2 += dx[k]*dx[k];
      }
      if(r2 <= cutoff2)
      {
        if((i < InternalAtoms_) && (i != j))
        {
          a = a+1;
          neigh_.nListAtom_.push_back(j);
        }
      }
    }
    // part i has a neighbors
    neigh_.numNeigh_.push_back(a);
  }

  KIM_API_setm_data(
      pkim_, &status, 4 * 4,
      "numberOfParticles",                     1,              &numberOfParticles_, 1,
      "particleSpecies",      numberOfParticles_, ((int *) &(particleSpecies_[0])), 1,
      "coordinates",      3 * numberOfParticles_,  ((double *) &(coordinates_[0])), 1,
      "forces",           3 * numberOfParticles_,        ((double*) &(forces_[0])), 1);

  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_data", status);
  }
}

void CBK_KIM::ComputeAndUpdate(int StiffnessYes)
{
  int status;
  StiffnessYes_ = StiffnessYes;

  // Initialize for E1 and E2
  for (int i = 0; i < CBK_->DOFS(); i++)
  {
    ME1_static[i] = 0.0;
    if (StiffnessYes_ == 1)
    {
      for (int j = 0; j < CBK_->DOFS(); j++)
      {
        ME2_static[i][j] = 0.0;
      }
    }
  }
  for (int i = 0; i < CBK_F_->DOFS(); i++)
  {
    ME1_F_static[i] = 0.0;
    if (StiffnessYes_ == 1)
    {
      for (int j = 0; j < CBK_F_->DOFS(); j++)
      {
        ME2_F_static[i][j] = 0.0;
      }
    }
  }

  UpdateCoordinatesAndKIMValues();

  // Make sure the correct compute flags are set
  KIM_API_setm_compute(pkim_, &status, 4*3,
    "process_d2Edr2", StiffnessYes_, 1,
    "process_dEdr", KIM_COMPUTE_TRUE, 1,
    "forces", KIM_COMPUTE_TRUE, 1,
    "energy", KIM_COMPUTE_TRUE, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_compute", status);
  }

  status = KIM_API_model_compute(pkim_);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_compute", status);
  }

  // clear out the old BodyForces, and then combine the calcualted forces for periodic BC's
  for (int i = 0; i < InternalAtoms_; i ++)
  {
    for (int j = 0 ; j < 3; j ++)
    {
      BodyForce_[i][j] = forces_[i * DIM3 + j];
    }
  }
  for (int i = InternalAtoms_; i < (numberOfParticles_); i++)
  {
    for (int j = 0; j < DIM3; j++)
    {
      BodyForce_[(particleMap_[i])][j] += forces_[i * DIM3 + j];
    }
  }
}

int CBK_KIM::get_neigh(void *kimmdl, int *mode, int *request, int *part,
               int *numnei, int **nei1part, double **pRij)
{
  int ier;
  int N;
  int partToReturn;

  intptr_t* pkim = *((intptr_t**) kimmdl);
  neighObject* neigh = (neighObject*) KIM_API_get_data(pkim, "neighObject", &ier);

  // unpack number of particles
  int* numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles",
                                                   &ier);
  if (ier < KIM_STATUS_OK)
  {
    KIM_API_report_error(__LINE__,(char*) __FILE__, (char*) "get_data", ier);
    return ier;
  }
  N = *numberOfParticles;
  // check mode and request
  if(*mode == 0)
  {
    if(*request == 0)
    {
      // reset iterator
      neigh->resetIterator();
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    }
    else if(*request == 1)
    {
      // incriment iterator
      neigh->incrementIterator();
      if(neigh->getIterator() > N)
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      else
      {
        partToReturn = neigh->getIterator() - 1;
      }
    }
    else
    {
      ier = KIM_STATUS_NEIGH_INVALID_REQUEST;
      KIM_API_report_error(__LINE__,(char*) __FILE__,(char*) "Invalid request in get_neigh", ier);
      return ier;
    }
  }
  else if(*mode == 1)
  {
    if( *request >= N || *request < 0)
    {
      ier = KIM_STATUS_PARTICLE_INVALID_ID;
      KIM_API_report_error(__LINE__,(char*) __FILE__,(char*) "Invalid part ID in get_neigh", ier);
      return ier;
    }
    else
      partToReturn = *request;
  }
  else
  {
    ier = KIM_STATUS_NEIGH_INVALID_MODE;
    KIM_API_report_error(__LINE__,(char*) __FILE__,(char*) "Invalid mode in get_neigh", ier);
    return ier;
  }

  // set the returned part
  *part = partToReturn;

  // set the returned number of neighbors for the returned part
  *numnei = neigh->numNeigh_[(*part)];

  // set the location of the returned neighbor list
  int temp = 0;
  if (*part != 0)
  {
    for (int i = 1; i <= (*part); i++)
    {
      temp += neigh->numNeigh_[i - 1];
    }
  }
  *nei1part = ((int*) &(neigh->nListAtom_[temp]));

  pRij = NULL;

  return KIM_STATUS_OK;
}

int CBK_KIM::process_dEdr(void* kimmdl, double* dEdr, double* r,
                                  double** dx, int* iSpec, int* jSpec)
{
  int status;
  intptr_t* pkim = *((intptr_t**) kimmdl);
  CBK_KIM* obj;
  obj = (CBK_KIM*) KIM_API_get_sim_buffer(pkim, &status);
  double DX[DIM3];
  int i = obj->particleMap_[*iSpec];
  int j = obj->particleMap_[*jSpec];


  // @@ need to find a more efficient way to do this...
  Matrix InverseF = (obj->CBK_->F()).Inverse();

  double temp1, temp2;
  for (int rows = 0; rows < DIM3; rows++)
  {
    DX[rows] = 0.0;
    for (int cols = 0; cols < DIM3; cols++)
    {
      DX[rows] += InverseF[rows][cols] * (*dx)[cols];
    }
  }

  // dEdUij or dEdFij
  for (int rows = 0; rows < DIM3; rows++)
  {
    for (int cols = 0; cols < DIM3; cols++)
    {
      obj->ME1_static[(obj->CBK_->INDF(rows, cols))]
          += *dEdr * (obj->CBK_->DyDF(*dx, DX, rows, cols)) / (2.0 * (*r));
      obj->ME1_F_static[(obj->CBK_F_->INDF(rows, cols))]
          += *dEdr * (obj->CBK_F_->DyDF(*dx, DX, rows, cols)) / (2.0 * (*r));
    }
  }

  // dEdSij
  for (int atom = 0; atom < (obj->CBK_->InternalAtoms()); atom++)
  {
    for (int k = 0; k < DIM3; k++)
    {
      obj->ME1_static[(obj->CBK_->INDS(atom, k))]
          += *dEdr * (obj->CBK_->DyDS(*dx, i, j, atom, k)) / (2.0 * (*r));
      obj->ME1_F_static[(obj->CBK_F_->INDS(atom, k))]
          += *dEdr * (obj->CBK_F_->DyDS(*dx, i, j, atom, k)) / (2.0 * (*r));
    }
  }
  if ((obj->StiffnessYes_) == 1)
  {
    double DyDF[DIM3][DIM3];
    int i1, j1, k1, l1;

    // Upper Diag Block (CBK_->Fsize(),CBK_->Fsize())
    for (i1 = 0; i1 < DIM3; i1++)
    {
      for (j1 = 0; j1 < DIM3; j1++)
      {
        for (k1 = 0; k1 < DIM3; k1++)
        {
          for (l1 = 0; l1 < DIM3; l1++)
          {
            obj->ME2_static[obj->CBK_->INDF(i1, j1)][obj->CBK_->INDF(k1, l1)]
                += (*dEdr) * ((0.5 / (*r))
                              * (obj->CBK_->D2yDFF(DX, i1, j1, k1, l1))
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_->DyDF((*dx), DX, i1, j1))
                              * (obj->CBK_->DyDF((*dx), DX, k1, l1)));
            obj->ME2_F_static[obj->CBK_F_->INDF(i1, j1)][obj->CBK_F_->INDF(k1, l1)]
                += (*dEdr) * ((0.5 / (*r))
                              * (obj->CBK_F_->D2yDFF(DX, i1, j1, k1, l1))
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_F_->DyDF((*dx), DX, i1, j1))
                              * (obj->CBK_F_->DyDF((*dx), DX, k1, l1)));
          }
        }
      }
    }

    // Lower Diagonal blocks
    for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
    {
      for (int atom1 = 0; atom1 < (obj->CBK_->InternalAtoms()); atom1++)
      {
        for (k1 = 0; k1 < DIM3; k1++)
        {
          for (l1 = 0; l1 < DIM3; l1++)
          {
            obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDS(atom1, l1)]
                += (*dEdr) * ((0.5 / (*r))
                              * obj->CBK_->D2yDSS(i, j, atom0, k1, atom1, l1)
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_->DyDS(*dx, i, j, atom0, k1))
                              * (obj->CBK_->DyDS(*dx, i, j, atom1, l1)));
            obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDS(atom1, l1)]
                += (*dEdr) * ((0.5 / (*r))
                              * obj->CBK_F_->D2yDSS(i, j, atom0, k1, atom1, l1)
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_F_->DyDS(*dx, i, j, atom0, k1))
                              * (obj->CBK_F_->DyDS(*dx, i, j, atom1, l1)));
          }
        }
      }
    }

    // Off-diagonal blocks
    for (i1 = 0; i1 < DIM3; i1++)
    {
      for (j1 = 0; j1 < DIM3; j1++)
      {
        for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
        {
          for (k1 = 0; k1 < DIM3; k1++)
          {
            double temp=(*dEdr)
                * ((0.5 / (*r))
                   * obj->CBK_->D2yDFS(*dx, DX, i, j, i1, j1, atom0, k1)
                   - (0.25 / ((*r) * (*r) * (*r)))
                   * (obj->CBK_->DyDS(*dx, i, j, atom0, k1))
                   * (obj->CBK_->DyDF((*dx), DX, i1, j1)));

            obj->ME2_static[obj->CBK_->INDF(i1, j1)][obj->CBK_->INDS(atom0, k1)]
                += temp;

            obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDF(i1, j1)]
                += temp;

            double Ftemp=(*dEdr)
                * ((0.5 / (*r))
                   * obj->CBK_F_->D2yDFS(*dx, DX, i, j, i1, j1, atom0, k1)
                   - (0.25 / ((*r) * (*r) * (*r)))
                   * (obj->CBK_F_->DyDS(*dx, i, j, atom0, k1))
                   * (obj->CBK_F_->DyDF((*dx), DX, i1, j1)));

            obj->ME2_F_static[obj->CBK_F_->INDF(i1, j1)][obj->CBK_F_->INDS(atom0, k1)]
                += temp;

            obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDF(i1, j1)]
                += temp;
          }
        }
      }
    }
  }

  return KIM_STATUS_OK;
}

int CBK_KIM::process_d2Edr2(void* kimmdl, double* d2Edr2, double** r,
                                    double** dx, int** iSpecs, int** jSpecs)
{
  int status;
  intptr_t* pkim = *((intptr_t**) kimmdl);
  CBK_KIM* obj;
  obj = (CBK_KIM*) KIM_API_get_sim_buffer(pkim, &status);

  // map periodic padding particle
  int const i[2] = {obj->particleMap_[(*iSpecs)[0]], obj->particleMap_[(*iSpecs)[1]]};
  int const j[2] = {obj->particleMap_[(*jSpecs)[0]], obj->particleMap_[(*jSpecs)[1]]};

  double DX[2][DIM3];
  double Dx[2][DIM3];

  // @@ need to find a more efficient way to do this...
  Matrix InverseF = (obj->CBK_->F()).Inverse();

  for (int k = 0; k < DIM3; k++)
  {
    Dx[0][k] = (*dx)[k];
    Dx[1][k] = (*dx)[k + DIM3];
  }

  for (int atoms = 0; atoms < 2; atoms++)
  {
    for (int rows = 0; rows < DIM3; rows++)
    {
      DX[atoms][rows] = 0.0;
      for (int cols = 0; cols < DIM3; cols++)
      {
        DX[atoms][rows] += InverseF[rows][cols] * Dx[atoms][cols];
      }
    }
  }
  int i1, j1, k1, l1;

  // Upper Diagonal
  for (i1 = 0; i1 < DIM3; i1++)
  {
    for (j1 = 0; j1 < DIM3; j1++)
    {
      for (k1 = 0; k1 < DIM3; k1++)
      {
        for (l1 = 0; l1 < DIM3; l1++)
        {
          obj->ME2_static[obj->CBK_->INDF(i1, j1)][obj->CBK_->INDF(k1, l1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1])) * (*d2Edr2)
              * (obj->CBK_->DyDF(Dx[0], DX[0], i1, j1))
              * (obj->CBK_->DyDF(Dx[1], DX[1], k1, l1));
          obj->ME2_F_static[obj->CBK_F_->INDF(i1, j1)][obj->CBK_F_->INDF(k1, l1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1])) * (*d2Edr2)
              * (obj->CBK_F_->DyDF(Dx[0], DX[0], i1, j1))
              * (obj->CBK_F_->DyDF(Dx[1], DX[1], k1, l1));
        }
      }
    }
  }

  // Lower Diagonal blocks
  for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
  {
    // @@ This can be made more efficient
    if (((atom0 == i[0]) || (atom0 == j[0])))
    {
      for (int atom1 = 0; atom1 < (obj->CBK_->InternalAtoms()); atom1++)
      {
        if (((atom1 == i[1]) || (atom1 == j[1])))
        {
          for (k1 = 0; k1 < DIM3; k1++)
          {
            for (l1 = 0; l1 < DIM3; l1++)
            {
              obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDS(atom1, l1)]
                  += (0.5 / (*r)[0]) * (0.5 / (*r)[1]) * (*d2Edr2)
                  * (obj->CBK_->DyDS(Dx[0], i[0], j[0], atom0, k1))
                  * (obj->CBK_->DyDS(Dx[1], i[1], j[1], atom1, l1));
              obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDS(atom1, l1)]
                  += (0.5 / (*r)[0]) * (0.5 / (*r)[1]) * (*d2Edr2)
                  * (obj->CBK_F_->DyDS(Dx[0], i[0], j[0], atom0, k1))
                  * (obj->CBK_F_->DyDS(Dx[1], i[1], j[1], atom1, l1));
            }
          }
        }
      }
    }
  }

  // Off Diagonal
  for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
  {
    for (i1 = 0; i1 < DIM3; i1++)
    {
      for (j1 = 0; j1 < DIM3; j1++)
      {
        for (k1 = 0; k1 < DIM3; k1++)
        {
          obj->ME2_static[obj->CBK_->INDF(i1, j1)][obj->CBK_->INDS(atom0, k1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_->DyDF(Dx[0], DX[0], i1, j1)
                 * obj->CBK_->DyDS(Dx[1], i[1], j[1], atom0, k1));

          obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDF(i1, j1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_->DyDF(Dx[1], DX[1], i1, j1)
                 * obj->CBK_->DyDS(Dx[0], i[0], j[0], atom0, k1));

          obj->ME2_F_static[obj->CBK_F_->INDF(i1, j1)][obj->CBK_F_->INDS(atom0, k1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_F_->DyDF(Dx[0], DX[0], i1, j1)
                 * obj->CBK_F_->DyDS(Dx[1], i[1], j[1], atom0, k1));

          obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDF(i1, j1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_F_->DyDF(Dx[1], DX[1], i1, j1)
                 * obj->CBK_F_->DyDS(Dx[0], i[0], j[0], atom0, k1));
        }
      }
    }
  }

  return KIM_STATUS_OK;
}