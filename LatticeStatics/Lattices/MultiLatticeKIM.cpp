#include "MultiLatticeKIM.h"
#include "UtilityFunctions.h"
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>

using namespace std;

int const MultiLatticeKIM::DIM3 = 3;

double const RoEig_[MultiLatticeKIM::DIM3] = {1.0, 2.0, 3.0};
double const TrEig_[MultiLatticeKIM::DIM3] = {4.0, 5.0, 6.0};



MultiLatticeKIM::~MultiLatticeKIM()
{
  delete[] BodyForce_;
  delete CBK_;
  int status;

  status = KIM_API_model_destroy(pkim_);
  KIM_API_free(&pkim_, &status);
  // delete memory
  delete[] particleSpecies_;
  delete[] coords_;
  delete[] forces_;

}

MultiLatticeKIM::MultiLatticeKIM(PerlInput const& Input, int const& Echo = 1,
                                 int const& Width = 20) :
    Lattice(Input, Echo)
{
  // Get Lattice definition
  PerlInput::HashStruct Hash = Input.getHash("Lattice", "MultiLatticeKIM");
  // Set default values
  KillTranslations_ = 1; // 1-true, 0-false
  StiffnessYes_ = 0;
  int needKillRotations = 1;
  KillRotations_ = 0; // 0-do nothing, 1-kill one rotation, 2-kill 3 rotations
  if (Input.ParameterOK(Hash, "FastPrint"))
  {
    const char* FastPrnt = Input.getString(Hash, "FastPrint");
    if ((!strcmp("Yes", FastPrnt)) || (!strcmp("yes", FastPrnt)))
    {
      FastPrint_ = 1;
    }
    else
    {
      FastPrint_ = 0;
    }
  }
  else
  {
    FastPrint_ = 0;
    Input.useString("No", Hash, "FastPrint");
  }

  PerlInput::HashStruct CBKHash = Input.getHash(Hash, "CBKinematics");
  const char* CBKin = Input.getString(CBKHash, "Type");
  if (!strcmp("SymLagrangeWTransCB", CBKin))
  {
    CBK_ = new SymLagrangeWTransCB(Input, &Hash);
    CBK_F_ = new LagrangeCB(Input, &Hash);
    needKillRotations = 0;
  }
  else if (!strcmp("LagrangeCB", CBKin))
  {
    CBK_ = new LagrangeCB(Input, &Hash);
    CBK_F_ = CBK_;
  }
  else
  {
    cerr << "Error Unknown/unsupported MultiLatticeKIM{CBKinematics}{Type} "
         << "specified\n";
    exit(-9);
  }
  InternalAtoms_ = CBK_->InternalAtoms();

  // Update KillRotations_ if needed
  if (needKillRotations)
  {
    Vector R(DIM3), r(DIM3);
    double norm;
    const char* KillRot = Input.getString(Hash, "RotationConstraint", 0);
    if (!strcmp("FullRotationConstraint", KillRot))
    {
      KillRotations_ = 2;
    }
    else if (!strcmp("OneRotationConstraint", KillRot))
    {
      KillRotations_ = 1;
      Input.getVector(r, Hash, "RotationConstraint", 1);
      Input.getVector(R, Hash, "RotationConstraint", 2);
      KillOneRotation_.Resize(CBK_->DOFS(), 0.0);
      for (int i = 0; i < DIM3; ++i)
      {
        for (int j = 0; j < DIM3; ++j)
        {
          KillOneRotation_[CBK_->INDF(i, j)] = r[i] * R[j];
        }
      }

      norm = KillOneRotation_.Norm();
      for (int i = 0; i < KillOneRotation_.Dim(); ++i)
      {
        KillOneRotation_[i] /= norm;
      }
    }
    else if (!strcmp("NoRotationConstraint", KillRot))
    {
      KillRotations_ = 0;
      for (int i=0 ; i<3; i++)
      {
        Rsq_static[i] = 0.0;
      }
    }
    else
    {
      cerr << "Error (MultiLatticeKIM()): Unknown RotationConstraint type"
           << "\n";
      exit(-2);
    }
  }
  else
  {
    for (int i=0 ; i<3; i++)
    {
      Rsq_static[i] = 0.0;
    }
  }

  // Setup Bodyforce_
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

  numberOfParticles_ = InternalAtoms_; // From CBKHash

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

  particleSpecies_ = new int[InternalAtoms_];

  coords_ = new double[InternalAtoms_ * 3];
  forces_ = new double[InternalAtoms_ * 3];
  KIM_API_setm_data(
      pkim_, &status, 7 * 4,
      "numberOfParticles", 1, &numberOfParticles_, 1,
      "numberOfSpecies", 1, &numberOfSpecies_, 1,
      "particleSpecies", InternalAtoms_, &(particleSpecies_[0]), 1,
      "coordinates", 3 * InternalAtoms_, &(coords_[0]), 1,
      "forces", 3 * InternalAtoms_, &(forces_[0]), 1,
      "cutoff", 1, &cutoff_, 1,
      "energy", 1, &energy_, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_data", status);
  }
  KIM_API_setm_method(
      pkim_, & status, 3 * 4,
      "get_neigh", 1, (func_ptr) &MultiLatticeKIM::get_neigh, 1,
      "process_dEdr", 1, (func_ptr) &MultiLatticeKIM::process_dEdr, 1,
      "process_d2Edr2", 1, (func_ptr) &MultiLatticeKIM::process_d2Edr2, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_method", status);
  }

  status = KIM_API_model_init(pkim_);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_model_init", status);
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

  InfluenceDist_ = cutoff_;
  LatSum_(CBK_, InternalAtoms_, &InfluenceDist_);

  KIM_API_setm_data(pkim_, &status, 1 * 4,
                    "neighObject", 1, &LatSum_, 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_setm_data", status);
  }
  // figures out which atoms corresponds to each element in the species list
  for (int i = 0; i < InternalAtoms_; i++)
  {
    AtomSpeciesList[i] = Input.getString(CBKHash, "AtomSpeciesKIM", i);
    for (int j = 0; j < numberOfSpecies_; j++)
    {
      if (!strcmp(AtomSpeciesList[i], SpeciesList[j]))
      {
        particleSpecies_[i] = KIM_API_get_species_code(
            pkim_, (char*) SpeciesList[j], &status);
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

  // Get Lattice parameters
  if (Input.ParameterOK(Hash, "Density"))
  {
    Density_ = Input.getInt(Hash, "Density");
  }
  else
  {
    Density_ = Input.useInt(1, Hash, "Density"); // Default Value
  }
  ConvexityDX_ = Input.getDouble(Hash, "ConvexityDX");
  // Set Loading parameter
  LoadParameter_ = Load;

  Lambda_ = 0.0;
  Matrix Tractions(DIM3,DIM3);
  Matrix Normals(DIM3,DIM3);
  // Tractions = [ t^(1) | t^(2) | t^(3) ]
  Input.getMatrix(Tractions, Hash, "LoadTractions");
  // Normals = [ N^(1) | N^(2) | N^(3) ]
  Input.getMatrix(Normals, Hash, "LoadNormals");
  Loading_.Resize(DIM3, DIM3, 0.0);
  Loading_ = Tractions*(Normals.Inverse());

  // Check that Loading_ is symmetric
  //
  // Make the restriction that for 1st PK loading (with full F kinematics),
  // the reference config. is the equilibrium config
  // obtained in the limit as lambda->0.
  Matrix Dev = (Loading_ - Loading_.Transpose());
  if ((abs(Dev.MaxElement()) > 1.0e-12) ||
      (abs(Dev.MinElement()) > 1.0e-12))
  {
    cerr << "Error (MultiLatticeKIM()): Unsymmetric Loading_ found; "
        "Symmetric Loading_ required!" << "\n";
    exit(-9);
  }

  // needed to initialize reference length
  int iter;
  iter = Input.getPosInt(Hash, "MaxIterations");

  // values to identifiy REFERENCE CONFIGURATION
  if (Input.ParameterOK(Hash, "ReferenceLambda"))
  {
    REFLambda_ = Input.getDouble(Hash, "ReferenceLambda");
  }
  else
  {
    // Default Value
    REFLambda_ = Input.useDouble(0.0, Hash, "ReferenceLambda");
  }

  PerlInput::HashStruct TFHash = Input.getHash(Hash, "ExtraTestFunctions");
  const char* TFtyp = Input.getString(TFHash, "Type");
  if ((!strcmp("None", TFtyp)) || (!strcmp("none", TFtyp)))
  {
    TFType_ = 0;
    NumExtraTFs_ = 0;
  }
  else if ((!strcmp("LoadingParameters", TFtyp))
           || (!strcmp("loadingparameters", TFtyp)))
  {
    TFType_ = 1;
    NumExtraTFs_ = Input.getArrayLength(TFHash, "LoadingParameters");

    TFLoad_.Resize(NumExtraTFs_);
    Input.getVector(TFLoad_, TFHash, "LoadingParameters");
  }
  else
  {
    cerr << "Error (MultiLatticeKIM()): Unknown TestFunctions{Type}" << "\n";
    exit(-3);
  }

  // Initialize various data storage space
  ME1_static.Resize(CBK_->DOFS(), 0.0);
  ME1_F_static.Resize(CBK_F_->DOFS(), 0.0);
  ME2_static.Resize(CBK_->DOFS(), CBK_->DOFS(), 0.0);
  ME2_F_static.Resize(CBK_F_->DOFS(), CBK_F_->DOFS(), 0.0);
  str_static.Resize(CBK_->DOFS());
  stiff_static.Resize(CBK_->DOFS(), CBK_->DOFS());
  CondEV_static.Resize(1, CBK_F_->Fsize());
  CondModuli_static.Resize(CBK_F_->Fsize(), CBK_F_->Fsize());
  OPM_static.Resize(CBK_F_->Ssize(), CBK_F_->Ssize());
  UnitCellShiftModuli_static.Resize(CBK_F_->Ssize(), CBK_F_->Ssize());
  OpticEV_static.Resize(1, CBK_F_->Ssize());
  OpticEV_Print.Resize(CBK_F_->Ssize());
  TestFunctVals_static.Resize(NumTestFunctions());
  if (TFType_ == 1) // only print stiffness eigenvalues
  {
    TestFunctVals_Print.Resize(CBK_->DOFS());
  }
  else // print everything
  {
    TestFunctVals_Print.Resize(TestFunctVals_static.Dim());
  }

  // Initialize Cached_ values
  for (int i=0; i<cachesize; ++i) Cached_[i] = 0;
  if (Input.ParameterOK(Hash, "InitialEqbm"))
  {
    const char* init_equil = Input.getString(Hash, "InitialEqbm");
    if (!strcmp("Yes", init_equil) || !strcmp("yes", init_equil))
    {
      int err = 0;
      if (Input.ParameterOK(Hash, "InitialEqbmCubic"))
      {
        const char* cubic_equil = Input.getString(Hash, "InitialEqbmCubic");
        if (!strcmp("Yes", cubic_equil) || !strcmp("yes", cubic_equil))
        {
          err = FindLatticeSpacing(iter, true);
        }
        else
        {
          err = FindLatticeSpacing(iter, false);
        }
      }
      else
      {
        err = FindLatticeSpacing(iter, false);
      }
      if (err)
      {
        cerr << "unable to find initial lattice spacing!" << "\n";
        exit(-1);
      }
    }
  }
  else
  {
    Input.useString("Yes", Hash, "InitialEqbm");
    int err = 0;
    err = FindLatticeSpacing(iter,false);
    if (err)
    {
      cerr << "unable to find initial lattice spacing!" << "\n";
      exit(-1);
    }
  }

  // Setup initial status for parameters
  Lambda_ = Input.getDouble(Hash, "Lambda");

  Input.EndofInputSection();
}

int MultiLatticeKIM::FindLatticeSpacing(int const& iter, bool cubicEqbm)
{
  const double Tol = DOF().Dim() * 1.0e-13;

  Lambda_ = REFLambda_;

  CBK_->SetReferenceDOFs();
  if (CBK_F_ != CBK_) CBK_F_->SetReferenceDOFs();
  LatSum_.Recalc();

  if (cubicEqbm)
  {
    if (Echo_)
    {
      RefineCubicEqbm(Tol, iter, &cout);
    }
    else
    {
      RefineCubicEqbm(Tol, iter, 0);
    }
  }
  else
  {
    if (Echo_)
    {
      RefineEqbm(Tol, iter, &cout);
    }
    else
    {
      RefineEqbm(Tol, iter, 0);
    }
  }

  // Clean up numerical round off (at least for zero values)
  Vector doftmp(CBK_->DOFS(), 0.0);
  for (int i = 0; i < CBK_->DOFS(); ++i)
  {
    if (fabs(CBK_->DOF()[i]) < Tol)
    {
      doftmp[i] = 0.0;
    }
    else
    {
      doftmp[i] = CBK_->DOF()[i];
    }
  }
  CBK_->SetDOF(doftmp);

  CBK_->SetReferenceToCurrent();
  if (CBK_F_ != CBK_) CBK_F_->SetReferenceToCurrent();

  LatSum_.Recalc();

  return 0;
}

void MultiLatticeKIM::SetParameters(double const* const Vals,
                                    int const& ResetRef = 1)
{
  cout << "ERROR IN MultiLatticeKIM::SetParameters. Exiting" << endl;
  exit(1);
}

void MultiLatticeKIM::UpdateKIMValues() const
{
  int status;

  // Initialize for E1 and E2
  for (int i = 0; i < CBK_->DOFS(); i++)
  {
    ME1_static[i] = 0.0;
    if (StiffnessYes_==1)
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
    if (StiffnessYes_==1)
    {
      for (int j = 0; j < CBK_F_->DOFS(); j++)
      {
        ME2_F_static[i][j] = 0.0;
      }
    }
  }

  Vector coordsTemp = CBK_->CBKtoCoords();
  for (int i = 0; i < (DIM3 * InternalAtoms_); i++)
  {
    coords_[i] = coordsTemp[i];
  }

  LatSum_.Recalc();

  // Make sure the correct compute flags are set
  KIM_API_set_compute(pkim_, "process_d2Edr2", StiffnessYes_, &status);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_set_compute", status);
  }

  status = KIM_API_model_compute(pkim_);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "KIM_API_compute", status);
  }

  for (int i = 0; i < (InternalAtoms_); i++)
  {
    for (int j = 0; j < DIM3; j++)
    {
      BodyForce_[i][j] = forces_[i * DIM3 + j];
    }
  }
}

// Lattice Routines
double MultiLatticeKIM::E0() const
{
  if (!Cached_[0])
  {
    E0CachedValue_ = energy();

    if (KillTranslations_)
    {
      for (int j = 0; j < 3; ++j)
      {
        Tsq_static[j] = 0.0;
        for (int i = 0; i < InternalAtoms_; ++i)
        {
          Tsq_static[j] += CBK_->DOF()[CBK_->INDS(i, j)];
        }
        Tsq_static[j] = (Tsq_static[j] * Tsq_static[j]) / InternalAtoms_;
      }
    }

    switch (KillRotations_)
    {
      case 2:
        // Kill three rotations
        Rsq_static[0] = (CBK_->DOF()[CBK_->INDF(0, 1)]
                         - CBK_->DOF()[CBK_->INDF(1, 0)]);
        Rsq_static[1] = (CBK_->DOF()[CBK_->INDF(1, 2)]
                         - CBK_->DOF()[CBK_->INDF(2, 1)]);
        Rsq_static[2] = (CBK_->DOF()[CBK_->INDF(2, 0)]
                         - CBK_->DOF()[CBK_->INDF(0, 2)]);
        for (int i = 0; i < DIM3; ++i)
        {
          Rsq_static[i] *= 0.5 * Rsq_static[i];
        }
        break;
      case 1:
        // Kill one rotation
        for (int i = 0; i < DIM3; ++i)
        {
          Rsq_static[i] = 0.0;
          for (int j = 0; j < DIM3; ++j)
          {
            Rsq_static[0] += KillOneRotation_[CBK_->INDF(i, j)]
                * CBK_->DOF()[CBK_->INDF(i, j)];
          }
        }
        Rsq_static[0] *= Rsq_static[0];
        break;
    }

    E0CachedValue_ += 0.5 * (TrEig_[0] * Tsq_static[0] + TrEig_[1]
                             * Tsq_static[1] + TrEig_[2] * Tsq_static[2])
        + 0.5 * (RoEig_[0] * Rsq_static[0] + RoEig_[1] * Rsq_static[1]
                 + RoEig_[2] * Rsq_static[2]);

    // Update Cached_[0] value
    Cached_[0] = 1;
  }

  return E0CachedValue_;
}


double MultiLatticeKIM::energy(LDeriv const& dl) const
{
  double phi = 0.0;
  double Vr = Density_ ? CBK_->RefVolume() : 1.0;

  if (dl == L0)
  {
    StiffnessYes_ = 0;
    UpdateKIMValues();
    //compute energy per volume
    phi = energy_ / Vr;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        phi -= Lambda_ * Loading_[i][j] * ((CBK_->DOF())[CBK_->INDF(i,j)]
                                           - Del(i,j));
      }
    }
  }
  else
  {
    cerr << "Unknown LDeriv dl in MultiLatticeKIM::energy()" << "\n";
    exit(-1);
  }

  return phi;
}

double MultiLatticeKIM::ConjugateToLambda() const
{
  double conj = 0.0;

  for (int i = 0; i < DIM3; ++i)
  {
    for (int j = 0; j < DIM3; ++j)
    {
      conj += Loading_[i][j] * ((CBK_->DOF())[CBK_->INDF(i, j)] - Del(i, j));
    }
  }

  return conj;
}

Vector const& MultiLatticeKIM::E1() const
{
  if (!Cached_[1])
  {
    stress();
    if (KillTranslations_)
    {
      for (int j = 0; j < DIM3; ++j)
      {
        T_static[j] = 0.0;
        for (int i = 0; i < InternalAtoms_; ++i)
        {
          T_static[j] += CBK_->DOF()[CBK_->INDS(i, j)];
        }
        T_static[j] /= InternalAtoms_;
      }
      for (int i = 0; i < InternalAtoms_; ++i)
      {
        for (int j = 0; j < DIM3; ++j)
        {
          ME1_static[CBK_->INDS(i, j)] += TrEig_[j] * T_static[j];
        }
      }
    }

    switch (KillRotations_)
    {
      case 2:
        // Kill three rotations
        R_static[0] = (CBK_->DOF()[CBK_->INDF(0, 1)]
                       - CBK_->DOF()[CBK_->INDF(1, 0)]) / 2.0;
        ME1_static[CBK_->INDF(0, 1)] += RoEig_[0] * R_static[0];
        ME1_static[CBK_->INDF(1, 0)] -= RoEig_[0] * R_static[0];
        R_static[1] = (CBK_->DOF()[CBK_->INDF(1, 2)]
                       - CBK_->DOF()[CBK_->INDF(2, 1)]) / 2.0;
        ME1_static[CBK_->INDF(1, 2)] += RoEig_[1] * R_static[1];
        ME1_static[CBK_->INDF(2, 1)] -= RoEig_[1] * R_static[1];
        R_static[2] = (CBK_->DOF()[CBK_->INDF(2, 0)]
                       - CBK_->DOF()[CBK_->INDF(0, 2)]) / 2.0;
        ME1_static[CBK_->INDF(2, 0)] += RoEig_[2] * R_static[2];
        ME1_static[CBK_->INDF(0, 2)] -= RoEig_[2] * R_static[2];
        break;
      case 1:
        // Kill one rotation
        for (int i = 0; i < DIM3; ++i)
        {
          R_static[i] = 0.0;
          for (int j = 0; j < DIM3; ++j)
          {
            R_static[0] += KillOneRotation_[CBK_->INDF(i, j)]
                * CBK_->DOF()[CBK_->INDF(i, j)];
          }
        }

        for (int i = 0; i < ME1_static.Dim(); ++i)
        {
          ME1_static[i] += RoEig_[0] * R_static[0] * KillOneRotation_[i];
        }
        break;
    }
    // Update Cached_[1] value
    Cached_[1] = 1;
  }

  return ME1_static;
}

Vector const& MultiLatticeKIM::stress(LDeriv const& dl) const
{
  double Vr = Density_ ? CBK_->RefVolume() : 1.0;

  if (dl == L0)
  {
    StiffnessYes_ = 0;
    UpdateKIMValues();
    ME1_static *= 1.0 / Vr;

    for (int i = 0; i < DIM3; ++i)
    {
      for (int j = 0; j < DIM3; ++j)
      {
        ME1_static[(CBK_->INDF(i, j))] -= Lambda_ * Loading_[i][j];
      }
    }
  }
  else if (dl == DL)
  {
    for (int i = 0; i<ME1_static.Dim(); ++i)
    {
      ME1_static[i] = 0.0;
    }

    for (int i = 0; i < DIM3; ++i)
    {
      for (int j = 0; j < DIM3; ++j)
      {
        ME1_static[(CBK_->INDF(i, j))] -= Loading_[i][j];
      }
    }
  }
  else
  {
    cerr << "Unknown LDeriv dl in MultiLatticeKIM::stress()" << "\n";
    exit(-1);
  }
  return ME1_static;
}

Matrix const& MultiLatticeKIM::E2() const
{
  if (!Cached_[2])
  {
    stiffness();
    if (KillTranslations_)
    {
      for (int i = 0; i < InternalAtoms_; ++i)
      {
        for (int j = 0; j < DIM3; ++j)
        {
          for (int k = 0; k < InternalAtoms_; ++k)
          {
            ME2_static[CBK_->INDS(i, j)][CBK_->INDS(k, j)]
                += TrEig_[j] / InternalAtoms_;
          }
        }
      }
    }

    switch (KillRotations_)
    {
      case 2:
        // Kill three rotations
        ME2_static[CBK_->INDF(0, 1)][CBK_->INDF(0, 1)] += RoEig_[0] / 2.0;
        ME2_static[CBK_->INDF(0, 1)][CBK_->INDF(1, 0)] -= RoEig_[0] / 2.0;
        ME2_static[CBK_->INDF(1, 0)][CBK_->INDF(1, 0)] += RoEig_[0] / 2.0;
        ME2_static[CBK_->INDF(1, 0)][CBK_->INDF(0, 1)] -= RoEig_[0] / 2.0;

        ME2_static[CBK_->INDF(1, 2)][CBK_->INDF(1, 2)] += RoEig_[1] / 2.0;
        ME2_static[CBK_->INDF(1, 2)][CBK_->INDF(2, 1)] -= RoEig_[1] / 2.0;
        ME2_static[CBK_->INDF(2, 1)][CBK_->INDF(2, 1)] += RoEig_[1] / 2.0;
        ME2_static[CBK_->INDF(2, 1)][CBK_->INDF(1, 2)] -= RoEig_[1] / 2.0;

        ME2_static[CBK_->INDF(2, 0)][CBK_->INDF(2, 0)] += RoEig_[2] / 2.0;
        ME2_static[CBK_->INDF(2, 0)][CBK_->INDF(0, 2)] -= RoEig_[2] / 2.0;
        ME2_static[CBK_->INDF(0, 2)][CBK_->INDF(0, 2)] += RoEig_[2] / 2.0;
        ME2_static[CBK_->INDF(0, 2)][CBK_->INDF(2, 0)] -= RoEig_[2] / 2.0;
        break;
      case 1:
        // Kill one rotation
        for (int i = 0; i < ME2_static.Rows(); ++i)
        {
          for (int j = 0; j < ME2_static.Cols(); ++j)
          {
            ME2_static[i][j] += RoEig_[0] * KillOneRotation_[i]
                * KillOneRotation_[j];
          }
        }
        break;
    }

    // Update Cached_[2] value
    Cached_[2] = 1;
  }

  return ME2_static;
}

Matrix const& MultiLatticeKIM::stiffness(LDeriv const& dl) const
{
  double Vr = Density_ ? CBK_->RefVolume() : 1.0;

  if (dl == L0)
  {
    StiffnessYes_ = 1;
    UpdateKIMValues();
    ME2_static *= 1.0 / Vr;
    ME2_F_static *= 1.0 /Vr;
  }
  else if (dl == DL)
  {
    for (int i = 0; i < ME2_static.Rows(); ++i)
    {
      for (int j = 0; j < ME2_static.Cols(); ++i)
      {
        ME2_static[i][j] = 0.0;
      }
    }
    for (int i = 0; i < ME2_F_static.Rows(); ++i)
    {
      for (int j = 0; j < ME2_F_static.Cols(); ++i)
      {
        ME2_F_static[i][j] = 0.0;
      }
    }
  }
  else
  {
    cerr << "Unknown LDeriv dl in MultiLatticeTKIM::stiffness()" << "\n";
    exit(-1);
  }

  return ME2_static;
}

Matrix const& MultiLatticeKIM::E3() const
{
  cerr << "Error in MultiLatticeKIM::E3() empty function \n ";
  exit(-45);
}

Matrix const& MultiLatticeKIM::E4() const
{
  cerr << "Error in MultiLatticeKIM::E4() empty function \n ";
  exit(-45);
}

int MultiLatticeKIM::CriticalPointInfo(int* const CPCrossingNum,
                                       int const& TFIndex, Vector const& DrDt,
                                       int const& CPorBif,
                                       int const& NumZeroEigenVals,
                                       double const& Tolerance,
                                       int const& Width,
                                       PerlInput const& Input, ostream& out)
{
  int Bif;

  // do standard CPInfo stuff and output bfb restart file
  Bif = Lattice::CriticalPointInfo(
      CPCrossingNum, TFIndex, DrDt, CPorBif, NumZeroEigenVals, Tolerance,
      Width, Input, out);

  // @@ provide additional stuff?

  return Bif;
}

void MultiLatticeKIM::ExtraTestFunctions(Vector& TF) const
{
  if (TFType_ == 1) // LoadingParameter
  {
    for (int i = 0; i < NumExtraTFs_; ++i)
    {
      TF[i] = (TFLoad_[i] - Lambda());
    }
  }
}

Matrix const& MultiLatticeKIM::CondensedModuli() const
{
  Matrix const& stiff = E2();  // Make sure ME2_F_static is up to date
  int intrn = CBK_F_->Ssize();
  double factor = 1.0 / (intrn / DIM3);
  int fsz = CBK_F_->Fsize();
  Matrix IM(intrn, intrn);
  CM_static.Resize(fsz, fsz);

  for (int i = 0; i < fsz; i++)
  {
    for (int j = 0; j < fsz; j++)
    {
      CM_static[i][j] = ME2_F_static[i][j];
    }
  }

  // Make sure there are internal DOF's
  if (intrn)
  {
    for (int i = 0; i < intrn; i++)
    {
      for (int j = 0; j < intrn; j++)
      {
        IM[i][j] = ME2_F_static[fsz + i][fsz + j];

        // add translational stiffness to regularize IM, if needed
        if ((!CBK_F_->NoTrans()) && (i % DIM3 == j % DIM3))
        {
          IM[i][j] += factor;
        }
      }
    }
    IM = IM.Inverse();

    // Set up Condensed Moduli
    for (int i = 0; i < fsz; i++)
    {
      for (int j = 0; j < fsz; j++)
      {
        for (int m = 0; m < intrn; m++)
        {
          for (int n = 0; n < intrn; n++)
          {
            CM_static[i][j] -= ME2_F_static[i][fsz + m] * IM[m][n]
                * ME2_F_static[fsz + n][j];
          }
        }
      }
    }
  }

  return CM_static;
}

Matrix const& MultiLatticeKIM::UnitCellShiftModuli() const
{
  Matrix const& stiff = E2();  // Make sure ME2_F_static is up to date
  int fsz = CBK_F_->Fsize();
  int intrn = CBK_F_->Ssize();
  double factor = 1.0 / (intrn / DIM3);

  // Make sure there are internal DOF's
  if (intrn)
  {
    for (int i = 0; i < intrn; i++)
    {
      for (int j = 0; j < intrn; j++)
      {
        OPM_static[i][j] = ME2_F_static[fsz + i][fsz + j];

        // add translational stiffness to regularize IM, if needed
        if ((!CBK_F_->NoTrans()) && (i % DIM3 == j % DIM3))
        {
          OPM_static[i][j] += factor;
        }
      }
    }
  }

  return OPM_static;
}


// @@ provide definition
void MultiLatticeKIM::LongWavelengthModuli(double const& dk,
                                           int const& gridsize,
                                           char const* const prefix,
                                           ostream& out) const
{
  cout << "ERROR IN MultiLatticeKIM::LongWavelengthModuli. Exiting" << endl;
  exit(1);
}

void MultiLatticeKIM::Print(ostream& out, PrintDetail const& flag,
                            PrintPathSolutionType const& SolType = RegularPt)
{
  int W;
  int NoNegTestFunctions = 0;
  double engy;
  int RankOneConvex = 0;
  double minRK1 = 0.0;
  int UnitCellShiftStable = 0;
  double minUnitCellShift = 0.0;
  double mintestfunct;
  double conj;
  int NoFP = !FastPrint_;

  W = out.width();

  out.width(0);
  if (Echo_)
  {
    cout.width(0);
  }
  LatSum_.Recalc();
  LatSum_.Reset();

  engy = E0();
  conj = ConjugateToLambda();
  str_static = stress();

  if (NoFP)
  {
    stiff_static = E2();

    TestFunctions(TestFunctVals_static, LHS);
    mintestfunct = TestFunctVals_static[0];

    if (TFType_ == 1) // LoadingParameters
    {
      for (int i = 0; i < CBK_->DOFS(); ++i)
      {
        TestFunctVals_Print[i] = TestFunctVals_static[i];
        if (TestFunctVals_static[i] < 0.0)
        {
          ++NoNegTestFunctions;
        }
        if (mintestfunct > TestFunctVals_static[i])
        {
          mintestfunct = TestFunctVals_static[i];
        }
      }
    }
    else // None
    {
      for (int i = 0; i < TestFunctVals_static.Dim(); ++i)
      {
        TestFunctVals_Print[i] = TestFunctVals_static[i];
        if (TestFunctVals_static[i] < 0.0)
        {
          ++NoNegTestFunctions;
        }
        if (mintestfunct > TestFunctVals_static[i])
        {
          mintestfunct = TestFunctVals_static[i];
        }
      }
    }

    UnitCellShiftModuli_static = UnitCellShiftModuli();
    OpticEV_static = SymEigVal(UnitCellShiftModuli_static);
    minUnitCellShift = OpticEV_static[0][0];
    for (int i = 0; i < OpticEV_static.Cols(); ++i)
    {
      if (minUnitCellShift > OpticEV_static[0][i])
      {
        minUnitCellShift = OpticEV_static[0][i];
      }
      OpticEV_Print[i] = OpticEV_static[0][i];
    }
    UnitCellShiftStable = (minUnitCellShift > 0.0);

    CondModuli_static = CondensedModuli();
    CondEV_static = SymEigVal(CondModuli_static);
    minRK1 = FullScanRank1Convex3D(CBK_F_, CondModuli_static,
                                   ConvexityDX_);
    RankOneConvex = (minRK1 > 0.0);
  }

  switch (flag)
  {
    case PrintLong:
      out << "MultiLatticeKIM:" << "\n" << "\n";
      out << "Density_ = " << Density_ << "\n";
      out << "Using: " << (*CBK_) << " Kinematics" << "\n";
      out << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
      for (int i = 0; i < InternalAtoms_; ++i)
      {
        out << "Atom_" << i << (i > 9 ? "" : " ") << "          "
            << "          Position : " << setw(W) << CBK_->AtomPositions(i)
            << "\n";
      }
      out << "REFLambda_ : " << setw(W) << REFLambda_ << "\n";
      out << "Influence Distance : " << setw(W) << InfluenceDist_ << "\n";
      out << "Loading : " << setw(W) << Loading_;

      // also send to cout
      if (Echo_)
      {
        cout << "MultiLatticeKIM:" << "\n" << "\n";
        cout << "Density_ = " << Density_ << "\n";
        cout << "Using: " << (*CBK_) << " Kinematics" << "\n";
        cout << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
        for (int i = 0; i < InternalAtoms_; ++i)
        {
          cout << "Atom_" << i << (i > 9 ? "" : " ") << "          "
               << "          Position : " << setw(W)
               << CBK_->AtomPositions(i) << "\n";
        }
        cout << "REFLambda_ : " << setw(W) << REFLambda_ << "\n";
        cout << "Influence Distance : " << setw(W) << InfluenceDist_
             << "\n";
        cout << "Loading : " << setw(W) << Loading_;
      }
      // passthrough to short

    case PrintShort:
      out << "Lambda : " << setw(W) << Lambda_ << "\n"
          << "ConjugateToLambda : " << setw(W) << conj << "\n"
          << "DOF's :" << "\n" << setw(W) << CBK_->DOF() << "\n"
          << "Potential Value (eV) :" << setw(W) << engy << "\n";
      for (int i = 0; i < InternalAtoms_; ++i)
      {
        out << "BodyForce Value (KIM Output) " << i << setw(W)
            << BodyForce_[i] << "\n";
      }
      out << "Stress (eV/A^3):" << "\n" << setw(W) << str_static << "\n";
      if (NoFP)
      {
        out << "\nStiffness (eV/A^3):" << setw(W) << stiff_static
            << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n"
            << setw(W)
            << TestFunctVals_Print << "\n"
            << "Bifurcation Info:" << setw(W) << mintestfunct
            << setw(W) << NoNegTestFunctions << "\n"
            << "Condensed Moduli (eV/A^3):" << setw(W) << CondModuli_static
            << "CondEV Info:" << setw(W) << CondEV_static
            << "Condensed Moduli Rank1Convex:" << setw(W)
            << RankOneConvex << setw(W) << minRK1 << "\n"
            << "UnitCellShift EigenValues:" << setw(W) << OpticEV_Print
            << "\n"
            << "UnitCellShiftStability:" << setw(W) << UnitCellShiftStable
            << setw(W) << minUnitCellShift << "\n";

      }

      if (Echo_)
      {
        cout << "Lambda : " << setw(W) << Lambda_ << "\n"
             << "ConjugateToLambda: " << setw(W) << conj << "\n"
             << "DOF's :" << "\n" << setw(W) << CBK_->DOF() << "\n"
             << "Potential Value (eV):" << setw(W) << engy << "\n";
        for (int i = 0; i < InternalAtoms_; ++i)
        {
          cout << "BodyForce Value (KIM Output) " << i << setw(W)
               << BodyForce_[i] << "\n";
        }
        cout << "Stress (eV/A^3):" << "\n" << setw(W) << str_static
             << "\n";
        if (NoFP)
        {
          cout << "\nStiffness (eV/A^3):" << setw(W) << stiff_static
               << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n"
               << setw(W)
               << TestFunctVals_static << "\n"
               << "Bifurcation Info:" << setw(W) << mintestfunct
               << setw(W) << NoNegTestFunctions << "\n"
               << "Condensed Moduli (eV/A^3):" << setw(W)
               << CondModuli_static
               << "CondEV Info:" << setw(W) << CondEV_static
               << "Condensed Moduli Rank1Convex:" << setw(W)
               << RankOneConvex << setw(W) << minRK1 << "\n"
               << "UnitCellShift EigenValues:" << setw(W) << OpticEV_Print
               << "\n"
               << "UnitCellShiftStability:" << setw(W) << UnitCellShiftStable
               << setw(W) << minUnitCellShift << "\n";

        }
      }
      break;
  }
}

ostream& operator<<(ostream& out, MultiLatticeKIM& A)
{
  A.Print(out, Lattice::PrintShort);
  return out;
}

void MultiLatticeKIM::RefineEqbm(double const& Tol, int const& MaxItr,
                                 ostream* const out)
{
  Vector dx(CBK_->DOFS(), 0.0);
  Vector Stress = E1();
  int itr = 0;
  while ((itr < MaxItr) && Stress.Norm() > Tol)
  {
    ++itr;

#ifdef SOLVE_SVD
    dx = SolveSVD(E2(), Stress, MAXCONDITION, Echo_);
#else
    dx = SolvePLU(E2(), Stress);
#endif

    SetDOF(CBK_->DOF() - dx);

    Stress = E1();

    if (out != 0)
    {
      *out << setw(20) << Stress;

      *out << "\t" << itr << "\tdx "
           << dx.Norm() << "\tstress " << Stress.Norm()
           << "\n";
    }
  }
}

void MultiLatticeKIM::RefineCubicEqbm(double const& Tol, int const& MaxItr,
                                      ostream* const out)
{
  Vector dx(CBK_->DOFS(), 0.0);
  Vector Stress = E1();
  Matrix Stiff = E2();
  double f;
  double k;
  double d;

  int itr = 0;
  f = 0.0;
  k = 0.0;
  for (int i=0; i<DIM3; ++i)
  {
    f += Stress[CBK_->INDF(i,i)];
    for (int j=0; j<DIM3; ++j)
    {
      k += Stiff[CBK_->INDF(i,i)][CBK_->INDF(j,j)];
    }
  }
  while ((itr < MaxItr) && (fabs(f) > Tol))
  {
    ++itr;

    d = f/k;

    for (int i=0; i<DIM3; ++i)
    {
      dx[CBK_->INDF(i,i)] = d;
    }

    SetDOF(CBK_->DOF() - dx);

    Stress = E1();
    Stiff = E2();
    f = 0.0;
    k = 0.0;
    for (int i=0; i<DIM3; ++i)
    {
      f += Stress[CBK_->INDF(i,i)];
      for (int j=0; j<DIM3; ++j)
      {
        k += Stiff[CBK_->INDF(i,i)][CBK_->INDF(j,j)];
      }
    }

    if (out != 0)
    {
      *out << setw(20) << Stress;

      *out << "\t" << itr << "\tdx "
           << dx.Norm() << "\tstress " << Stress.Norm()
           << "\n";
    }
  }
}

void MultiLatticeKIM::Write_KIM_descriptor_file(const char** SpeciesList,
                                                int numberOfSpecies_)
{
  descriptor_file_ << "#\n"
      "# BEGINNING OF KIM DESCRIPTOR FILE\n"
      "#\n"
      "# This file is automatically generated from MultiLatticeKIM.\n"
      "#################################################" << endl;
  descriptor_file_ << "KIM_API_Version  := 1.6.0" << endl;
  descriptor_file_ << "Unit_length      := A" << endl;
  descriptor_file_ << "Unit_energy      := eV" << endl;
  descriptor_file_ << "Unit_charge      := e" << endl;
  descriptor_file_ << "Unit_temperature := K" << endl;
  descriptor_file_ << "Unit_time        := ps" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "SUPPORTED_ATOM/PARTICLES_TYPES:" << endl;
  for (int i = 0; i < numberOfSpecies_; i++)
  {
    descriptor_file_ << SpeciesList[i] << " spec  1" << endl;
  }

  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "CONVENTIONS:" << endl;
  descriptor_file_ << "# Name                      Type" << endl;
  descriptor_file_ << "ZeroBasedLists              flag" << endl;
  descriptor_file_ << "Neigh_BothAccess            flag" << endl;
  //      descriptor_file_ << "CLUSTER                    flag" << endl;
  descriptor_file_ << "NEIGH_RVEC_F                flag" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "MODEL_INPUT:" << endl;
  descriptor_file_ << "# Name                      Type         Unit"
      "                Shape" << endl;
  descriptor_file_ << "numberOfParticles           integer      none"
      "                []" << endl;
  descriptor_file_ << "numberOfSpecies             integer      none"
      "                []" << endl;
  descriptor_file_ << "particleSpecies             integer      none"
      "                [numberOfParticles]" << endl;
  descriptor_file_ << "coordinates                 double       length"
      "              [numberOfParticles,3]" << endl;
  descriptor_file_ << "get_neigh                   method       none"
      "                []" << endl;
  descriptor_file_ << "neighObject                 pointer      none"
      "                []" << endl;
  descriptor_file_ << "process_dEdr                method       none"
      "                []" << endl;
  descriptor_file_ << "process_d2Edr2              method       none"
      "                []" << endl;
  descriptor_file_ << "#################################" << endl;
  descriptor_file_ << "MODEL_OUTPUT:" << endl;
  descriptor_file_ << "# Name                      Type         Unit"
      "                Shape" << endl;
  descriptor_file_ << "destroy                     method       none"
      "                []" << endl;
  descriptor_file_ << "compute                     method       none"
      "                []" << endl;
  descriptor_file_ << "cutoff                      double       length"
      "              []" << endl;
  descriptor_file_ << "energy                      double       energy"
      "              []" << endl;
  descriptor_file_ << "forces                      double       force"
      "               [numberOfParticles,3]" << endl;
  descriptor_file_ <<
      "#\n"
      "# END OF KIM DESCRIPTOR FILE\n"
      "#\n"
                   << endl;
}


int MultiLatticeKIM::get_neigh(void* kimmdl, int* mode, int* request, int* atom,
                               int* numnei, int** nei1atom, double** Rij)
{
  intptr_t* pkim = *((intptr_t**) kimmdl);
  int atomToReturn;
  int status;

  PPSumKIM* LatSum;
  LatSum = (PPSumKIM*) KIM_API_get_data(pkim, "neighObject", &status);

  int* numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles",
                                                   &status);

  if (0 == *mode) // iterator mode
  {
    if (0 == *request) // reset iterator
    {
      LatSum->Reset();
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    }
    else if (1 == *request) // increment iterator
    {
      if ((LatSum->CurrentPOS()) >= *numberOfParticles)
      {
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      }
      else
      {
        *numnei = LatSum->numNeigh();
        *atom = LatSum->CurrentPOS();
        *nei1atom = LatSum->nListAtom();
        *Rij = LatSum->nListRVec();
      }
      LatSum->operator++();
      return KIM_STATUS_OK;
    }
    else // invalid request value
    {
      KIM_API_report_error(__LINE__, (char*) __FILE__,
                           (char*) "Invalid request in get_periodic_neigh",
                           KIM_STATUS_NEIGH_INVALID_REQUEST);
      return KIM_STATUS_NEIGH_INVALID_REQUEST;
    }
  }
  else if (1 == *mode) // locator mode
  {
    if ((*request >= *numberOfParticles) || (*request < 0)) // invalid id
    {
      KIM_API_report_error(__LINE__, (char*) __FILE__,
                           (char*) "Invalid atom ID in get_neigh",
                           KIM_STATUS_PARTICLE_INVALID_ID);
      return KIM_STATUS_PARTICLE_INVALID_ID;
    }
    else
    {
      LatSum->Reset();
      for (int i = 0; i < *request; i++)
      {
        LatSum->operator++();
      }
      *numnei = LatSum->numNeigh();
      *atom = LatSum->CurrentPOS();
      *nei1atom = LatSum->nListAtom();
      *Rij = LatSum->nListRVec();
      return KIM_STATUS_OK;
    }
  }
  else // invalid mode
  {
    KIM_API_report_error(__LINE__, (char*) __FILE__,
                         (char*) "Invalid mode in get_periodic_neigh",
                         KIM_STATUS_NEIGH_INVALID_MODE);
    return KIM_STATUS_NEIGH_INVALID_MODE;
  }
}


int MultiLatticeKIM::process_dEdr(void* kimmdl, double* dEdr, double* r,
                                  double** dx, int* i, int* j)
{
  int status;
  intptr_t* pkim = *((intptr_t**) kimmdl);
  MultiLatticeKIM* obj;
  obj = (MultiLatticeKIM*) KIM_API_get_sim_buffer(pkim, &status);
  double DX[DIM3];

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
          += *dEdr * (obj->CBK_->DyDS(*dx, *i, *j, atom, k)) / (2.0 * (*r));
      obj->ME1_F_static[(obj->CBK_F_->INDS(atom, k))]
          += *dEdr * (obj->CBK_F_->DyDS(*dx, *i, *j, atom, k)) / (2.0 * (*r));
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
                              * obj->CBK_->D2yDSS(*i, *j, atom0, k1, atom1, l1)
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_->DyDS(*dx, *i, *j, atom0, k1))
                              * (obj->CBK_->DyDS(*dx, *i, *j, atom1, l1)));
            obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDS(atom1, l1)]
                += (*dEdr) * ((0.5 / (*r))
                              * obj->CBK_F_->D2yDSS(*i, *j, atom0, k1, atom1, l1)
                              - (0.25 / ((*r) * (*r) * (*r)))
                              * (obj->CBK_F_->DyDS(*dx, *i, *j, atom0, k1))
                              * (obj->CBK_F_->DyDS(*dx, *i, *j, atom1, l1)));
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
                   * obj->CBK_->D2yDFS(*dx, DX, *i, *j, i1, j1, atom0, k1)
                   - (0.25 / ((*r) * (*r) * (*r)))
                   * (obj->CBK_->DyDS(*dx, *i, *j, atom0, k1))
                   * (obj->CBK_->DyDF((*dx), DX, i1, j1)));

            obj->ME2_static[obj->CBK_->INDF(i1, j1)][obj->CBK_->INDS(atom0, k1)]
                += temp;

            obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDF(i1, j1)]
                += temp;

            double Ftemp=(*dEdr)
                * ((0.5 / (*r))
                   * obj->CBK_F_->D2yDFS(*dx, DX, *i, *j, i1, j1, atom0, k1)
                   - (0.25 / ((*r) * (*r) * (*r)))
                   * (obj->CBK_F_->DyDS(*dx, *i, *j, atom0, k1))
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

int MultiLatticeKIM::process_d2Edr2(void* kimmdl, double* d2Edr2, double** r,
                                    double** dx, int** i, int** j)
{
  int status;
  intptr_t* pkim = *((intptr_t**) kimmdl);
  MultiLatticeKIM* obj;
  obj = (MultiLatticeKIM*) KIM_API_get_sim_buffer(pkim, &status);

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
    if (((atom0 == (*i)[0]) || (atom0 == (*j)[0])))
    {
      for (int atom1 = 0; atom1 < (obj->CBK_->InternalAtoms()); atom1++)
      {
        if (((atom1 == (*i)[1]) || (atom1 == (*j)[1])))
        {
          for (k1 = 0; k1 < DIM3; k1++)
          {
            for (l1 = 0; l1 < DIM3; l1++)
            {
              obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDS(atom1, l1)]
                  += (0.5 / (*r)[0]) * (0.5 / (*r)[1]) * (*d2Edr2)
                  * (obj->CBK_->DyDS(Dx[0], (*i)[0], (*j)[0], atom0, k1))
                  * (obj->CBK_->DyDS(Dx[1], (*i)[1], (*j)[1], atom1, l1));
              obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDS(atom1, l1)]
                  += (0.5 / (*r)[0]) * (0.5 / (*r)[1]) * (*d2Edr2)
                  * (obj->CBK_F_->DyDS(Dx[0], (*i)[0], (*j)[0], atom0, k1))
                  * (obj->CBK_F_->DyDS(Dx[1], (*i)[1], (*j)[1], atom1, l1));
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
                 * obj->CBK_->DyDS(Dx[1], (*i)[1], (*j)[1], atom0, k1));

          obj->ME2_static[obj->CBK_->INDS(atom0, k1)][obj->CBK_->INDF(i1, j1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_->DyDF(Dx[1], DX[1], i1, j1)
                 * obj->CBK_->DyDS(Dx[0], (*i)[0], (*j)[0], atom0, k1));

          obj->ME2_F_static[obj->CBK_F_->INDF(i1, j1)][obj->CBK_F_->INDS(atom0, k1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_F_->DyDF(Dx[0], DX[0], i1, j1)
                 * obj->CBK_F_->DyDS(Dx[1], (*i)[1], (*j)[1], atom0, k1));

          obj->ME2_F_static[obj->CBK_F_->INDS(atom0, k1)][obj->CBK_F_->INDF(i1, j1)]
              += (0.5 / ((*r)[0])) * (0.5 / ((*r)[1]))
              * (*d2Edr2)
              * (obj->CBK_F_->DyDF(Dx[1], DX[1], i1, j1)
                 * obj->CBK_F_->DyDS(Dx[0], (*i)[0], (*j)[0], atom0, k1));
        }
      }
    }
  }

  return KIM_STATUS_OK;
}
