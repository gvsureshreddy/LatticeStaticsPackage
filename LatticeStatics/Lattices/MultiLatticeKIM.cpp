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
  delete CBK_;
  if(CBK_ != CBK_F_)
    delete CBK_F_;
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

  // set up CBK_kim object
  KIM_CBK_(CBK_, CBK_F_, Input);
  InfluenceDist_ = KIM_CBK_.get_cutoff();
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
  str_static.Resize(CBK_->DOFS());
  stiff_static.Resize(CBK_->DOFS(), CBK_->DOFS());
  CondEV_static.Resize(1, CBK_F_->Fsize());
  CondModuli_static.Resize(CBK_F_->Fsize(), CBK_F_->Fsize());
  OPM_static.Resize(CBK_F_->Ssize(), CBK_F_->Ssize());
  UnitCellShiftModuli_static.Resize(CBK_F_->Ssize(), CBK_F_->Ssize());
  OpticEV_static.Resize(1, CBK_F_->Ssize());
  OpticEV_Print.Resize(CBK_F_->Ssize());
  TestFunctVals_static.Resize(NumTestFunctions());
  ME1_static.Resize(CBK_->DOFS(), 0.0);
  ME1_F_static.Resize(CBK_F_->DOFS(), 0.0);
  ME2_static.Resize(CBK_->DOFS(), CBK_->DOFS(), 0.0);
  ME2_F_static.Resize(CBK_F_->DOFS(), CBK_F_->DOFS(), 0.0);
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
  KIM_CBK_.UpdateCoordinatesAndKIMValues();

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

  KIM_CBK_.UpdateCoordinatesAndKIMValues();

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
  KIM_CBK_.ComputeAndUpdate(StiffnessYes_);
  energy_ = KIM_CBK_.Energy();
  BodyForce_ = KIM_CBK_.get_BodyForce();
  ME1_static = KIM_CBK_.get_ME1_static();
  ME1_F_static = KIM_CBK_.get_ME1_F_static();
  if (StiffnessYes_ == 1)
  {
    ME2_static = KIM_CBK_.get_ME2_static();
    ME2_F_static = KIM_CBK_.get_ME2_F_static();
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
  KIM_CBK_.UpdateCoordinatesAndKIMValues();

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