#include "DisplacementControl.h"
#include "KnownLattices.h"
#include <fstream>

using namespace std;

DisplacementControl::~DisplacementControl()
{
  if (NULL != DispIndex_) delete [] DispIndex_;
  if (NULL != FreeIndex_) delete [] FreeIndex_;

  cout << "DisplacementControl Function Calls:\n"
       << "\tE0 calls - " << CallCount_[0] << "\n"
       << "\tE1 calls - " << CallCount_[1] << "\n"
       << "\tE1DLoad calls - " << CallCount_[2] << "\n"
       << "\tE2 calls - " << CallCount_[3] << "\n";
}

DisplacementControl::DisplacementControl(PerlInput const& Input,
                                         int const& Echo, int const& Width)
    :
    Lattice(Input, Echo),
    Lambda_(0.0),
    Width_(Width)
{
  LoadParameter_ = Load;
  for (int i = 0; i < cachesize; ++i)
  {
    Cached_[i] = 0;
    CallCount_[i] = 0;
  }

  // Get Lattice definition
  PerlInput::HashStruct Hash = Input.getHash("Lattice", "DisplacementControl");

  if (Input.ParameterOK(Hash, "WrappedLattice"))
  {
    char const* const LatticeType = Input.getString(Hash, "WrappedLattice");
    Lat_ = InitializeLattice(LatticeType, Input, Echo, Width);
  }
  else
  {
    cerr << "Error missing DisplacementControl{WrappedLattice} parameter\n";
    exit(-9);
  }

  LatDOFS_ = (Lat_->DOF()).Dim();

  if (Input.ParameterOK(Hash, "ProportionalDevice"))
  {
    DispDOFS_ = Input.getArrayLength(Hash,"ProportionalDevice");
    DispIndex_ = new int[DispDOFS_];
    Proportion_.Resize(DispDOFS_);
    Reference_.Resize(DispDOFS_);
    for (int i = 0; i < DispDOFS_; ++i)
    {
      if (3 != Input.getArrayLength(Hash, "ProportionalDevice", i))
      {
        cerr << "Error DisplacementControl{ProportionalDevice}[" << i
             << "] is not of length 3.\n";
        exit(-9);
      }
      DispIndex_[i] = Input.getPosInt(Hash, "ProportionalDevice", i, 0);
      Reference_[i] = Input.getDouble(Hash, "ProportionalDevice", i, 1);
      Proportion_[i] = Input.getDouble(Hash, "ProportionalDevice",i, 2);
    }

    DOFS_ = LatDOFS_ - DispDOFS_;
    FreeIndex_ = new int[DOFS_];
    int id = 0;
    int count = 0;
    while (count < DOFS_)
    {
      bool found = false;
      for (int i = 0; i < DispDOFS_; ++i)
      {
        if (id == DispIndex_[i]) found = true;
      }
      if (!found) FreeIndex_[count++] = id;
      ++id;
    }
  }
  else
  {
    cerr << "Error missing DisplacementControl{ProportionalDevice} parameter\n";
    exit(-9);
  }

  DOF_.Resize(DOFS_, 0.0);
  if (Input.ParameterOK(Hash, "InitialDOFs"))
  {
    Input.getVector(DOF_, Hash, "InitialDOFs");
  }

  E1CachedValue_.Resize(DOFS_);
  E1DLoadCachedValue_.Resize(DOFS_);
  E2CachedValue_.Resize(DOFS_, DOFS_);
  SetLatDOF();
}

void DisplacementControl::SetLatDOF()
{
  Vector LatDOF(LatDOFS_, 0.0);

  for (int i = 0; i < DispDOFS_; ++i)
  {
    LatDOF[DispIndex_[i]] = Reference_[i] + Lambda_*Proportion_[i];
  }
  for (int i = 0; i < DOFS_; ++i)
  {
    LatDOF[FreeIndex_[i]] = DOF_[i];
  }

  Lat_->SetDOF(LatDOF);
}


double DisplacementControl::E0() const
{
  if (!Cached_[0])
  {
    E0CachedValue_ = Lat_->E0();

    Cached_[0] = 1;
    CallCount_[0]++;
  }

  return E0CachedValue_;
}

Vector const& DisplacementControl::E1() const
{
  if (!Cached_[1])
  {
    Vector const& LatE1 = Lat_->E1();
    for (int i=0; i < DOFS_; ++i)
    {
      E1CachedValue_[i] = LatE1[FreeIndex_[i]];
    }

    Cached_[1] = 1;
    CallCount_[1]++;
  }

  return E1CachedValue_;
}

Vector const& DisplacementControl::E1DLoad() const
{
  if (!Cached_[2])
  {
    Matrix const& LatE2 = Lat_->E2();
    for (int i = 0; i < DOFS_; ++i)
    {
      E1DLoadCachedValue_[i] = 0.0;
      for (int j = 0; j < DispDOFS_; ++j)
      {
        E1DLoadCachedValue_[i]
            += LatE2[FreeIndex_[i]][DispIndex_[j]]*Proportion_[j];
      }
    }

    Cached_[2] = 1;
    CallCount_[2]++;
  }

  return E1DLoadCachedValue_;
}

Matrix const& DisplacementControl::E2() const
{
  if (!Cached_[3])
  {
    Matrix const& LatE2 = Lat_->E2();
    for (int i=0; i < DOFS_; ++i)
    {
      for (int j=0; j < DOFS_; ++j)
      {
        E2CachedValue_[i][j] = LatE2[FreeIndex_[i]][FreeIndex_[j]];
      }
    }

    Cached_[3] = 1;
    CallCount_[3]++;
  }

  return E2CachedValue_;
}

void DisplacementControl::Print(ostream& out, PrintDetail const& flag,
                                PrintPathSolutionType const& SolType)
{
  int W;
  W = out.width();

  double conjugateForce = 0.0;
  Vector LatE1 = Lat_->E1();
  for (int i=0; i<DispDOFS_; ++i)
  {
    conjugateForce += Proportion_[i]*LatE1[DispIndex_[i]];
  }

  out.width(0);
  if (Echo_)
  {
    cout.width(0);
  }

  switch (flag)
  {
    case PrintLong:
      out << "DisplacementControl:" << "\n";
      out << "  Controled DOFS:" << setw(3) << DispIndex_[0];
      for (int i=1; i < DispDOFS_; ++i)
      {
        out << "," << setw(3) << DispIndex_[i];
      }
      out << "\n";
      out << "  Free DOFS:     " << setw(3) << FreeIndex_[0];
      for (int i=1; i < DOFS_; ++i)
      {
        out << "," << setw(3) << FreeIndex_[i];
      }
      out << "\n\n";

      if (Echo_)
      {
        cout << "DisplacementControl:" << "\n";
        cout << "  Controled DOFS:" << setw(3) << DispIndex_[0];
        for (int i=1; i < DispDOFS_; ++i)
        {
          cout << "," << setw(3) << DispIndex_[i];
        }
        cout << "\n";
        cout << "  Free DOFS:     " << setw(3) << FreeIndex_[0];
        for (int i=1; i < DOFS_; ++i)
        {
          cout << "," << setw(3) << FreeIndex_[i];
        }
        cout << "\n\n";
      }
      // passthrough to short
    case PrintShort:
      out << "DispLambda :" << setw(W) << Lambda_ << "\n";
      out << "DispConjugaeToLambda :" << setw(W) << conjugateForce << "\n";
      out << "DispDOF's  :" << setw(W) << DOF_ << "\n";
      out << "DispE1     :" << setw(W) << E1() << "\n";
      out << "--\n";
      if (Echo_)
      {
        cout << "DispLambda :" << setw(W) << Lambda_ << "\n";
        cout << "DispConjugaeToLambda :" << setw(W) << conjugateForce << "\n";
        cout << "DispDOF's  :" << setw(W) << DOF_ << "\n";
        cout << "DispE1     :" << setw(W) << E1() << "\n";
        cout << "--\n";
      }
      out << setw(W);
      Lat_->Print(out, flag, SolType);
      break;
  }
}

ostream& operator<<(ostream& out, DisplacementControl& A)
{
  A.Print(out, Lattice::PrintShort);
  return out;
}
