#include "SCLDQMultiChainTPP.h"
#include "UtilityFunctions.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "time.h"
#include <Vector.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//  If you want to use OpenMP then make the following changes:
//
//  1) Change "CC = g++" to "CC = g++ -fopenmp" in Makefile
//  2) Recompile the entire code
//  3) Specify the number of processors by typing "export OMP_NUM_THREADS=8" in the terminal
//  4) Run the code: nice nohup bin/LatticeStatics input output >& term &
//  5) You can check the performance of CPU by typing "top"...the LatticeStatics code should take around 780 to 800% of CPU performance since we are using 8 processors


using namespace std;

const int SCLDQMultiChainTPP::DIM1 = 1;

SCLDQMultiChainTPP::~SCLDQMultiChainTPP()
{
   delete[] BodyForce_;
   delete[] SpeciesMass_;
   delete[] AtomicMass_;
   for (int i = 0; i < NumberofSpecies_; ++i)
   {
      for (int j = i; j < NumberofSpecies_; ++j)
      {
         delete SpeciesPotential_[i][j];
      }
   }
   delete[] SpeciesPotential_[0];
   delete[] SpeciesPotential_;
   delete[] Potential_[0];
   delete[] Potential_;
   delete[] AtomPositions_;

   delete[] HEigVals_static;

   for (int x = 0; x < DOFS; ++x)
   {
      for (int y = x; y < DOFS; ++y)
      {
         delete[] (HEigValsDOFDOF_static[x])[y];
      }
      delete[] HEigValsDOF_static[x];
      delete[] HEigValsDOFDOF_static[x];
   }
   delete[] HEigValsDOF_static;
   delete[] HEigValsDOFDOF_static;

   delete[] HEigVec_static;

   for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
   {
      for (int k = 0; k < INTERNAL_ATOMS; ++k)
      {
         for (int x = 0; x < DOFS; ++x)
         {
            delete[] ((HEigVecDOFDOF_static[i])[k])[x];
         }
         delete[] (HEigVecDOF_static[i])[k];
         delete[] (HEigVecDOFDOF_static[i])[k];
      }
      delete[] HEigVec_static[i];
      delete[] HEigVecDOF_static[i];
      delete[] HEigVecDOFDOF_static[i];
   }

   delete[] V4DOF_static;

   for (int x = 0; x < DOFS; ++x)
   {
      delete[] V4DOFDOF_static[x];
   }
   delete[] V4DOFDOF_static;

   delete[] EigVals_previous_static;

   delete[] Z_static;

   delete[] EigVals_static;
   delete[] EigValsT_static;
   delete[] EigValsTT_static;

   for (int x = 0; x < DOFS; ++x)
   {
      for (int y = x; y < DOFS; ++y)
      {
         delete[] (EigValsDOFDOF_static[x])[y];
      }
      delete[] EigValsDOF_static[x];
      delete[] EigValsTDOF_static[x];
      delete[] EigValsDOFDOF_static[x];
   }
   delete[] EigValsDOF_static;
   delete[] EigValsTDOF_static;
   delete[] EigValsDOFDOF_static;
}

SCLDQMultiChainTPP::SCLDQMultiChainTPP(PerlInput const& Input,
                                       int const& Echo, int const& Width,
                                       int const& Debug) :
   Lattice(Input, Echo)
{
   dbg_ = Debug;
   // Get Lattice definition
   stringstream tmp;

   PerlInput::HashStruct Hash = Input.getHash("Lattice", "SCLDQMultiChainTPP");
   INTERNAL_ATOMS = Input.getPosInt(Hash, "InternalAtoms");
   DOFS = 1 + INTERNAL_ATOMS;

   // Set RefLattice_
   RefLattice_.Resize(DIM1, DIM1);
   RefLattice_[0][0] = Input.getDouble(Hash, "LatticeBasis");

   // First Size DOF
   DOF_.Resize(DOFS, 0.0);
   DOF_[0] = 1.0;

   // Set AtomPositions_
   AtomPositions_ = new Vector[INTERNAL_ATOMS];
   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      AtomPositions_[i].Resize(DIM1);
      tmp.str("");
      tmp << "AtomPosition_" << i;
      AtomPositions_[i][0] = Input.getDouble(Hash, tmp.str().c_str());
   }

   // Setup Bodyforce_
   BodyForce_ = new Vector[INTERNAL_ATOMS];
   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      BodyForce_[i].Resize(DIM1, 0.0);
   }

   // Get Boltzmann constant
   kB_ = Input.getDouble(Hash, "kB");

   // Get Planck constant
   h_ = Input.getDouble(Hash, "h");

   // Get Thermo parameters
   Tref_ = Input.getDouble(Hash, "Tref");
   // PhiRef_ = Input.getDouble(Hash,"PhiRef");
   // EntropyRef_ = Input.getDouble(Hash,"EntropyRef");
   // HeatCapacityRef_ = Input.getDouble(Hash,"HeatCapacityRef");

   Input.getIntVector(AtomSpecies_, INTERNAL_ATOMS, Hash, "AtomSpecies");
   NumberofSpecies_ = AtomSpecies_[0];
   for (int i = 1; i < INTERNAL_ATOMS; ++i)
   {
      if (NumberofSpecies_ < AtomSpecies_[i])
      {
         NumberofSpecies_ = AtomSpecies_[i];
      }
   }
   NumberofSpecies_++;

   // Get Potential Parameters
   SpeciesPotential_ = new PairPotentials * *[NumberofSpecies_];
   SpeciesPotential_[0] = new PairPotentials *[NumberofSpecies_ * NumberofSpecies_];
   for (int i = 1; i < NumberofSpecies_; ++i)
   {
      SpeciesPotential_[i] = SpeciesPotential_[i - 1] + NumberofSpecies_;
   }
   Potential_ = new PairPotentials * *[INTERNAL_ATOMS];
   Potential_[0] = new PairPotentials *[INTERNAL_ATOMS * INTERNAL_ATOMS];
   for (int i = 1; i < INTERNAL_ATOMS; ++i)
   {
      Potential_[i] = Potential_[i - 1] + INTERNAL_ATOMS;
   }

   SpeciesMass_ = new double[NumberofSpecies_];
   AtomicMass_ = new double[INTERNAL_ATOMS];

   for (int i = 0; i < NumberofSpecies_; ++i)
   {
      for (int j = i; j < NumberofSpecies_; ++j)
      {
         SpeciesPotential_[i][j] =
            SpeciesPotential_[j][i]
               = InitializePairPotential(Hash, Input, i, j);
      }
      tmp.str("");
      tmp << "AtomicMass_" << i;
      SpeciesMass_[i] = Input.getDouble(Hash, tmp.str().c_str());
   }

   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      for (int j = i; j < INTERNAL_ATOMS; ++j)
      {
         Potential_[i][j] = Potential_[j][i]
                               = SpeciesPotential_[AtomSpecies_[i]][AtomSpecies_[j]];
      }

      AtomicMass_[i] = SpeciesMass_[AtomSpecies_[i]];
   }

   // Get Lattice parameters
   // NTemp_ = 1.0;
   RefNTemp_ = Input.getDouble(Hash, "RefNTemp");
   NTemp_ = RefNTemp_;
   InfluenceDist_ = Input.getDouble(Hash, "InfluenceDist");
   if (Input.ParameterOK(Hash, "Density"))
   {
      Density_ = Input.getInt(Hash, "Density");
   }
   else
   {
      Density_ = Input.useInt(1, Hash, "Density"); // Default Value
   }
   NormModulus_ = Input.getDouble(Hash, "NormModulus");

   // Get Loading parameters
   const char* loadparam = Input.getString(Hash, "LoadingParameter");
   if (!strcmp("Temperature", loadparam))
   {
      LoadParameter_ = Temperature;
   }
   else if (!strcmp("Load", loadparam))
   {
      LoadParameter_ = Load;
   }
   else
   {
      cerr << "Unknown Loading Parameter" << "\n"; exit(-1);
   }
   // Lambda_ = 0.0;
   RefLambda_ = Input.getDouble(Hash, "RefLambda");
   Lambda_ = RefLambda_;

   // needed to initialize reference length
   int iter;
   iter = Input.getPosInt(Hash, "MaxIterations");
   GridSize_ = Input.getPosInt(Hash, "BlochWaveGridSize");

   if (GridSize_ % 2 == 1)
   {
      cout << "GridSize is odd. Please keep it some even number. Even number is good for converged results when GridSize is increased." << endl;
      exit(-1);
   }

   // set LagrangeCB_
   const char* CBKin = Input.getString(Hash, "CBKinematics");
   if (!strcmp("LagrangeCB", CBKin))
   {
      LagrangeCB_ = 1;
   }
   else if (!strcmp("MixedCB", CBKin))
   {
      LagrangeCB_ = 0;
   }
   else
   {
      LagrangeCB_ = 1;
   }

   // Initiate the Lattice Sum object
   SCLDChainSum_(&DOF_, LagrangeCB_, 1, &RefLattice_, INTERNAL_ATOMS, AtomPositions_, Potential_,
                 &InfluenceDist_, &NTemp_);

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   ChainIter_(GridSize_, 1, 0);
   ChainIterNew_(GridSize_, 1, 0);


   HEigVals_static = new Matrix[GridSize_ / 2 + 1];
   HEigValsDOF_static = new Matrix *[DOFS];
   HEigValsDOFDOF_static = new Matrix * *[DOFS];

   for (int x = 0; x < DOFS; ++x)
   {
      HEigValsDOF_static[x] = new Matrix[GridSize_ / 2 + 1];
      HEigValsDOFDOF_static[x] = new Matrix *[DOFS];

      for (int y = x; y < DOFS; ++y)
      {
         (HEigValsDOFDOF_static[x])[y] = new Matrix[GridSize_ / 2 + 1];
      }
   }


   HEigVec_static = new CMatrix *[GridSize_ / 2 + 1];
   HEigVecDOF_static = new CMatrix * *[GridSize_ / 2 + 1];
   HEigVecDOFDOF_static = new CMatrix * **[GridSize_ / 2 + 1];

   for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
   {
      HEigVec_static[i] = new CMatrix[INTERNAL_ATOMS];
      HEigVecDOF_static[i] = new CMatrix *[INTERNAL_ATOMS];
      HEigVecDOFDOF_static[i] = new CMatrix * *[INTERNAL_ATOMS];

      for (int k = 0; k < INTERNAL_ATOMS; ++k)
      {
         (HEigVecDOF_static[i])[k] = new CMatrix[DOFS];
         (HEigVecDOFDOF_static[i])[k] = new CMatrix *[DOFS];

         for (int x = 0; x < DOFS; ++x)
         {
            ((HEigVecDOFDOF_static[i])[k])[x] = new CMatrix[DOFS];
         }
      }
   }


   V4DOF_static = new Matrix[DOFS];
   V4DOFDOF_static = new Matrix *[DOFS];

   for (int x = 0; x < DOFS; ++x)
   {
      V4DOFDOF_static[x] = new Matrix[DOFS];
   }

   EigVals_previous_static = new Matrix[GridSize_ / 2 + 1];

   Z_static = new Vector[GridSize_ / 2 + 1];

   EigVals_static = new Matrix[GridSize_ / 2 + 1];
   EigValsT_static = new Matrix[GridSize_ / 2 + 1];
   EigValsTT_static = new Matrix[GridSize_ / 2 + 1];

   EigValsDOF_static = new Matrix *[DOFS];
   EigValsTDOF_static = new Matrix *[DOFS];
   EigValsDOFDOF_static = new Matrix * *[DOFS];

   for (int x = 0; x < DOFS; ++x)
   {
      EigValsDOF_static[x] = new Matrix[GridSize_ / 2 + 1];
      EigValsTDOF_static[x] = new Matrix[GridSize_ / 2 + 1];
      EigValsDOFDOF_static[x] = new Matrix *[DOFS];

      for (int y = x; y < DOFS; ++y)
      {
         (EigValsDOFDOF_static[x])[y] = new Matrix[GridSize_ / 2 + 1];
      }
   }

   if (Input.ParameterOK(Hash, "InitialEigVals"))
   {
      InitialEigVals_available = 1;
      InitialEigVals_static.Resize((GridSize_ / 2 + 1), INTERNAL_ATOMS);
      Input.getMatrix(InitialEigVals_static, Hash, "InitialEigVals");

      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         EigVals_previous_static[i].Resize(1, INTERNAL_ATOMS, 0.0);

         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            (EigVals_previous_static[i])[0][k] = InitialEigVals_static[i][k];
         }
      }
   }
   else
   {
      InitialEigVals_available = 0;
   }

   int err = 0;
   err = FindLatticeSpacing(iter);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << "\n";
      exit(-1);
   }

   // Setup initial status for parameters

   // NTemp_ = Input.getDouble(Hash,"NTemp");
   // Lambda_ = Input.getDouble(Hash,"Lambda");
   NTemp_ = RefNTemp_;
   Lambda_ = RefLambda_;

   // Make any changes to atomic potentials that might be required
   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      for (int j = i; j < INTERNAL_ATOMS; ++j)
      {
         if (AtomSpecies_[i] < AtomSpecies_[j])
         {
            UpdatePairPotential(Hash, Input,
                                AtomSpecies_[i], AtomSpecies_[j], Potential_[i][j]);
         }
         else
         {
            UpdatePairPotential(Hash, Input,
                                AtomSpecies_[j], AtomSpecies_[i], Potential_[j][i]);
         }
      }
   }
   SCLDChainSum_.Recalc();

   Input.EndofInputSection();
}

int SCLDQMultiChainTPP::FindLatticeSpacing(int const& iter)
{
   // Lambda_=0.0;
   // NTemp_=1.0;
   Lambda_ = RefLambda_;
   NTemp_ = RefNTemp_;
   DOF_[0] = 1.0;
   for (int i = 1; i < DOFS; i++)
   {
      DOF_[i] = 0.0;
   }

   SCLDChainSum_.Recalc();

   FreqCached = 0;
   FirstConverged = 0;
   FirstPrintLong = 0;

   if (Echo_)
   {
      RefineEqbm(1.0e-13, iter, &cout);
   }
   else
   {
      RefineEqbm(1.0e-13, iter, 0);
   }

   FreqCached = 0;

   // Clean up numerical round off (at least for zero values)
   for (int i = 0; i < DOFS; ++i)
   {
      if (fabs(DOF_[i]) < 1.0e-13)
      {
         DOF_[i] = 0.0;
      }
   }

   // Update RefLattice_
   RefLattice_ = RefLattice_ * DOF_[0];

   // update atom pos
   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      if (LagrangeCB_)
      {
         AtomPositions_[i][0] += DOF_[i + 1];
      }
      else
      {
         AtomPositions_[i][0] += DOF_[i + 1] / DOF_[0];
      }
   }

   // reset DOF
   DOF_[0] = 1.0;
   for (int i = 1; i < DOFS; i++)
   {
      DOF_[i] = 0.0;
   }

   SCLDChainSum_.Recalc();
   return 0;
}

void SCLDQMultiChainTPP::SetParameters(double const* const Vals, int const& ResetRef)
{
   int no = SpeciesPotential_[0][0]->GetNoParameters();
   int cur = 0;
   for (int i = 0; i < NumberofSpecies_; ++i)
   {
      for (int j = i; j < NumberofSpecies_; ++j)
      {
         SpeciesPotential_[i][j]->SetParameters(&(Vals[cur]));
         cur += no;
      }
   }

   SCLDChainSum_.Recalc();

   if (ResetRef)
   {
      FindLatticeSpacing(50);
   }
}

// Lattice Routines

double SCLDQMultiChainTPP::PI(double const* const Dx, double const* const DX) const
{
   return 2.0 * Dx[0] * DX[0];
}

double SCLDQMultiChainTPP::PSI(double const* const DX) const
{
   return 2.0 * DX[0] * DX[0];
}

double SCLDQMultiChainTPP::OMEGA(double const* const Dx, int const& p, int const& q,
                                 int const& i)
const
{
   return LagrangeCB_ ?
          2.0 * DOF_[0] * RefLattice_[0][0] * DELTA(i, p, q) * Dx[0] :
          2.0 * DELTA(i, p, q) * Dx[0];
}

double SCLDQMultiChainTPP::SIGMA(int const& p, int const& q, int const& i, int const& j)
const
{
   return LagrangeCB_ ?
          2.0 * DOF_[0] * DOF_[0] * RefLattice_[0][0] * DELTA(i, p, q) * RefLattice_[0][0] * DELTA(j, p, q) :
          2.0 * DELTA(i, p, q) * DELTA(j, p, q);
}

double SCLDQMultiChainTPP::GAMMA(double const* const Dx, double const* const DX,
                                 int const& p, int const& q, int const& i) const
{
   return LagrangeCB_ ?
          4.0 * RefLattice_[0][0] * DELTA(i, p, q) * Dx[0] :
          2.0 * DELTA(i, p, q) * DX[0];
}

double SCLDQMultiChainTPP::THETA(double const* const DX, int const& p, int const& q,
                                 int const& i) const
{
   return LagrangeCB_ ?
          4.0 * RefLattice_[0][0] * DELTA(i, p, q) * DX[0] :
          0.0;
}

double SCLDQMultiChainTPP::XI(int const& p, int const& q, int const& i, int const& j)
const
{
   return LagrangeCB_ ?
          4.0 * DOF_[0] * RefLattice_[0][0] * DELTA(i, p, q) * RefLattice_[0][0] * DELTA(j, p, q) :
          0.0;
}

double SCLDQMultiChainTPP::LAMDA(int const& p, int const& q, int const& i, int const& j)
const
{
   return LagrangeCB_ ?
          4.0 * RefLattice_[0][0] * DELTA(i, p, q) * RefLattice_[0][0] * DELTA(j, p, q) :
          0.0;
}

double SCLDQMultiChainTPP::E0() const
{
   double E0, Tsq;

   E0 = FreeEnergy();

   Tsq = 0;

   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      Tsq += DOF_[i + 1];
   }
   Tsq = (Tsq * Tsq) / INTERNAL_ATOMS;

   return E0 + 0.5 * Tsq;
}

double SCLDQMultiChainTPP::energy(PairPotentials::TDeriv const& dt) const
{
   double Phi = 0.0;
   double Vr;

   for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
   {
      // Calculate Phi
      Phi += Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
         NTemp_, SCLDChainSum_.r2(), PairPotentials::Y0, dt);
   }

   // Phi = Phi/(2*Vr*NormModulus)
   Vr = Density_ ? RefLattice_.Det() : 1.0;
   Phi *= 1.0 / (2.0 * (Vr * NormModulus_));

   // Apply loading potential and Thermal term
   if (dt == PairPotentials::T0)
   {
      // Loading
      Phi -= Lambda_ * (DOF_[0] - 1.0);

      // Thermal term
   }
   else if (dt == PairPotentials::DT)
   {
      // Loading

      // Thermal term
      // Phi += (-EntropyRef_ - HeatCapacityRef_*log(NTemp_*Tref_))/NormModulus_;
   }
   else if (dt == PairPotentials::D2T)
   {
      // Loading

      // Thermal term
      // Phi += (-HeatCapacityRef_/(NTemp_*Tref_))/NormModulus_;
   }
   else
   {
      cerr << "Error in SCLDQMultiChainTPP::energy" << "\n";
      exit(-1);
   }

   return Phi;
}

double SCLDQMultiChainTPP::FreeEnergy(PairPotentials::TDeriv const& dt) const
{
   double FPhi = 0.0;
   double Vr;
   double Term, TermT, TermTT;
   double Tol = 1.0e-13;


   ReferenceHarmonic();
   ReferenceV4();
   ReferencePseudoHarmonic();

   // Adding vibrational energy
   int i = 0;
   for (ChainIter_.Reset(); !ChainIter_.Done(); ++ChainIter_)
   {
      for (int k = 0; k < INTERNAL_ATOMS; ++k)
      {
         if ((EigVals_static[i])[0][k] > Tol)
         {
            if (dt == PairPotentials::T0)
            {
               Term = sqrt((EigVals_static[i])[0][k]);

               FPhi += 2.0 * (0.5 * h_ * Term + kB_ * NTemp_ * log(1.0 - exp(-h_ * Term / (kB_ * NTemp_))));
            }
            else if (dt == PairPotentials::DT)
            {
               Term = sqrt((EigVals_static[i])[0][k]);
               TermT = (EigValsT_static[i])[0][k] / (2 * Term);

               FPhi += 2.0 * (
                  0.5 * h_ * TermT + kB_ * log(1.0 - exp(-h_ * Term / (kB_ * NTemp_)))
                  + h_ * (TermT - Term / NTemp_) / (exp(h_ * Term / (kB_ * NTemp_)) - 1.0)
                             );
            }
            else if (dt == PairPotentials::D2T)
            {
               Term = sqrt((EigVals_static[i])[0][k]);
               TermT = (EigValsT_static[i])[0][k] / (2 * Term);
               TermTT = ((EigValsTT_static[i])[0][k] - 2 * TermT * TermT) / (2 * Term);

               FPhi += 2.0 * (
                  0.5 * h_ * TermTT
                  + h_ * TermTT / (exp(h_ * Term / (kB_ * NTemp_)) - 1.0)
                  - h_ * (h_ / (kB_ * NTemp_)) * (TermT - Term / NTemp_) * (TermT - Term / NTemp_) * exp(h_ * Term / (kB_ * NTemp_))
                  / ((exp(h_ * Term / (kB_ * NTemp_)) - 1.0) * (exp(h_ * Term / (kB_ * NTemp_)) - 1.0))
                             );
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::FreeEnergy" << "\n";
               exit(-1);
            }
         }
      }
      i = i + 1;
   }

   // (GridSize_) number of unit cells when GridSize_ is odd in 1D
   Vr = Density_ ? (GridSize_ + 2) * RefLattice_.Det() : 1.0;
   FPhi *= 1.0 / (Vr * NormModulus_); // Note that there won't be 1/2 here


   // Adding potential energy
   if (dt == PairPotentials::T0)
   {
      FPhi += energy(PairPotentials::T0);
   }
   else if (dt == PairPotentials::DT)
   {
      FPhi += energy(PairPotentials::DT);
   }
   else if (dt == PairPotentials::D2T)
   {
      FPhi += energy(PairPotentials::D2T);
   }
   else
   {
      cerr << "Error in SCLDQMultiChainTPP::FreeEnergy" << "\n";
      exit(-1);
   }

   return FPhi;
}

Vector const& SCLDQMultiChainTPP::E1() const
{
   Phi1_static.Resize(DOFS, 0.0);
   double T = 0.0;

   Phi1_static = Fstress();

   for (int i = 0; i < INTERNAL_ATOMS; ++i)
   {
      T += DOF_[i + 1];
   }
   T = T / INTERNAL_ATOMS;
   for (int i = 1; i < DOFS; ++i)
   {
      Phi1_static[i] += T;
   }

   return Phi1_static;
}

Vector const& SCLDQMultiChainTPP::stress(PairPotentials::TDeriv const& dt,
                                         LDeriv const& dl) const
{
   double ForceNorm = 0.0;
   double phi, Vr;
   double Tol = 1.0e-13;
   int i;

   stress_static.Resize(DOFS, 0.0);

   Vr = Density_ ? RefLattice_.Det() : 1.0;

   if (dl == L0)
   {
      for (i = 0; i < INTERNAL_ATOMS; ++i)
      {
         BodyForce_[i][0] = 0.0;
      }

      for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
      {
         // Calculate bodyforce
         // NOTE: phi1 = d(phi)/d(r2)
         // We need d(phi)/dr = 2*r*d(phi)/d(r2)
         phi = 2.0 * sqrt(SCLDChainSum_.r2()) * SCLDChainSum_.phi1();
         if (ForceNorm < fabs(-phi / 2.0))
         {
            ForceNorm = fabs(-phi / 2.0);
         }
         BodyForce_[SCLDChainSum_.Atom(0)][0] += -phi* SCLDChainSum_.Dx(0) / (2.0 * sqrt(SCLDChainSum_.r2()));


         // Claculate the Stress
         if (dt == PairPotentials::T0)
         {
            phi = SCLDChainSum_.phi1();
         }
         else if (dt == PairPotentials::DT)
         {
            phi = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
               NTemp_, SCLDChainSum_.r2(), PairPotentials::DY, dt);
         }
         else
         {
            cerr << "Error in SCLDQMultiChainTPP::stress" << "\n";
            exit(-1);
         }

         stress_static[0] += phi * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX());
         for (i = 0; i < INTERNAL_ATOMS; i++)
         {
            stress_static[i + 1]
               += phi * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i);

            //               cout << "DOF:\t" << i+1 << "OMEGA:\t";
            //               cout << OMEGA(SCLDChainSum_.pDx(),SCLDChainSum_.Atom(0),SCLDChainSum_.Atom(1),i) <<endl;
         }
      }

      // BodyForce[i] = BodyForce[i] / ForceNorm
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         BodyForce_[i][0] /= ForceNorm;
      }

      // stress_static = stress_static/(2*Vr*NormModulus)
      stress_static *= 1.0 / (2.0 * (Vr * NormModulus_));

      // Load terms
      if (dt == PairPotentials::T0)
      {
         stress_static[0] -= Lambda_;
      }
   }
   else if (dl == DL)
   {
      // dl=DL
      stress_static[0] -= 1.0;
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiChainTpp::stress()" << "\n";
      exit(-1);
   }

   // Clean up numerical round off (at least for zero values)
   for (int k = 0; k < DOFS; ++k)
   {
      if (fabs(stress_static[k]) < Tol)
      {
         stress_static[k] = 0.0;
      }
   }

   //   cout << "stress:" << endl;
   //   cout << setw(20) << stress_static << endl;

   return stress_static;
}

Vector const& SCLDQMultiChainTPP::Fstress(PairPotentials::TDeriv const& dt, LDeriv const& dl) const
{
   double Vr;
   double Term, TermT, TermDOF, TermTDOF;
   double Tol = 1.0e-13;

   Fstress_static.Resize(DOFS, 0.0);

   ReferenceHarmonic();
   ReferenceV4();
   ReferencePseudoHarmonic();

   if (dl == L0)
   {
      int i = 0;
      for (ChainIter_.Reset(); !ChainIter_.Done(); ++ChainIter_)
      {
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            if ((EigVals_static[i])[0][k] > Tol)
            {
               // Claculate the Stress due to vibrational energy
               if (dt == PairPotentials::T0)
               {
                  Term = sqrt((EigVals_static[i])[0][k]);

                  for (int x = 0; x < DOFS; ++x)
                  {
                     TermDOF = ((EigValsDOF_static[x])[i])[0][k] / (2.0 * Term);

                     Fstress_static[x] += 2.0 * (h_ * TermDOF * (0.5 + 1.0 / (exp(h_ * Term / (kB_ * NTemp_)) - 1.0)));
                  }
               }
               else if (dt == PairPotentials::DT)
               {
                  Term = sqrt((EigVals_static[i])[0][k]);
                  TermT = (EigValsT_static[i])[0][k] / (2.0 * Term);

                  for (int x = 0; x < DOFS; ++x)
                  {
                     TermDOF = ((EigValsDOF_static[x])[i])[0][k] / (2.0 * Term);
                     TermTDOF = (((EigValsTDOF_static[x])[i])[0][k] - 2.0 * TermDOF * TermT) / (2.0 * Term);

                     Fstress_static[x] += 2.0 * (
                        h_ * TermTDOF * (0.5 + 1.0 / (exp(h_ * Term / (kB_ * NTemp_)) - 1.0))
                        - h_ * (h_ / (kB_ * NTemp_)) * (TermT - Term / NTemp_) * TermDOF * exp(h_ * Term / (kB_ * NTemp_))
                        / ((exp(h_ * Term / (kB_ * NTemp_)) - 1.0) * (exp(h_ * Term / (kB_ * NTemp_)) - 1.0))
                                                );
                  }
               }
               else
               {
                  cerr << "Error in SCLDQMultiChainTPP::Fstress, D2T" << "\n";
                  exit(-1);
               }
            }
         }
         i = i + 1;
      }

      // (GridSize_) number of unit cells when GridSize_ is odd in 1D
      Vr = Density_ ? (GridSize_ + 2) * RefLattice_.Det() : 1.0;

      // Fstress_static = Fstress_static/(Vr*NormModulus)
      Fstress_static *= 1.0 / (Vr * NormModulus_);
   }
   else if (dl == DL)
   {
      // dl=DL
      //      Nothing to do ?? Ask elliott;
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiChainTpp::Fstress()" << "\n";
      exit(-1);
   }

   // Clean up numerical round off (at least for zero values)
   for (int k = 0; k < DOFS; ++k)
   {
      if (fabs(Fstress_static[k]) < Tol)
      {
         Fstress_static[k] = 0.0;
      }
   }

   //    cout << "Vib. stress:" << endl;
   //    cout << setw(20) << Fstress_static << endl;

   Fstress_static += stress(dt, dl);

   return Fstress_static;
}

Matrix const& SCLDQMultiChainTPP::E2() const
{
   Phi2_static.Resize(DOFS, DOFS, 0.0);
   static int i, j;

   Phi2_static = Fstiffness();

   // Add translational stiffness
   double factor = 1.0 / INTERNAL_ATOMS;
   for (i = 0; i < INTERNAL_ATOMS; ++i)
   {
      for (j = 0; j < INTERNAL_ATOMS; ++j)
      {
         Phi2_static[i + 1][j + 1] += factor;
      }
   }

   return Phi2_static;
}

Matrix const& SCLDQMultiChainTPP::stiffness(PairPotentials::TDeriv const& dt,
                                            LDeriv const& dl) const
{
   double phi, phi1;
   double Vr;
   double Tol = 1.0e-13;
   int i, j;

   stiff_static.Resize(DOFS, DOFS, 0.0);

   Vr = Density_ ? RefLattice_.Det() : 1.0;

   if (dl == L0)
   {
      for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
      {
         if (dt == PairPotentials::T0)
         {
            phi = SCLDChainSum_.phi2();
            phi1 = SCLDChainSum_.phi1();
         }
         else if (dt == PairPotentials::DT)
         {
            phi = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
               NTemp_, SCLDChainSum_.r2(), PairPotentials::D2Y, dt);
            phi1 = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
               NTemp_, SCLDChainSum_.r2(), PairPotentials::DY, dt);
         }
         else
         {
            cerr << "Error in SCLDQMultiChainTPP::stiffness" << "\n";
            exit(-1);
         }

         // Upper Diag Block (1,1)
         stiff_static[0][0] += phi * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                      * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()))
                               + phi1* PSI(SCLDChainSum_.pDX());

         // Lower Diag Block (INTERNAL_ATOMS,INTERNAL_ATOMS)
         for (i = 0; i < INTERNAL_ATOMS; i++)
         {
            for (j = 0; j < INTERNAL_ATOMS; j++)
            {
               stiff_static[i + 1][j + 1] +=
                  phi * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                         * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j))
                  + phi1* SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j);
            }
         }

         // Off Diag Blocks
         for (i = 0; i < INTERNAL_ATOMS; i++)
         {
            stiff_static[0][i + 1] = stiff_static[i + 1][0] +=
                                        phi * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                               * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i))
                                        + phi1* GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i);
         }
      }
      stiff_static *= 1.0 / (2.0 * (Vr * NormModulus_));
   }
   else if (dl == DL)
   {
      // Nothing to do: Phi2_static is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiChainTpp::stiffness()" << "\n";
      exit(-1);
   }


   // Clean up numerical round off (at least for zero values)
   for (int k = 0; k < DOFS; ++k)
   {
      for (int l = 0; l < DOFS; ++l)
      {
         if (fabs(stiff_static[k][l]) < Tol)
         {
            stiff_static[k][l] = 0.0;
         }
      }
   }


   //   cout << "stiffness:" << endl;
   //   cout << setw(20) << stiff_static << endl;

   return stiff_static;
}

Matrix const& SCLDQMultiChainTPP::Fstiffness(PairPotentials::TDeriv const& dt, LDeriv const& dl) const
{
   double Vr;
   double Term, TermDOF, TermDOF2, TermDOFDOF;
   double Tol = 1.0e-13;

   FPhi2_static.Resize(DOFS, DOFS, 0.0);

   ReferenceHarmonic();
   ReferenceV4();
   ReferencePseudoHarmonic();

   if (dl == L0)
   {
      int i = 0;
      for (ChainIter_.Reset(); !ChainIter_.Done(); ++ChainIter_)
      {
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            if ((EigVals_static[i])[0][k] > Tol)
            {
               if (dt == PairPotentials::T0)
               {
                  Term = sqrt((EigVals_static[i])[0][k]);

                  for (int x = 0; x < DOFS; ++x)
                  {
                     TermDOF = ((EigValsDOF_static[x])[i])[0][k] / (2.0 * Term);

                     for (int y = x; y < DOFS; ++y)
                     {
                        TermDOF2 = ((EigValsDOF_static[y])[i])[0][k] / (2.0 * Term);
                        TermDOFDOF = ((((EigValsDOFDOF_static[x])[y])[i])[0][k] - 2.0 * TermDOF * TermDOF2) / (2.0 * Term);

                        FPhi2_static[x][y] =
                           FPhi2_static[y][x] += 2.0 * (
                              h_ * TermDOFDOF * (0.5 + 1.0 / (exp(h_ * Term / (kB_ * NTemp_)) - 1.0))
                              - h_ * (h_ / (kB_ * NTemp_)) * TermDOF * TermDOF2 * exp(h_ * Term / (kB_ * NTemp_))
                              / ((exp(h_ * Term / (kB_ * NTemp_)) - 1.0) * (exp(h_ * Term / (kB_ * NTemp_)) - 1.0))
                                                       );
                     }
                  }
               }
               else
               {
                  cerr << "Error in SCLDQMultiChainTPP::Fstiffness, DT" << "\n";
                  exit(-1);
               }
            }
         }
         i = i + 1;
      }

      // (GridSize_) number of unit cells when GridSize_ is odd in 1D
      Vr = Density_ ? (GridSize_ + 2) * RefLattice_.Det() : 1.0;

      // Fstiffness_static = Fstiffness_static/(Vr*NormModulus)
      FPhi2_static *= 1.0 / (Vr * NormModulus_);
   }
   else if (dl == DL)
   {
      // Nothing to do: Phi2_static is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in SCLDQMultiChainTPP::Fstiffness()" << "\n";
      exit(-1);
   }

   // Clean up numerical round off (at least for zero values)
   for (int k = 0; k < DOFS; ++k)
   {
      for (int l = 0; l < DOFS; ++l)
      {
         if (fabs(FPhi2_static[k][l]) < Tol)
         {
            FPhi2_static[k][l] = 0.0;
         }
      }
   }


   //    cout << "Vib. stiffness:" << endl;
   //    cout << setw(20) << FPhi2_static << endl;

   FPhi2_static += stiffness(dt, dl);

   return FPhi2_static;
}

Matrix const& SCLDQMultiChainTPP::E3() const
{
   double phi, phi1, phi2;
   int i, j, k;

   Phi3_static.Resize(DOFS * DOFS, DOFS, 0.0);

   for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
   {
      phi = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
         NTemp_, SCLDChainSum_.r2(), PairPotentials::D3Y, PairPotentials::T0);
      phi1 = SCLDChainSum_.phi2();
      phi2 = SCLDChainSum_.phi1();

      // DF^3 block
      Phi3_static[0][0] +=
         phi * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()))
         + phi1 * (3.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                   * PSI(SCLDChainSum_.pDX()));
      // DV^3 block
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         for (j = 0; j < INTERNAL_ATOMS; j++)
         {
            for (k = 0; k < INTERNAL_ATOMS; k++)
            {
               Phi3_static[(i + 1) * DOFS + (j + 1)][k + 1] +=
                  phi * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                         * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                         * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k))
                  + phi1 * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                            + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                            * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, k)
                            + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                            * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j));
            }
         }
      }
      // DU^2DV blocks
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         Phi3_static[0][i + 1] =
            Phi3_static[(i + 1) * DOFS][0] =
               Phi3_static[i + 1][0] += (
                  phi * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                         * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                         * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i))
                  + phi1 * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                            * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            + PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                            * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            * PSI(SCLDChainSum_.pDX()))
                  + phi2 * THETA(SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i));
      }
      // DV^2DU blocks
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         for (j = 0; j < INTERNAL_ATOMS; j++)
         {
            Phi3_static[(i + 1) * DOFS + (j + 1)][0] =
               Phi3_static[(i + 1) * DOFS][j + 1] =
                  Phi3_static[i + 1][j + 1] += (
                     phi * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                            * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()))
                     + phi1 * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                               * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                       SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                               + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                               * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                       SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                               + PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                               * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j))
                     + phi2 * XI(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j));
         }
      }
   }


   // Phi3_static = Phi3_static/(2*Vr*NormModulus)
   Phi3_static *= 1.0 / (2.0 * ((Density_ ? RefLattice_.Det() : 1.0) * NormModulus_));

   return Phi3_static;
}

Matrix const& SCLDQMultiChainTPP::E4() const
{
   double phi, phi1, phi2, phi3;
   int i, j, k, m;

   Phi4_static.Resize(DOFS * DOFS, DOFS * DOFS, 0.0);

   for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
   {
      phi = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
         NTemp_, SCLDChainSum_.r2(), PairPotentials::D4Y, PairPotentials::T0);
      phi1 = Potential_[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)]->PairPotential(
         NTemp_, SCLDChainSum_.r2(), PairPotentials::D3Y, PairPotentials::T0);
      phi2 = SCLDChainSum_.phi2();
      phi3 = SCLDChainSum_.phi1();

      // DU^4 block
      Phi4_static[0][0] +=
         phi * (PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()))
         + phi1 * (
            6.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
            * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
            * PSI(SCLDChainSum_.pDX()))
         + phi2 * (
            3.0 * PSI(SCLDChainSum_.pDX())
            * PSI(SCLDChainSum_.pDX()));
      // DV^4 block
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         for (j = 0; j < INTERNAL_ATOMS; j++)
         {
            for (k = 0; k < INTERNAL_ATOMS; k++)
            {
               for (m = 0; m < INTERNAL_ATOMS; m++)
               {
                  Phi4_static[(i + 1) * DOFS + (j + 1)][(k + 1) * DOFS + (m + 1)] +=
                     phi * (OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), m))
                     + phi1 * (
                        OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k, m)
                        + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, m)
                        + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), m)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                        + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, m)
                        + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), m)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, k)
                        + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), m)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j))
                     + phi2 * (
                        SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k, m)
                        + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, m)
                        + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, m)
                        * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k));
               }
            }
         }
      }

      // DU^3DV blocks
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         Phi4_static[0][i + 1] =
            Phi4_static[0][(i + 1) * DOFS] =
               Phi4_static[i + 1][0] =
                  Phi4_static[(i + 1) * DOFS][0] += (
                     phi * (
                        PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i))
                     + phi1 * (
                        3.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        + 3.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        * PSI(SCLDChainSum_.pDX()))
                     + phi2 * (
                        3.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                        * THETA(SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        + 3.0 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                        * PSI(SCLDChainSum_.pDX())));
      }
      // DV^3DU blocks
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         for (j = 0; j < INTERNAL_ATOMS; j++)
         {
            for (k = 0; k < INTERNAL_ATOMS; k++)
            {
               Phi4_static[(i + 1) * DOFS + (j + 1)][(k + 1) * DOFS] =
                  Phi4_static[(i + 1) * DOFS + (j + 1)][k + 1] =
                     Phi4_static[(i + 1) * DOFS][(j + 1) * DOFS + (k + 1)] =
                        Phi4_static[i + 1][(j + 1) * DOFS + (k + 1)] += (
                           phi * (
                              OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()))
                           + phi1 * (
                              OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                              * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                              * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, k)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                              * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j))
                           + phi2 * (
                              OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                              * XI(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              * XI(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, k)
                              + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              * XI(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                              + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), k)
                              + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, k)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                              + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j, k)
                              * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                      SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)));
            }
         }
      }
      // DU^2DV^2 blocks
      for (i = 0; i < INTERNAL_ATOMS; i++)
      {
         for (j = 0; j < INTERNAL_ATOMS; j++)
         {
            Phi4_static[0][(i + 1) * DOFS + (j + 1)] =
               Phi4_static[(j + 1) * DOFS][i + 1] =
                  Phi4_static[(i + 1) * DOFS + (j + 1)][0] =
                     Phi4_static[i + 1][(j + 1) * DOFS] += (
                        phi * (
                           PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j))
                        + phi1 * (
                           PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                           + 2.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                   SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           + 2.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                   SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           * PSI(SCLDChainSum_.pDX()))
                        + phi2 * (
                           2.0 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                           * XI(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                           + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           * THETA(SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           + OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           * THETA(SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           + 2.0 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                         SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i)
                           * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(),
                                   SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), j)
                           + SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j)
                           * PSI(SCLDChainSum_.pDX()))
                        + phi3 * LAMDA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), i, j));
         }
      }
   }


   // Phi4_static = Phi4_static/(2*Vr*NormModulus)
   Phi4_static *= 1.0 / (2.0 * ((Density_ ? RefLattice_.Det() : 1.0) * NormModulus_));

   return Phi4_static;
}

Matrix const& SCLDQMultiChainTPP::CondensedModuli() const
{
   Matrix stiff = E2();
   int intrn = DOFS - 1;
   Matrix IM(intrn, intrn);
   CM_static.Resize(1, 1);

   CM_static[0][0] = stiff[0][0];

   // Make sure there are internal DOF's
   if (intrn)
   {
      for (int i = 0; i < intrn; i++)
      {
         for (int j = 0; j < intrn; j++)
         {
            IM[i][j] = stiff[1 + i][1 + j];
         }
      }
      IM = IM.Inverse();

      // Set up Condensed Moduli
      for (int m = 0; m < intrn; m++)
      {
         for (int n = 0; n < intrn; n++)
         {
            CM_static[0][0] -= stiff[0][1 + m] * IM[m][n] * stiff[1 + n][0];
         }
      }
   }

   return CM_static;
}

Vector const& SCLDQMultiChainTPP::ThermalExpansion() const
{
   ThermalExp_static.Resize(DOFS);
#ifdef SOLVE_SVD
   return ThermalExp_static = SolveSVD(E2(), -StressDT());
#else
   return ThermalExp_static = SolvePLU(E2(), -StressDT());
#endif
}

int SCLDQMultiChainTPP::comp(void const* const a, void const* const b)
{
   double t;
   if (*((double*) a) == *((double*) b))
   {
      return 0;
   }
   else
   {
      t = *((double*) a) - *((double*) b);
      t /= fabs(t);
      return int(t);
   }
}

int SCLDQMultiChainTPP::abscomp(void const* const a, void const* const b)
{
   double t;
   if (fabs(*((double*) a)) == fabs(*((double*) b)))
   {
      return 0;
   }
   else
   {
      t = fabs(*((double*) a)) - fabs(*((double*) b));
      t /= fabs(t);
      return int(t);
   }
}

void SCLDQMultiChainTPP::interpolate(Matrix* const EigVals, int const& zero,
                                     int const& one, int const& two)
{
   // Calculate expected value for eigvals and store in zero position
   EigVals[zero] = 2.0 * EigVals[one] - EigVals[zero];

   double delta, dtmp;
   int i, j, pos;

   for (i = 0; i < EigVals[0].Cols(); ++i)
   {
      pos = i;
      delta = fabs(EigVals[zero][0][i] - EigVals[two][0][i]);
      for (j = i + 1; j < EigVals[0].Cols(); ++j)
      {
         dtmp = fabs(EigVals[zero][0][i] - EigVals[two][0][j]);
         if (dtmp < delta)
         {
            delta = dtmp;
            pos = j;
         }
      }
      // move correct eigval to current pos
      dtmp = EigVals[two][0][i];
      EigVals[two][0][i] = EigVals[two][0][pos];
      EigVals[two][0][pos] = dtmp;
   }
}

CMatrix const& SCLDQMultiChainTPP::ReferenceDynamicalStiffness(Vector const& K, PairPotentials::TDeriv const& dt, int const& DOFderiv, int const& x, int const& y) const
{
   double pi = 4.0 * atan(1.0);
   double InverseLat;
   InverseLat = 1.0 / RefLattice_[0][0];
   MyComplexDouble Ic(0, 1);
   MyComplexDouble A = (2.0 * pi - 2.0 * pi / 2.5) * InverseLat;
   MyComplexDouble B = (2.0 * pi / 4.0) * InverseLat;
   MyComplexDouble Z = (A * K[0] + B) * Ic;


   double Term1;
   double Term2;

   Dk_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);

   for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
   {
      // Calculate Dk_static
      if (SCLDChainSum_.Atom(0) != SCLDChainSum_.Atom(1))
      {
         // y != y' terms (i.e., off diagonal terms)
         if (DOFderiv == 0)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2());
            }
            else if (dt == PairPotentials::DT)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1T()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2T());
            }
            else if (dt == PairPotentials::D2T)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1TT()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2TT());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 0,D3T" << "\n";
               exit(-1);
            }

            Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * exp(Z * SCLDChainSum_.DXref(0));

            // y==y' components (i.e., Phi(0,y,y) term)
            Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -Term1;
         }
         else if (DOFderiv == 1)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-6.0 * SCLDChainSum_.phi2()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3());
            }
            else if (dt == PairPotentials::DT)
            {
               Term1 = (-6.0 * SCLDChainSum_.phi2T()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3T());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 1,D2T" << "\n";
               exit(-1);
            }

            if (x == 0)
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -Term1* PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX());
            }
            else if (x > 0)
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                          * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -Term1* OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1);
            }
         }
         else if (DOFderiv == 2)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-10.0 * SCLDChainSum_.phi3()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi4());
               Term2 = (-6.0 * SCLDChainSum_.phi2()
                        - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 2,DT,D2T" << "\n";
               exit(-1);
            }

            if ((x == 0) && (y == 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           + Term2 * PSI(SCLDChainSum_.pDX())) * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -(Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                            + Term2 * PSI(SCLDChainSum_.pDX()));
            }
            else if ((x == 0) && (y > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                           + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1))
                                                                          * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -(Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                            + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1));
            }
            else if ((y == 0) && (x > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                           + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1))
                                                                          * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -(Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                            + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1));
            }
            else if ((x > 0) && (y > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                           + Term2 * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1, y - 1))
                                                                          * exp(Z * SCLDChainSum_.DXref(0));

               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(0)] += -(Term1 * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                            * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                            + Term2 * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1, y - 1));
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 2" << "\n";
               exit(-1);
            }
         }
         else
         {
            cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 3" << "\n";
            exit(-1);
         }
      }
      else
      {
         // y = y' terms (i.e., diagonal terms)
         if (DOFderiv == 0)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2());
            }
            else if (dt == PairPotentials::DT)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1T() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2T());
            }
            else if (dt == PairPotentials::D2T)
            {
               Term1 = (-2.0 * SCLDChainSum_.phi1TT() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi2TT());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 0,D3T" << "\n";
               exit(-1);
            }

            Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
         }
         else if (DOFderiv == 1)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-6.0 * SCLDChainSum_.phi2() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3());
            }
            else if (dt == PairPotentials::DT)
            {
               Term1 = (-6.0 * SCLDChainSum_.phi2T() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3T());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 1,D2T" << "\n";
               exit(-1);
            }

            if (x == 0)
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
            else
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += Term1 * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                          * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
         }
         else if (DOFderiv == 2)
         {
            if (dt == PairPotentials::T0)
            {
               Term1 = (-10.0 * SCLDChainSum_.phi3() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi4());
               Term2 = (-6.0 * SCLDChainSum_.phi2() - 4.0 * SCLDChainSum_.Dx(0) * SCLDChainSum_.Dx(0) * SCLDChainSum_.phi3());
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 2,DT,D2T" << "\n";
               exit(-1);
            }

            if ((x == 0) && (y == 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           + Term2 * PSI(SCLDChainSum_.pDX())) * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
            else if ((x == 0) && (y > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                           + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1))
                                                                          * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
            else if ((y == 0) && (x > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                           + Term2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1))
                                                                          * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
            else if ((x > 0) && (y > 0))
            {
               Dk_static[SCLDChainSum_.Atom(0)][SCLDChainSum_.Atom(1)] += (Term1 * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1)
                                                                           * OMEGA(SCLDChainSum_.pDx(), SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), y - 1)
                                                                           + Term2 * SIGMA(SCLDChainSum_.Atom(0), SCLDChainSum_.Atom(1), x - 1, y - 1))
                                                                          * (exp(Z * SCLDChainSum_.DXref(0)) - 1.0);
            }
            else
            {
               cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 2" << "\n";
               exit(-1);
            }
         }
         else
         {
            cerr << "Error in SCLDQMultiChainTPP::ReferenceDynamicalStiffness - DOFderiv = 3" << "\n";
            exit(-1);
         }
      }
   }

   // Normalize through the Mass Matrix
   for (int p = 0; p < INTERNAL_ATOMS; ++p)
   {
      for (int q = 0; q < INTERNAL_ATOMS; ++q)
      {
         Dk_static[p][q] /= sqrt(AtomicMass_[p] * AtomicMass_[q]);
      }
   }

   return Dk_static;
}

void SCLDQMultiChainTPP::ReferenceDispersionCurves(Vector const& K, int const& NoPTS,
                                                   char const* const prefix,
                                                   ostream& out) const
{
   int w = out.width();
   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }

   double InverseLat;
   InverseLat = 1.0 / RefLattice_[0][0];

   Matrix EigVal[3];
   for (int i = 0; i < 3; ++i)
   {
      EigVal[i].Resize(1, INTERNAL_ATOMS);
   }

   double Z1, Z2;
   Z1 = K[0];
   Z2 = K[3];

   Z1 = InverseLat * Z1;
   Z2 = InverseLat * Z2;

   Vector Z(1);
   double
      DZ = Z2 - Z1;
   double dz = 1.0 / (NoPTS - 1);
   for (int k = 0; k < 2; ++k)
   {
      Z[0] = Z1 + (k * dz) * DZ;
      EigVal[k] = HermiteEigVal(ReferenceDynamicalStiffness(Z, PairPotentials::T0, 0, 0, 0));
      qsort(EigVal[k][0], INTERNAL_ATOMS, sizeof(double), &comp);

      out << prefix << setw(w) << k * dz;
      if (Echo_)
      {
         cout << prefix << setw(w) << k * dz;
      }
      for (int i = 0; i < INTERNAL_ATOMS; ++i)
      {
         out << setw(w) << EigVal[k][0][i];
         if (Echo_)
         {
            cout << setw(w) << EigVal[k][0][i];
         }
      }
      out << "\n";
      if (Echo_)
      {
         cout << "\n";
      }
   }
   int zero = 0, one = 1, two = 2;
   for (int k = 2; k < NoPTS; ++k)
   {
      Z[0] = Z1 + (k * dz) * DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z, PairPotentials::T0, 0, 0, 0));
      qsort(EigVal[two][0], INTERNAL_ATOMS, sizeof(double), &comp);
      interpolate(EigVal, zero, one, two);

      out << prefix << setw(w) << k * dz;
      if (Echo_)
      {
         cout << prefix << setw(w) << k * dz;
      }
      for (int i = 0; i < INTERNAL_ATOMS; ++i)
      {
         out << setw(w) << EigVal[two][0][i];
         if (Echo_)
         {
            cout << setw(w) << EigVal[two][0][i];
         }
      }
      out << "\n";
      if (Echo_)
      {
         cout << "\n";
      }

      zero = (zero + 1) % 3; one = (zero + 1) % 3; two = (one + 1) % 3;
   }
}

void SCLDQMultiChainTPP::ReferenceHarmonic() const
{
   if (FreqCached == 0)
   {
#ifdef _OPENMP
      double start_omp, end_omp;
      start_omp = omp_get_wtime();
#endif

      double Tol = 1.0e-12;

      double temp1, temp;
      MyComplexDouble temp2;

      int i;

      Matrix I;
      I.SetIdentity(INTERNAL_ATOMS);

      CMatrix A_static, A_DOF_static, A_DOF2_static, A_DOFDOF_static;

      Vector EigVec_norm;
      EigVec_norm.Resize(INTERNAL_ATOMS);

      CMatrix Eig, J, F, X;
      Eig.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
      J.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
      F.Resize(INTERNAL_ATOMS, 1, 0.0);
      X.Resize(INTERNAL_ATOMS, 1, 0.0);

      i = 0;
      for (ChainIter_.Reset(); !ChainIter_.Done(); ++ChainIter_)
      {
         Z_static[i].Resize(1);
         (Z_static[i])[0] = ChainIter_[0];
         i = i + 1;
      }

#pragma omp parallel for private(Eig,temp,temp1,temp2,EigVec_norm,J,F,X,A_static,A_DOF_static,A_DOF2_static,A_DOFDOF_static) schedule(dynamic)
      for (i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         // HEigVals_static and HEigVec_static
         A_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
         HEigVals_static[i].Resize(1, INTERNAL_ATOMS, 0.0);
         Eig.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
         EigVec_norm.Resize(INTERNAL_ATOMS);

#pragma omp critical
         {
            A_static = ReferenceDynamicalStiffness(Z_static[i], PairPotentials::T0, 0, 0, 0);
         }
         HEigVals_static[i] = HermiteEigVal(A_static, &(Eig));

         // sorting the eigenvalues in increasing order
         for (int k = 0; k < (INTERNAL_ATOMS - 1); ++k)
         {
            for (int p = k + 1; p < INTERNAL_ATOMS; ++p)
            {
               if ((HEigVals_static[i])[0][k] < (HEigVals_static[i])[0][p])
               {
                  temp1 = (HEigVals_static[i])[0][k];
                  (HEigVals_static[i])[0][k] = (HEigVals_static[i])[0][p];
                  (HEigVals_static[i])[0][p] = temp1;

                  for (int q = 0; q < INTERNAL_ATOMS; ++q)
                  {
                     temp2 = Eig[q][k];
                     Eig[q][k] = Eig[q][p];
                     Eig[q][p] = temp2;
                  }
               }
            }
         }

         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            (HEigVec_static[i])[k].Resize(INTERNAL_ATOMS, 1, 0.0);

            for (int p = 0; p < INTERNAL_ATOMS; ++p)
            {
               ((HEigVec_static[i])[k])[p][0] = Eig[p][k];
            }
         }


         // Clean up numerical round off (at least for zero values)
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            if (fabs((HEigVals_static[i])[0][k]) < Tol)
            {
               (HEigVals_static[i])[0][k] = 0.0;
            }
         }

         // checking whether any eigen values are repeated
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            for (int p = k + 1; p < INTERNAL_ATOMS; ++p)
            {
               if (fabs((HEigVals_static[i])[0][k] - (HEigVals_static[i])[0][p]) < Tol)
               {
                  cout << "Error in ReferenceHarmonic: Repeated eigenvalues" << endl;
                  cout << "HEigVals_static i = " << i << endl;
                  cout << setw(20) << HEigVals_static[i] << endl;
                  exit(-1);
               }
            }
         }

         // storing index of max abs. value in each HEigVec
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            temp = 0.0;
            for (int p = 0; p < INTERNAL_ATOMS; ++p)
            {
               if (temp < (((HEigVec_static[i])[k])[p][0]).mod())
               {
                  temp = (((HEigVec_static[i])[k])[p][0]).mod();
                  EigVec_norm[k] = p;
               }
            }
         }

         // Derivatives w.r.t DOF
         for (int x = 0; x < DOFS; ++x)
         {
            A_DOF_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
            (HEigValsDOF_static[x])[i].Resize(1, INTERNAL_ATOMS, 0.0);

#pragma omp critical
            {
               A_DOF_static = ReferenceDynamicalStiffness(Z_static[i], PairPotentials::T0, 1, x, 0);
            }

            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               J.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
               F.Resize(INTERNAL_ATOMS, 1, 0.0);
               X.Resize(INTERNAL_ATOMS, 1, 0.0);
               ((HEigVecDOF_static[i])[k])[x].Resize(INTERNAL_ATOMS, 1, 0.0);

               ((HEigValsDOF_static[x])[i])[0][k] = (((HEigVec_static[i])[k].ConjTrans() * A_DOF_static * (HEigVec_static[i])[k])[0][0]).real();

               F = -(A_DOF_static - ((HEigValsDOF_static[x])[i])[0][k] * I) * (HEigVec_static[i])[k];
               J = A_static - (HEigVals_static[i])[0][k] * I;

               for (int p = 0; p < INTERNAL_ATOMS; ++p)
               {
                  J[p][int(EigVec_norm[k])] = 0.0;
                  J[int(EigVec_norm[k])][p] = 0.0;
               }
               J[int(EigVec_norm[k])][int(EigVec_norm[k])] = 1.0;
               F[int(EigVec_norm[k])][0] = 0.0;

               X = SolvePLU(J, F);

               temp = -0.5 * ((X.ConjTrans() * (HEigVec_static[i])[k] + (HEigVec_static[i])[k].ConjTrans() * X)[0][0]).real();

               ((HEigVecDOF_static[i])[k])[x] = X + temp * (HEigVec_static[i])[k];
            }
         }


         // HEigValsDOFDOF_static and HEigVecDOFDOF_static
         for (int x = 0; x < DOFS; ++x)
         {
            A_DOF_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);

#pragma omp critical
            {
               A_DOF_static = ReferenceDynamicalStiffness(Z_static[i], PairPotentials::T0, 1, x, 0);
            }

            for (int y = x; y < DOFS; ++y)
            {
               A_DOF2_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);

#pragma omp critical
               {
                  A_DOF2_static = ReferenceDynamicalStiffness(Z_static[i], PairPotentials::T0, 1, y, 0);
               }

               A_DOFDOF_static.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
               ((HEigValsDOFDOF_static[x])[y])[i].Resize(1, INTERNAL_ATOMS, 0.0);

#pragma omp critical
               {
                  A_DOFDOF_static = ReferenceDynamicalStiffness(Z_static[i], PairPotentials::T0, 2, x, y);
               }

               for (int k = 0; k < INTERNAL_ATOMS; ++k)
               {
                  J.Resize(INTERNAL_ATOMS, INTERNAL_ATOMS, 0.0);
                  F.Resize(INTERNAL_ATOMS, 1, 0.0);
                  X.Resize(INTERNAL_ATOMS, 1, 0.0);
                  (((HEigVecDOFDOF_static[i])[k])[x])[y].Resize(INTERNAL_ATOMS, 1, 0.0);

                  (((HEigValsDOFDOF_static[x])[y])[i])[0][k]
                     = ((
                           (HEigVec_static[i])[k].ConjTrans() * A_DOFDOF_static * (HEigVec_static[i])[k]
                           + (HEigVec_static[i])[k].ConjTrans() * (A_DOF_static - ((HEigValsDOF_static[x])[i])[0][k] * I) * ((HEigVecDOF_static[i])[k])[y]
                           + (HEigVec_static[i])[k].ConjTrans() * (A_DOF2_static - ((HEigValsDOF_static[y])[i])[0][k] * I) * ((HEigVecDOF_static[i])[k])[x]
                        )[0][0]).real();

                  F = -(
                     (A_DOFDOF_static - (((HEigValsDOFDOF_static[x])[y])[i])[0][k] * I) * (HEigVec_static[i])[k]
                     + (A_DOF_static - ((HEigValsDOF_static[x])[i])[0][k] * I) * ((HEigVecDOF_static[i])[k])[y]
                     + (A_DOF2_static - ((HEigValsDOF_static[y])[i])[0][k] * I) * ((HEigVecDOF_static[i])[k])[x]
                       );

                  J = A_static - (HEigVals_static[i])[0][k] * I;

                  for (int p = 0; p < INTERNAL_ATOMS; ++p)
                  {
                     J[p][int(EigVec_norm[k])] = 0.0;
                     J[int(EigVec_norm[k])][p] = 0.0;
                  }
                  J[int(EigVec_norm[k])][int(EigVec_norm[k])] = 1.0;
                  F[int(EigVec_norm[k])][0] = 0.0;

                  X = SolvePLU(J, F);

                  temp = -0.5 * ((
                                    X.ConjTrans() * (HEigVec_static[i])[k] + (HEigVec_static[i])[k].ConjTrans() * X
                                    + ((HEigVecDOF_static[i])[k])[x].ConjTrans() * ((HEigVecDOF_static[i])[k])[y]
                                    + ((HEigVecDOF_static[i])[k])[y].ConjTrans() * ((HEigVecDOF_static[i])[k])[x]
                                 )[0][0]).real();

                  (((HEigVecDOFDOF_static[i])[k])[x])[y] = X + temp * (HEigVec_static[i])[k];
               }
            }
         }
      }


      //    cout <<"HEigVals_static:"<<endl;
      //    i = 0;
      //    for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //    {
      //        for (int k=0;k<INTERNAL_ATOMS;++k)
      //        {
      //            cout << setw(20) << (HEigVals_static[i])[0][k] ;
      //        }
      //        cout << endl;
      //        i = i+1;
      //    }

      //    for (int x=0;x<DOFS;++x)
      //    {
      //        cout <<"HEigValsDOF_static: x="<<x<<endl;
      //        i = 0;
      //        for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //        {
      //            for (int k=0;k<INTERNAL_ATOMS;++k)
      //            {
      //                cout << setw(20) << ((HEigValsDOF_static[x])[i])[0][k] ;
      //            }
      //            cout << endl;
      //
      //            i = i+1;
      //        }
      //        cout << endl;
      //    }

      //    for (int x=0;x<DOFS;++x)
      //    {
      //    for (int y=x;y<DOFS;++y)
      //    {
      //        cout <<"HEigValsDOFDOF_static: x="<<x<<" y="<<y<<endl;
      //        i = 0;
      //        for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //        {
      //            for (int k=0;k<INTERNAL_ATOMS;++k)
      //            {
      //                cout << setw(20) << (((HEigValsDOFDOF_static[x])[y])[i])[0][k] ;
      //            }
      //            cout << endl;
      //            i = i+1;
      //        }
      //        cout << endl;
      //    }
      //    }

#ifdef _OPENMP
      end_omp = omp_get_wtime();
      cout << "Harmonic Total time OMP: " << (start_omp - end_omp) * CLOCKS_PER_SEC << endl;
#endif
   }
}

void SCLDQMultiChainTPP::ReferenceV4() const
{
   if (FreqCached == 0)
   {
#ifdef _OPENMP
      double start_omp, end_omp;
      start_omp = omp_get_wtime();
#endif

      double pi = 4.0 * atan(1.0);
      MyComplexDouble Ic(0, 1);
      double InverseLat = 1.0 / RefLattice_[0][0];
      MyComplexDouble A = (2.0 * pi - 2.0 * pi / 2.5) * InverseLat;
      MyComplexDouble B = (2.0 * pi / 4.0) * InverseLat;

      double Mass0, Mass1;
      double M1, M2, M3, M4;

      double DXref, Dx;
      int Atom0, Atom1;


      double PTerm1, PTerm2, PTerm3;

      // Calculates V4_static, V4DOF_static, and V4DOFDOF_static for different values of k,\nu,k', and \nu'

      // Initialization
      V4_static.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);

      for (int x = 0; x < DOFS; ++x)
      {
         (V4DOF_static[x]).Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);

         for (int y = x; y < DOFS; ++y)
         {
            ((V4DOFDOF_static[x])[y]).Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);
         }
      }


      for (SCLDChainSum_.Reset(); !SCLDChainSum_.Done(); ++SCLDChainSum_)
      {
         DXref = SCLDChainSum_.DXref(0);
         Dx = SCLDChainSum_.Dx(0);

         Atom0 = SCLDChainSum_.Atom(0);
         Atom1 = SCLDChainSum_.Atom(1);

         Mass0 = AtomicMass_[Atom0];
         Mass1 = AtomicMass_[Atom1];

         M1 = 1 / (Mass0 * Mass0);
         M2 = 1 / (sqrt(Mass0 * Mass1) * Mass1);
         M3 = 1 / (sqrt(Mass0 * Mass1) * Mass0);
         M4 = 1 / (Mass0 * Mass1);

         PTerm1 = 16 * SCLDChainSum_.phi4() * Dx * Dx * Dx * Dx
                  + 48 * SCLDChainSum_.phi3() * Dx * Dx
                  + 12 * SCLDChainSum_.phi2();

         PTerm2 = 16 * SCLDChainSum_.phi5() * Dx * Dx * Dx * Dx
                  + 80 * SCLDChainSum_.phi4() * Dx * Dx
                  + 60 * SCLDChainSum_.phi3();

         PTerm3 = 16 * SCLDChainSum_.phi6() * Dx * Dx * Dx * Dx
                  + 112 * SCLDChainSum_.phi5() * Dx * Dx
                  + 140 * SCLDChainSum_.phi4();

         // #pragma omp parallel for schedule(dynamic) // collapse only works when gcc version is greater than 4.3 or something..it works on "bop.aem.umn.edu" but not "soul.aem.umn.edu" machine
#pragma omp parallel for schedule(dynamic) collapse(4)
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
               {
                  for (int l = 0; l < INTERNAL_ATOMS; ++l)
                  {
                     MyComplexDouble V1, V2, V3, V4, V5, V6, V7, V8;
                     MyComplexDouble V1DOF, V2DOF, V3DOF, V4DOF, V5DOF, V6DOF, V7DOF, V8DOF;
                     MyComplexDouble V1DOF2, V2DOF2, V3DOF2, V4DOF2, V5DOF2, V6DOF2, V7DOF2, V8DOF2;
                     MyComplexDouble V1DOFDOF, V2DOFDOF, V3DOFDOF, V4DOFDOF, V5DOFDOF, V6DOFDOF, V7DOFDOF, V8DOFDOF;

                     MyComplexDouble Term1, Term2, Term22, Term3;

                     MyComplexDouble Z, ZNew;
                     MyComplexDouble Exp1, Exp2;
                     double Iter, IterNew;

                     Iter = (Z_static[i])[0];
                     Z = (A * Iter + B) * Ic;
                     Exp1 = exp(-Z * DXref);

                     IterNew = (Z_static[j])[0];
                     ZNew = (A * IterNew + B) * Ic;
                     Exp2 = exp(-ZNew * DXref);

                     if ((i * INTERNAL_ATOMS + k) >= (j * INTERNAL_ATOMS + l))
                     {
                        V1 = ((HEigVec_static[i])[k])[Atom0][0];
                        V2 = ((HEigVec_static[i])[k])[Atom1][0];
                        V3 = ((HEigVec_static[j])[l])[Atom0][0];
                        V4 = ((HEigVec_static[j])[l])[Atom1][0];
                        V5 = V1.conj();
                        V6 = V2.conj();
                        V7 = V3.conj();
                        V8 = V4.conj();

                        Term1 = V1 *
                                (
                           V5 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                           + V6 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                );

                        V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                           = V4_static[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                += (PTerm1 * Term1.real());

                        for (int x = 0; x < DOFS; ++x)
                        {
                           V1DOF = (((HEigVecDOF_static[i])[k])[x])[Atom0][0];
                           V2DOF = (((HEigVecDOF_static[i])[k])[x])[Atom1][0];
                           V3DOF = (((HEigVecDOF_static[j])[l])[x])[Atom0][0];
                           V4DOF = (((HEigVecDOF_static[j])[l])[x])[Atom1][0];
                           V5DOF = V1DOF.conj();
                           V6DOF = V2DOF.conj();
                           V7DOF = V3DOF.conj();
                           V8DOF = V4DOF.conj();

                           Term2 = V1DOF *
                                   (
                              V5 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                              + V6 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                   )
                                   +
                                   V1 *
                                   (
                              V5DOF * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                              + V5 * (
                                 V3DOF * (M1 * V7 - M3 * V8 * Exp2) + V4DOF * (M4 * V8 - M3 * V7 * Exp2.conj())
                                 + V3 * (M1 * V7DOF - M3 * V8DOF * Exp2) + V4 * (M4 * V8DOF - M3 * V7DOF * Exp2.conj())
                                     )
                              + V6DOF * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                              + V6 * (
                                 V3DOF * (M4 * V8 * Exp2 - M3 * V7) + V4DOF * (M4 * V7 * Exp2.conj() - M2 * V8)
                                 + V3 * (M4 * V8DOF * Exp2 - M3 * V7DOF) + V4 * (M4 * V7DOF * Exp2.conj() - M2 * V8DOF)
                                     ) * Exp1
                                   );


                           if (x == 0)
                           {
                              (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                 = (V4DOF_static[x])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                      += PTerm2 * Term1.real() * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) + PTerm1* Term2.real();
                           }
                           else
                           {
                              (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                 = (V4DOF_static[x])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                      += PTerm2 * Term1.real() * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, x - 1) + PTerm1* Term2.real();
                           }

                           for (int y = x; y < DOFS; ++y)
                           {
                              V1DOF2 = (((HEigVecDOF_static[i])[k])[y])[Atom0][0];
                              V2DOF2 = (((HEigVecDOF_static[i])[k])[y])[Atom1][0];
                              V3DOF2 = (((HEigVecDOF_static[j])[l])[y])[Atom0][0];
                              V4DOF2 = (((HEigVecDOF_static[j])[l])[y])[Atom1][0];
                              V5DOF2 = V1DOF2.conj();
                              V6DOF2 = V2DOF2.conj();
                              V7DOF2 = V3DOF2.conj();
                              V8DOF2 = V4DOF2.conj();

                              Term22 = V1DOF2 *
                                       (
                                 V5 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V6 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                       )
                                       +
                                       V1 *
                                       (
                                 V5DOF2 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V5 * (
                                    V3DOF2 * (M1 * V7 - M3 * V8 * Exp2) + V4DOF2 * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3 * (M1 * V7DOF2 - M3 * V8DOF2 * Exp2) + V4 * (M4 * V8DOF2 - M3 * V7DOF2 * Exp2.conj())
                                        )
                                 + V6DOF2 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                 + V6 * (
                                    V3DOF2 * (M4 * V8 * Exp2 - M3 * V7) + V4DOF2 * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3 * (M4 * V8DOF2 * Exp2 - M3 * V7DOF2) + V4 * (M4 * V7DOF2 * Exp2.conj() - M2 * V8DOF2)
                                        ) * Exp1
                                       );

                              V1DOFDOF = ((((HEigVecDOFDOF_static[i])[k])[x])[y])[Atom0][0];
                              V2DOFDOF = ((((HEigVecDOFDOF_static[i])[k])[x])[y])[Atom1][0];
                              V3DOFDOF = ((((HEigVecDOFDOF_static[j])[l])[x])[y])[Atom0][0];
                              V4DOFDOF = ((((HEigVecDOFDOF_static[j])[l])[x])[y])[Atom1][0];
                              V5DOFDOF = V1DOFDOF.conj();
                              V6DOFDOF = V2DOFDOF.conj();
                              V7DOFDOF = V3DOFDOF.conj();
                              V8DOFDOF = V4DOFDOF.conj();

                              Term3 = V1DOFDOF *
                                      (
                                 V5 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V6 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                      )
                                      +
                                      V1DOF2 *
                                      (
                                 V5DOF * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V5 * (
                                    V3DOF * (M1 * V7 - M3 * V8 * Exp2) + V4DOF * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3 * (M1 * V7DOF - M3 * V8DOF * Exp2) + V4 * (M4 * V8DOF - M3 * V7DOF * Exp2.conj())
                                        )
                                 + V6DOF * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                 + V6 * (
                                    V3DOF * (M4 * V8 * Exp2 - M3 * V7) + V4DOF * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3 * (M4 * V8DOF * Exp2 - M3 * V7DOF) + V4 * (M4 * V7DOF * Exp2.conj() - M2 * V8DOF)
                                        ) * Exp1
                                      )
                                      + V1DOF *
                                      (
                                 V5DOF2 * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V5 * (
                                    V3DOF2 * (M1 * V7 - M3 * V8 * Exp2) + V4DOF2 * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3 * (M1 * V7DOF2 - M3 * V8DOF2 * Exp2) + V4 * (M4 * V8DOF2 - M3 * V7DOF2 * Exp2.conj())
                                        )
                                 + V6DOF2 * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                 + V6 * (
                                    V3DOF2 * (M4 * V8 * Exp2 - M3 * V7) + V4DOF2 * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3 * (M4 * V8DOF2 * Exp2 - M3 * V7DOF2) + V4 * (M4 * V7DOF2 * Exp2.conj() - M2 * V8DOF2)
                                        ) * Exp1
                                      )
                                      +
                                      V1 *
                                      (
                                 V5DOFDOF * (V3 * (M1 * V7 - M3 * V8 * Exp2) + V4 * (M4 * V8 - M3 * V7 * Exp2.conj()))
                                 + V5DOF * (
                                    V3DOF2 * (M1 * V7 - M3 * V8 * Exp2) + V4DOF2 * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3 * (M1 * V7DOF2 - M3 * V8DOF2 * Exp2) + V4 * (M4 * V8DOF2 - M3 * V7DOF2 * Exp2.conj())
                                           )
                                 + V5DOF2 * (
                                    V3DOF * (M1 * V7 - M3 * V8 * Exp2) + V4DOF * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3 * (M1 * V7DOF - M3 * V8DOF * Exp2) + V4 * (M4 * V8DOF - M3 * V7DOF * Exp2.conj())
                                            )
                                 + V5 * (
                                    V3DOFDOF * (M1 * V7 - M3 * V8 * Exp2) + V4DOFDOF * (M4 * V8 - M3 * V7 * Exp2.conj())
                                    + V3DOF * (M1 * V7DOF2 - M3 * V8DOF2 * Exp2) + V4DOF * (M4 * V8DOF2 - M3 * V7DOF2 * Exp2.conj())
                                    + V3DOF2 * (M1 * V7DOF - M3 * V8DOF * Exp2) + V4DOF2 * (M4 * V8DOF - M3 * V7DOF * Exp2.conj())
                                    + V3 * (M1 * V7DOFDOF - M3 * V8DOFDOF * Exp2) + V4 * (M4 * V8DOFDOF - M3 * V7DOFDOF * Exp2.conj())
                                        )
                                 + V6DOFDOF * (V3 * (M4 * V8 * Exp2 - M3 * V7) + V4 * (M4 * V7 * Exp2.conj() - M2 * V8)) * Exp1
                                 + V6DOF * (
                                    V3DOF2 * (M4 * V8 * Exp2 - M3 * V7) + V4DOF2 * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3 * (M4 * V8DOF2 * Exp2 - M3 * V7DOF2) + V4 * (M4 * V7DOF2 * Exp2.conj() - M2 * V8DOF2)
                                           ) * Exp1
                                 + V6DOF2 * (
                                    V3DOF * (M4 * V8 * Exp2 - M3 * V7) + V4DOF * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3 * (M4 * V8DOF * Exp2 - M3 * V7DOF) + V4 * (M4 * V7DOF * Exp2.conj() - M2 * V8DOF)
                                            ) * Exp1
                                 + V6 * (
                                    V3DOFDOF * (M4 * V8 * Exp2 - M3 * V7) + V4DOFDOF * (M4 * V7 * Exp2.conj() - M2 * V8)
                                    + V3DOF * (M4 * V8DOF2 * Exp2 - M3 * V7DOF2) + V4DOF * (M4 * V7DOF2 * Exp2.conj() - M2 * V8DOF2)
                                    + V3DOF2 * (M4 * V8DOF * Exp2 - M3 * V7DOF) + V4DOF2 * (M4 * V7DOF * Exp2.conj() - M2 * V8DOF)
                                    + V3 * (M4 * V8DOFDOF * Exp2 - M3 * V7DOFDOF) + V4 * (M4 * V7DOFDOF * Exp2.conj() - M2 * V8DOFDOF)
                                        ) * Exp1
                                      );


                              if ((x == 0) && (y == 0))
                              {
                                 ((V4DOFDOF_static[x])[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                    = ((V4DOFDOF_static[x])[y])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                         += (PTerm3 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX()) * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                             + PTerm2 * PSI(SCLDChainSum_.pDX())) * Term1.real()
                                            + PTerm2* Term2.real() * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                            + PTerm2* Term22.real() * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                            + PTerm1* Term3.real();
                              }
                              else if ((x == 0) && (y > 0))
                              {
                                 ((V4DOFDOF_static[x])[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                    = ((V4DOFDOF_static[x])[y])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                         += (PTerm3 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                             * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, y - 1)
                                             + PTerm2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), Atom0, Atom1, y - 1))
                                            * Term1.real()
                                            + PTerm2* Term2.real() * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, y - 1)
                                            + PTerm2* Term22.real() * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                            + PTerm1* Term3.real();
                              }
                              else if ((y == 0) && (x > 0))
                              {
                                 ((V4DOFDOF_static[x])[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                    = ((V4DOFDOF_static[x])[y])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                         += (PTerm3 * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                             * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, x - 1)
                                             + PTerm2 * GAMMA(SCLDChainSum_.pDx(), SCLDChainSum_.pDX(), Atom0, Atom1, x - 1))
                                            * Term1.real()
                                            + PTerm2* Term22.real() * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, x - 1)
                                            + PTerm2* Term2.real() * PI(SCLDChainSum_.pDx(), SCLDChainSum_.pDX())
                                            + PTerm1* Term3.real();
                              }
                              else if ((x > 0) && (y > 0))
                              {
                                 ((V4DOFDOF_static[x])[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                    = ((V4DOFDOF_static[x])[y])[j * INTERNAL_ATOMS + l][i * INTERNAL_ATOMS + k]
                                         += (PTerm3
                                             * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, y - 1)
                                             * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, x - 1)
                                             + PTerm2 * SIGMA(Atom0, Atom1, x - 1, y - 1))
                                            * Term1.real()
                                            + PTerm2* Term22.real() * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, x - 1)
                                            + PTerm2* Term2.real() * OMEGA(SCLDChainSum_.pDx(), Atom0, Atom1, y - 1)
                                            + PTerm1* Term3.real();
                              }
                              else
                              {
                                 cerr << "Error in SCLDQMultiChainTPP::ReferenceV4 - DOFderiv = 3" << "\n";
                                 exit(-1);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }


      // V4 = V4 / (GridSize_+2)
      V4_static /= (GridSize_ + 2);

      for (int x = 0; x < DOFS; ++x)
      {
         V4DOF_static[x] /= (GridSize_ + 2);

         for (int y = x; y < DOFS; ++y)
         {
            (V4DOFDOF_static[x])[y] /= (GridSize_ + 2);
         }
      }


#ifdef _OPENMP
      end_omp = omp_get_wtime();
      cout << "V4 Total time OMP: " << (start_omp - end_omp) * CLOCKS_PER_SEC << endl;
#endif
   }
}

void SCLDQMultiChainTPP::ReferencePseudoHarmonic() const
{
   if (FreqCached == 0)
   {
#ifdef _OPENMP
      double start_omp, end_omp;
      start_omp = omp_get_wtime();
#endif

      int converged = 0;
      double Tol = 1.0e-10;

      int MaxCounter = 2000;
      int counter;

      Matrix F, J, X, X_previous;

      //  Finding EigVals_static
      F.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);
      J.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);
      X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);
      X_previous.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

      // Copying HEigVals_static into X_previous to start the iteration
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         EigVals_static[i].Resize(1, INTERNAL_ATOMS, 0.0);

         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            if (FirstConverged == 1)
            {
               X_previous[i * INTERNAL_ATOMS + k][0] = (EigVals_previous_static[i])[0][k];
            }
            else if (FirstConverged == 0)
            {
               X_previous[i * INTERNAL_ATOMS + k][0] = (HEigVals_static[i])[0][k];
            }
            else
            {
               cout << "Error in initiating the X_previous: EigVals" << endl;
               exit(-1);
            }
         }
      }

      // solving self-consistent equations using iterations
      counter = 0;
      converged = 0;
      while (converged == 0)
      {
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               F[i * INTERNAL_ATOMS + k][0] = -X_previous[i * INTERNAL_ATOMS + k][0] + (HEigVals_static[i])[0][k];

               for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
               {
                  for (int l = 0; l < INTERNAL_ATOMS; ++l)
                  {
                     if (X_previous[j * INTERNAL_ATOMS + l][0] > Tol)
                     {
                        F[i * INTERNAL_ATOMS + k][0] += 2.0 * 0.5 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / X_previous[j * INTERNAL_ATOMS + l][0];
                     }
                     else
                     {
                        F[i * INTERNAL_ATOMS + k][0] += 0.0;
                     }

                     if ((i == j) && (k == l))
                     {
                        if (X_previous[j * INTERNAL_ATOMS + l][0] > Tol)
                        {
                           J[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -1.0 - 2.0 * 0.5 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (X_previous[j * INTERNAL_ATOMS + l][0] * X_previous[j * INTERNAL_ATOMS + l][0]);
                        }
                        else
                        {
                           J[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -1.0;
                        }
                     }
                     else
                     {
                        if (X_previous[j * INTERNAL_ATOMS + l][0] > Tol)
                        {
                           J[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -2.0 * 0.5 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (X_previous[j * INTERNAL_ATOMS + l][0] * X_previous[j * INTERNAL_ATOMS + l][0]);
                        }
                        else
                        {
                           J[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = 0.0;
                        }
                     }
                  }
               }
            }
         }

         F *= -1.0;
         X = SolvePLU(J, F);
         X += X_previous;

         // checking convergence for next iteration
         for (int p = 0; p < (GridSize_ / 2 + 1) * INTERNAL_ATOMS; ++p)
         {
            if (fabs(X[p][0] - X_previous[p][0]) > Tol)
            {
               converged = 0;
               break;
            }
            else
            {
               converged = 1;
            }
         }

         // copying X to X-previous for next iteration
         for (int p = 0; p < (GridSize_ / 2 + 1) * INTERNAL_ATOMS; ++p)
         {
            X_previous[p][0] = X[p][0];
         }

         counter = counter + 1;
         if (counter > MaxCounter)
         {
            cout << "Error: ReferencePseudoHarmonic: EigVals_static convergence has surpassed Max. number of iterations" << endl;
            exit(-1);
         }
      }

      // Copying the converged X values into EigVals_static
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            (EigVals_static[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];

            // Clean up numerical round off (at least for zero values)
            if (fabs((EigVals_static[i])[0][k]) < Tol)
            {
               (EigVals_static[i])[0][k] = 0.0;
            }
         }
      }

      //    cout <<"EigVals_static:"<<endl;
      //    for (int i=0;i<(GridSize_/2+1);++i)
      //    {
      //        for (int k=0;k<INTERNAL_ATOMS;++k)
      //        {
      //            cout << setw(20) << (EigVals_static[i])[0][k] ;
      //        }
      //        cout << endl;
      //    }



      Matrix A, P, L, U, B;

      //  Finding EigValsT_static
      A.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);
      B.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

      // Finding A and B for the matrix equation A X = B
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            // B[i*INTERNAL_ATOMS+k][0] = - (HEigValsT_static[i])[0][k];
            B[i * INTERNAL_ATOMS + k][0] = -0.0;

            for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
            {
               for (int l = 0; l < INTERNAL_ATOMS; ++l)
               {
                  if ((EigVals_static[j])[0][l] > Tol)
                  {
                     B[i * INTERNAL_ATOMS + k][0] -= 2.0 * 0.5 * kB_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (EigVals_static[j])[0][l];
                  }
                  else
                  {
                     B[i * INTERNAL_ATOMS + k][0] -= 0.0;
                  }

                  if ((i == j) && (k == l))
                  {
                     if ((EigVals_static[j])[0][l] > Tol)
                     {
                        A[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -1.0 - 2.0 * 0.5 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                                                            / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l]);
                     }
                     else
                     {
                        A[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -1.0;
                     }
                  }
                  else
                  {
                     if ((EigVals_static[j])[0][l] > Tol)
                     {
                        A[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = -2.0 * 0.5 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l]
                                                                            / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l]);
                     }
                     else
                     {
                        A[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] = 0.0;
                     }
                  }
               }
            }
         }
      }

      P.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);
      L.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);
      U.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, (GridSize_ / 2 + 1) * INTERNAL_ATOMS, 0.0);

      PLU(A, P, L, U);

      X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

      X = SolvePLU(P, L, U, B);

      // Copying the converged X values into EigValsT_static
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         EigValsT_static[i].Resize(1, INTERNAL_ATOMS, 0.0);
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            (EigValsT_static[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];
         }
      }

      //    cout <<"EigValsT_static:"<<endl;
      //    for (int i=0;i<(GridSize_/2+1);++i)
      //    {
      //        for (int k=0;k<INTERNAL_ATOMS;++k)
      //        {
      //            cout << setw(20) << (EigValsT_static[i])[0][k] ;
      //        }
      //        cout << endl;
      //    }


      //  Finding EigValsTT_static
      //    A.Resize((GridSize_/2+1)*INTERNAL_ATOMS,(GridSize_/2+1)*INTERNAL_ATOMS,0.0);
      B.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

      // Finding B for the matrix equation A X = B
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            // B[i*INTERNAL_ATOMS+k][0] = - (HEigValsTT_static[i])[0][k];
            B[i * INTERNAL_ATOMS + k][0] = -0.0;

            for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
            {
               for (int l = 0; l < INTERNAL_ATOMS; ++l)
               {
                  if ((EigVals_static[j])[0][l] > Tol)
                  {
                     B[i * INTERNAL_ATOMS + k][0] -= 2.0 * (
                        -0.5 * kB_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * (EigValsT_static[j])[0][l]
                        / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                        - 0.5 * kB_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * (EigValsT_static[j])[0][l]
                        / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                        + 0.5 * 2.0 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * (EigValsT_static[j])[0][l] * (EigValsT_static[j])[0][l]
                        / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                                                           );
                  }
                  else
                  {
                     B[i * INTERNAL_ATOMS + k][0] -= 0.0;
                  }
               }
            }
         }
      }

      X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

      X = SolvePLU(P, L, U, B);

      // Copying the converged X values into EigValsTT_static
      for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
      {
         EigValsTT_static[i].Resize(1, INTERNAL_ATOMS, 0.0);
         for (int k = 0; k < INTERNAL_ATOMS; ++k)
         {
            (EigValsTT_static[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];
         }
      }

      //    cout <<"EigValsTT_static:"<<endl;
      //    i = 0;
      //    for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //    {
      //        for (int k=0;k<INTERNAL_ATOMS;++k)
      //        {
      //            cout << setw(20) << (EigValsTT_static[i])[0][k] ;
      //        }
      //        cout << endl;
      //        i = i+1;
      //    }



      //  Finding EigValsDOF_static
#pragma omp parallel for private(B,X) schedule(dynamic)
      for (int x = 0; x < DOFS; ++x)
      {
         B.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

         // Finding B for the matrix equation A X = B
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               B[i * INTERNAL_ATOMS + k][0] = -((HEigValsDOF_static[x])[i])[0][k];

               for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
               {
                  for (int l = 0; l < INTERNAL_ATOMS; ++l)
                  {
                     if ((EigVals_static[j])[0][l] > Tol)
                     {
                        B[i * INTERNAL_ATOMS + k][0] -= 2.0 * 0.5 * kB_ * NTemp_ * (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (EigVals_static[j])[0][l];
                     }
                     else
                     {
                        B[i * INTERNAL_ATOMS + k][0] -= 0.0;
                     }
                  }
               }
            }
         }

         X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

         X = SolvePLU(P, L, U, B);

         // Copying the converged X values into EigValsDOF_static
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            (EigValsDOF_static[x])[i].Resize(1, INTERNAL_ATOMS, 0.0);
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               ((EigValsDOF_static[x])[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];
            }
         }
      }

      //    for (int x=0;x<DOFS;++x)
      //    {
      //        cout <<"EigValsDOF_static: x="<<x<<endl;
      //        i = 0;
      //        for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //        {
      //            for (int k=0;k<INTERNAL_ATOMS;++k)
      //            {
      //                cout << setw(20) << ((EigValsDOF_static[x])[i])[0][k] ;
      //            }
      //            cout << endl;
      //            i = i+1;
      //        }
      //        cout << endl;
      //    }



      //  Finding EigValsTDOF_static
#pragma omp parallel for private(B,X) schedule(dynamic)
      for (int x = 0; x < DOFS; ++x)
      {
         B.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

         // Finding B for the matrix equation A X = B
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               // B[i*INTERNAL_ATOMS+k][0] = - ((HEigValsTDOF_static[x])[i])[0][k];
               B[i * INTERNAL_ATOMS + k][0] = -0.0;

               for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
               {
                  for (int l = 0; l < INTERNAL_ATOMS; ++l)
                  {
                     if ((EigVals_static[j])[0][l] > Tol)
                     {
                        B[i * INTERNAL_ATOMS + k][0] -= 2.0 * (
                           0.5 * kB_ * (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (EigVals_static[j])[0][l]
                           - 0.5 * kB_ * NTemp_ * (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * (EigValsT_static[j])[0][l]
                           / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                           - 0.5 * kB_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * ((EigValsDOF_static[x])[j])[0][l]
                           / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                           + 0.5 * 2.0 * kB_ * NTemp_ * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * ((EigValsDOF_static[x])[j])[0][l] * (EigValsT_static[j])[0][l]
                           / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                                                              );
                     }
                     else
                     {
                        B[i * INTERNAL_ATOMS + k][0] -= 0.0;
                     }
                  }
               }
            }
         }

         X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

         X = SolvePLU(P, L, U, B);

         // Copying the converged X values into EigValsDOF_static
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            (EigValsTDOF_static[x])[i].Resize(1, INTERNAL_ATOMS, 0.0);
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               ((EigValsTDOF_static[x])[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];
            }
         }
      }

      //    for (int x=0;x<DOFS;++x)
      //    {
      //        cout <<"EigValsTDOF_static: x="<<x<<endl;
      //        i = 0;
      //        for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
      //        {
      //            for (int k=0;k<INTERNAL_ATOMS;++k)
      //            {
      //                cout << setw(20) << ((EigValsTDOF_static[x])[i])[0][k] ;
      //            }
      //            cout << endl;
      //            i = i+1;
      //        }
      //        cout << endl;
      //    }




      //  Finding EigValsDOFDOF_static
#pragma omp parallel for private(B,X) schedule(dynamic)
      for (int x = 0; x < DOFS; ++x)
      {
         for (int y = x; y < DOFS; ++y)
         {
            B.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

            // Finding B for the matrix equation A X = B
            for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
            {
               for (int k = 0; k < INTERNAL_ATOMS; ++k)
               {
                  B[i * INTERNAL_ATOMS + k][0] = -(((HEigValsDOFDOF_static[x])[y])[i])[0][k];

                  for (int j = 0; j < (GridSize_ / 2 + 1); ++j)
                  {
                     for (int l = 0; l < INTERNAL_ATOMS; ++l)
                     {
                        if ((EigVals_static[j])[0][l] > Tol)
                        {
                           B[i * INTERNAL_ATOMS + k][0] -= 2.0 * (
                              0.5 * kB_ * NTemp_ * ((V4DOFDOF_static[x])[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] / (EigVals_static[j])[0][l]
                              - 0.5 * kB_ * NTemp_ * (V4DOF_static[x])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * ((EigValsDOF_static[y])[j])[0][l]
                              / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                              - 0.5 * kB_ * NTemp_ * (V4DOF_static[y])[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * ((EigValsDOF_static[x])[j])[0][l]
                              / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                              + 0.5 * kB_ * NTemp_ * 2.0 * V4_static[i * INTERNAL_ATOMS + k][j * INTERNAL_ATOMS + l] * ((EigValsDOF_static[x])[j])[0][l] * ((EigValsDOF_static[y])[j])[0][l]
                              / ((EigVals_static[j])[0][l] * (EigVals_static[j])[0][l] * (EigVals_static[j])[0][l])
                                                                 );
                        }
                        else
                        {
                           B[i * INTERNAL_ATOMS + k][0] -= 0.0;
                        }
                     }
                  }
               }
            }

            X.Resize((GridSize_ / 2 + 1) * INTERNAL_ATOMS, 1, 0.0);

            X = SolvePLU(P, L, U, B);

            // Copying the converged X values into EigValsDOFDOF_static
            for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
            {
               ((EigValsDOFDOF_static[x])[y])[i].Resize(1, INTERNAL_ATOMS, 0.0);

               for (int k = 0; k < INTERNAL_ATOMS; ++k)
               {
                  (((EigValsDOFDOF_static[x])[y])[i])[0][k] = X[i * INTERNAL_ATOMS + k][0];
               }
            }
         }
      }



#ifdef _OPENMP
      end_omp = omp_get_wtime();
      cout << "SC Total time OMP: " << (start_omp - end_omp) * CLOCKS_PER_SEC << endl;
#endif
   }
   FreqCached = 1;
}

int SCLDQMultiChainTPP::ReferenceBlochWave(Vector& K) const
{
   int i = 0;
   int NumOfNegEV = 0;
   double MinEigVal;
   double Tol = 1.0e-13;

   MinEigVal = (EigVals_static[0])[0][0];
   for (ChainIter_.Reset(); !ChainIter_.Done(); ++ChainIter_)
   {
      for (int k = 0; k < INTERNAL_ATOMS; ++k)
      {
         if ((EigVals_static[i])[0][k] < MinEigVal)
         {
            if (fabs((EigVals_static[i])[0][k]) > Tol)
            {
               MinEigVal = (EigVals_static[i])[0][k];
            }
         }
         if ((EigVals_static[i])[0][k] < -Tol)
         {
            NumOfNegEV += 1;
         }
      }
      i = i + 1;
   }

   K[0] = MinEigVal;
   return NumOfNegEV;
}

void SCLDQMultiChainTPP::LongWavelengthModuli(double const& dk, int const& gridsize,
                                              char const* const prefix, ostream& out)
const
{
}

void SCLDQMultiChainTPP::NeighborDistances(int const& cutoff, ostream& out) const
{
   Matrix NeighborDist =
      SCLDChainSum_.NeighborDistances(cutoff, pow(double(10), double(-(out.precision() - 1))));

   int W = out.width();
   int types = (INTERNAL_ATOMS * (INTERNAL_ATOMS + 1)) / 2;
   for (int i = 0; i < cutoff; ++i)
   {
      out << setw(W) << NTemp_ << setw(W) << NeighborDist[i][0];
      for (int j = 0; j < types; ++j)
      {
         out << setw(W / 4) << int(NeighborDist[i][1 + j]);
      }
      out << "\n";
   }
   out << "\n";
}

void SCLDQMultiChainTPP::Print(ostream& out, PrintDetail const& flag,
                               PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions = 0;
   double engy, entropy, heatcapacity;
   str_static.Resize(DOFS);
   TE_static.Resize(DOFS);
   pstiff_static.Resize(DOFS, DOFS);
   Matrix CondEV(1, 1);
   Matrix CondModuli(1, 1);
   TestFunctVals_static.Resize(NumTestFunctions());
   int RankOneConvex;
   Vector K(1);
   int BlochWaveStable;
   double mintestfunct;

   W = out.width();

   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }

   engy = FreeEnergy();

   entropy = Entropy();
   heatcapacity = HeatCapacity();

   str_static = Fstress();
   pstiff_static = Fstiffness();
   TE_static = ThermalExpansion();

   Matrix VibStiff, Stiff;
   VibStiff.Resize(DOFS, DOFS, 0.0);
   Stiff.Resize(DOFS, DOFS, 0.0);

   Stiff = stiffness();
   VibStiff = pstiff_static - Stiff;

   TestFunctions(TestFunctVals_static, LHS);
   mintestfunct = TestFunctVals_static[0];
   for (int i = 0; i < TestFunctVals_static.Dim(); ++i)
   {
      if (TestFunctVals_static[i] < 0.0)
      {
         ++NoNegTestFunctions;
      }
      if (mintestfunct > TestFunctVals_static[i])
      {
         mintestfunct = TestFunctVals_static[i];
      }
   }

   CondModuli = CondensedModuli();

   CondEV = CondModuli;
   RankOneConvex = (CondEV[0][0] > 0) ? 1 : 0;

   K.Resize(1, 0.0);
   //   if (RankOneConvex)
   //   {
   BlochWaveStable = BlochWave(K);
   //   }
   //   else
   //   {
   //      BlochWaveStable = -1;
   //   }


   switch (flag)
   {
      case PrintLong:
         out << "SCLDQMultiChainTPP:" << "\n" << "\n";
         out << "Density_ : " << Density_ << "\n";
         out << "LagrangeCB: " << LagrangeCB_ << "\n";
         out << "RefLattice_ : " << setw(W) << RefLattice_;
         for (int i = 0; i < INTERNAL_ATOMS; ++i)
         {
            out << "Atom_" << i << "          "
                << "Species : " << setw(5) << AtomSpecies_[i]
                << "          Position : " << setw(W) << AtomPositions_[i] << "\n";
         }
         out << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
         for (int i = 0; i < NumberofSpecies_; ++i)
         {
            out << "Atomic Mass " << i << "  : "
                << setw(W) << SpeciesMass_[i] << "\n";
         }
         out << "Tref = " << setw(W) << Tref_ << "\n";
         // << "PhiRef = " << setw(W) << PhiRef_ << "; "
         // << "EntropyRef = " << setw(W) << EntropyRef_ << "; "
         // << "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << "\n";
         out << "Potential Parameters : " << "\n";
         for (int i = 0; i < NumberofSpecies_; ++i)
         {
            for (int j = i; j < NumberofSpecies_; j++)
            {
               out << "[" << i << "][" << j << "] -- "
                   << setw(W) << *SpeciesPotential_[i][j] << "\n";
            }
         }
         out << "Normalization Modulus : " << setw(W) << NormModulus_ << "\n";
         // also send to cout
         if (Echo_)
         {
            cout << "SCLDQMultiChainTPP:" << "\n" << "\n";
            cout << "Density_ : " << Density_ << "\n";
            cout << "LagrangeCB: " << LagrangeCB_ << "\n";
            cout << "RefLattice_ : " << setw(W) << RefLattice_;
            for (int i = 0; i < INTERNAL_ATOMS; ++i)
            {
               cout << "Atom_" << i << "          "
                    << "Species : " << setw(5) << AtomSpecies_[i]
                    << "          Position : " << setw(W) << AtomPositions_[i] << "\n";
            }
            cout << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
            for (int i = 0; i < NumberofSpecies_; ++i)
            {
               cout << "Atomic Mass " << i << "  : "
                    << setw(W) << SpeciesMass_[i] << "\n";
            }
            cout << "Tref = " << setw(W) << Tref_ << "\n";
            // << "PhiRef = " << setw(W) << PhiRef_ << "; "
            // << "EntropyRef = " << setw(W) << EntropyRef_ << "; "
            // << "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << "\n";
            cout << "Potential Parameters : " << "\n";
            for (int i = 0; i < NumberofSpecies_; ++i)
            {
               for (int j = i; j < NumberofSpecies_; j++)
               {
                  cout << "[" << i << "][" << j << "] -- "
                       << setw(W) << *SpeciesPotential_[i][j] << "\n";
               }
            }
            cout << "Normalization Modulus : " << setw(W) << NormModulus_ << "\n";
         }

         FirstPrintLong = 1;
         cout << "FirstConvergedPrintLong: " << FirstConverged << endl;

      // passthrough to short
      case PrintShort:
         out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << "\n"
             << "Lambda (Normalized): " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
             << "Potential Value (Normalized):" << setw(W) << engy << "\n"
             << "Thermal Expansion:" << "\n" << setw(W) << TE_static << "\n\n"
             << "Entropy:" << setw(W) << entropy << "\n"
             << "HeatCapacity:" << setw(W) << heatcapacity << "\n";
         for (int i = 0; i < INTERNAL_ATOMS; ++i)
         {
            out << "BodyForce Value " << i << " (Inf Normalized):"
                << setw(W) << BodyForce_[i] << "\n";
         }
         out << "Stress (Normalized):" << "\n" << setw(W) << str_static << "\n\n"
             << "Stiffness (Normalized):" << setw(W) << pstiff_static
         //             << "Vib. Stiffness (Normalized):" << setw(W) << VibStiff
         //             << "Stat. Stiffness (Normalized):" << setw(W) << Stiff
         << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals_static << "\n"
         << "Bifurcation Info:" << setw(W) << mintestfunct
         << setw(W) << NoNegTestFunctions << "\n"
         << "Condensed Moduli (Normalized):" << setw(W) << CondModuli
         << "CondEV Info:" << setw(W) << CondEV
         << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
         << "BlochWave Stability:" << setw(W) << BlochWaveStable << ", "
         << setw(W) << K << "\n"
         << "EigVals: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(W) << (EigVals_static[i])[0][k];
            }
            out << "\n";
         }
         out << "HEigVals: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(W) << (HEigVals_static[i])[0][k];
            }
            out << "\n";
         }
         out << "EigValsDOF x=1: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(W) << ((EigValsDOF_static[1])[i])[0][k];
            }
            out << "\n";
         }
         out << "HEigValsDOF x=1: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(W) << ((HEigValsDOF_static[1])[i])[0][k];
            }
            out << "\n";
         }
         out << "EigValsDOFDOF x=1, y=1: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(25) << (((EigValsDOFDOF_static[1])[1])[i])[0][k];
            }
            out << "\n";
         }
         out << "HEigValsDOFDOF x=1, y=1: " << "\n";
         for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
         {
            for (int k = 0; k < INTERNAL_ATOMS; ++k)
            {
               out << setw(25) << (((HEigValsDOFDOF_static[1])[1])[i])[0][k];
            }
            out << "\n";
         }

         out << endl;
         // send to cout also
         if (Echo_)
         {
            cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << "\n"
                 << "Lambda (Normalized): " << setw(W) << Lambda_ << "\n"
                 << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                 << "Potential Value (Normalized):" << setw(W) << engy << "\n"
                 << "Thermal Expansion:" << "\n" << setw(W) << TE_static << "\n\n"
                 << "Entropy:" << setw(W) << entropy << "\n"
                 << "HeatCapacity:" << setw(W) << heatcapacity << "\n";
            for (int i = 0; i < INTERNAL_ATOMS; ++i)
            {
               cout << "BodyForce Value " << i << " (Inf Normalized):"
                    << setw(W) << BodyForce_[i] << "\n";
            }
            cout << "Stress (Normalized):" << "\n" << setw(W) << str_static << "\n\n"
                 << "Stiffness (Normalized):" << setw(W) << pstiff_static
                 << "Eigenvalue Info (Translations->1.0):" << "\n"
                 << setw(W) << TestFunctVals_static << "\n"
                 << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n"
                 << "Condensed Moduli (Normalized):" << setw(W) << CondModuli
                 << "CondEV Info:" << setw(W) << CondEV
                 << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
                 << "BlochWave Stability (GridSize=" << GridSize_ << "):"
                 << setw(W) << BlochWaveStable << ", "
                 << setw(W) << K << "\n";
            cout << endl;
         }

         cout << "FirstConvergedPrintShort: " << FirstConverged << endl;
         if ((FirstPrintLong == 2) && (FirstConverged == 1))
         {
            for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
            {
               EigVals_previous_static[i] = EigVals_static[i];
            }

            FirstConverged = 1;
         }
         else if ((FirstPrintLong == 2) && (FirstConverged == 0))
         {
            for (int i = 0; i < (GridSize_ / 2 + 1); ++i)
            {
               EigVals_previous_static[i].Resize(1, INTERNAL_ATOMS, 0.0);
               EigVals_previous_static[i] = EigVals_static[i];
            }

            FirstConverged = 1;
         }
         else if ((FirstPrintLong == 1) && (FirstConverged == 0))
         {
            if (InitialEigVals_available == 0)
            {
               FirstConverged = 0;
               FirstPrintLong = 2;
            }
            else if (InitialEigVals_available == 1)
            {
               FirstConverged = 1;
               FirstPrintLong = 2;
            }
         }
         else
         {
            cout << "gvError FirstConverged FirstPrintLong: Print: " << endl;
            exit(-1);
         }

         break;
   }
   // check for debug mode request
   if (dbg_)
   {
      if (EnterDebugMode())
      {
         cout << setw(W);
         DebugMode();
      }
   }
}

ostream& operator<<(ostream& out, SCLDQMultiChainTPP& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}


// ---------------------- Debug Mode Handler --------------------------


void SCLDQMultiChainTPP::DebugMode()
{
   const char* Commands[] = {
      "INTERNAL_ATOMS",
      "DOFS",
      "InfluenceDist_",
      "NTemp_",
      "DOF_",
      "RefLattice_",
      "Density_",
      "NormModulus_",
      "Lambda_",
      "BodyForce_",
      "AtomicMass_",
      "GridSize_",
      "Potential_",
      "stress",
      "stiffness",
      "CondensedModuli",
      "ReferenceDispersionCurves",
      "ReferenceBlochWave",
      "ReferenceDynamicalStiffness",
      "SetDOF",
      "StressDT",
      "StiffnessDT",
      "SetTemp",
      "SetInfluenceDist",
      "energy",
      "E0",
      "E1",
      "E2",
      "E3",
      "E4",
      "SetGridSize",
      "NeighborDistances",
      "Print-short",
      "Print-long",
      "SetLambda",
      "StressDL",
      "StiffnessDL",
      "FindLatticeSpacing",
      "ConsistencyCheck",
      "dbg_",
      "RefineEqbm",
      "Entropy",
      "SetParameters"
   };
   int NOcommands = 43;

   string response;
   char prompt[] = "Debug > ";
   int W = cout.width();

   cout << setw(0) << prompt;

   getline(cin, response);

   int indx;
   while (response != "q" && response != "quit" && response != "exit")
   {
      indx = 0;
      if (response == Commands[indx++])
      {
         cout << "INTERNAL_ATOMS = " << INTERNAL_ATOMS << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "DOFS = " << DOFS << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "InfluenceDist_ = " << InfluenceDist_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "NTemp_ = " << NTemp_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         for (int i = 0; i < DOFS; ++i)
         {
            cout << "DOF_[" << i << "] = " << DOF_[i] << "\n";
         }
      }
      else if (response == Commands[indx++])
      {
         cout << "RefLattice_= " << setw(W) << RefLattice_;
      }
      else if (response == Commands[indx++])
      {
         cout << "Density_= " << Density_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "NormModulus_= " << NormModulus_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "Lambda_= " << Lambda_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         for (int i = 0; i < INTERNAL_ATOMS; ++i)
         {
            cout << "BodyForce_[" << i << "]= " << setw(W)
                 << BodyForce_[i] << "\n";
         }
      }
      else if (response == Commands[indx++])
      {
         for (int i = 0; i < INTERNAL_ATOMS; ++i)
         {
            cout << "AtomicMass_[" << i << "]= " << setw(W)
                 << AtomicMass_[i] << "\n";
         }
      }
      else if (response == Commands[indx++])
      {
         cout << "GridSize_= " << GridSize_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         for (int i = 0; i < INTERNAL_ATOMS; ++i)
         {
            for (int j = i; j < INTERNAL_ATOMS; ++j)
            {
               cout << "Potential_[" << i << "][" << j << "]= "
                    << setw(W) << Potential_[i][j] << "\n";
            }
         }
      }
      else if (response == Commands[indx++])
      {
         cout << "stress= " << setw(W) << stress();
      }
      else if (response == Commands[indx++])
      {
         cout << "stiffness= " << setw(W) << stiffness();
      }
      else if (response == Commands[indx++])
      {
         cout << "CondensedModuli= " << setw(W) << CondensedModuli();
      }
      else if (response == Commands[indx++])
      {
         Vector K(6, 0.0);
         int NoPTS;
         string prefix;
         int oldEcho_ = Echo_;
         cout << "\tK > ";
         cin >> K;
         cin.sync(); // clear input
         cout << "\tNoPTS > ";
         cin >> NoPTS;
         cin.sync(); // clear input
         cout << "\tprefix > ";
         cin >> prefix;
         cin.sync(); // clear input
         Echo_ = 0;
         cout << "ReferenceDispersionCurves= ";
         ReferenceDispersionCurves(K, NoPTS, prefix.c_str(), cout);
         Echo_ = oldEcho_;
      }
      else if (response == Commands[indx++])
      {
         Vector K(1, 0.0);
         cout << "ReferenceBlochWave= " << ReferenceBlochWave(K) << "\t" << K << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "\tK > ";
         Vector K(1, 0.0);
         cin >> K;
         cin.sync(); // clear input
         cout << "ReferenceDynamicalStiffness= "
              << setw(W) << ReferenceDynamicalStiffness(K, PairPotentials::T0, 0, 0, 0) << "\n";
      }
      else if (response == Commands[indx++])
      {
         Vector DOF(DOFS, 0.0);
         cout << "\tDOF > ";
         cin >> DOF;
         cin.sync(); // clear input
         SetDOF(DOF);
      }
      else if (response == Commands[indx++])
      {
         cout << "StressDT= " << setw(W) << StressDT();
      }
      else if (response == Commands[indx++])
      {
         cout << "StiffnessDT= " << setw(W) << StiffnessDT();
      }
      else if (response == Commands[indx++])
      {
         double Temp;
         cout << "\tTemp > ";
         cin >> Temp;
         cin.sync(); // clear input
         SetTemp(Temp);
      }
      else if (response == Commands[indx++])
      {
         double dist;
         cout << "\tInfluenceDist > ";
         cin >> dist;
         cin.sync(); // clear input
         SetInfluenceDist(dist);
      }
      else if (response == Commands[indx++])
      {
         cout << "energy= " << energy() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "E0= " << setw(W) << E0();
      }
      else if (response == Commands[indx++])
      {
         cout << "E1= " << setw(W) << E1();
      }
      else if (response == Commands[indx++])
      {
         cout << "E2= " << setw(W) << E2();
      }
      else if (response == Commands[indx++])
      {
         cout << "E3= " << setw(W) << E3();
      }
      else if (response == Commands[indx++])
      {
         cout << "E4= " << setw(W) << E4();
      }
      else if (response == Commands[indx++])
      {
         int GridSize;
         cout << "\tGridSize > ";
         cin >> GridSize;
         cin.sync(); // clear input
         SetGridSize(GridSize);
      }
      else if (response == Commands[indx++])
      {
         int oldEcho_ = Echo_;
         int cutoff;
         cout << "\tcutoff > ";
         cin >> cutoff;
         cin.sync(); // clear input
         Echo_ = 0;
         NeighborDistances(cutoff, cout);
         Echo_ = oldEcho_;
      }
      else if (response == Commands[indx++])
      {
         int oldEcho_ = Echo_;
         Echo_ = 0;
         cout << setw(W) << *this;
         Echo_ = oldEcho_;
      }
      else if (response == Commands[indx++])
      {
         int oldEcho_ = Echo_;
         Echo_ = 0;
         cout << setw(W);
         Print(cout, PrintLong);
         Echo_ = oldEcho_;
      }
      else if (response == Commands[indx++])
      {
         double lambda;
         cout << "\tLambda > ";
         cin >> lambda;
         cin.sync(); // clear input
         SetLambda(lambda);
      }
      else if (response == Commands[indx++])
      {
         cout << "StressDL= " << setw(W) << StressDL();
      }
      else if (response == Commands[indx++])
      {
         cout << "StiffnessDL= " << setw(W) << StiffnessDL();
      }
      else if (response == Commands[indx++])
      {
         int iter;
         cout << "\titer > ";
         cin >> iter;
         cin.sync(); // clear input
         FindLatticeSpacing(iter);
      }
      else if (response == Commands[indx++])
      {
         int width;
         int oldEcho = Echo_;
         double epsilon;
         cout << "\tConsistencyEpsilon > ";
         cin >> epsilon;
         cout << "\tWidth > ";
         cin >> width;
         Echo_ = 0;
         ConsistencyCheck(epsilon, width, cout);
         Echo_ = oldEcho;
      }
      else if (response == Commands[indx++])
      {
         cout << "dbg_ = " << dbg_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         double Tol;
         int MaxItr;
         cout << "\tTolerence > ";
         cin >> Tol;
         cout << "\tMaxItr > ";
         cin >> MaxItr;
         RefineEqbm(Tol, MaxItr, &cout);
      }
      else if (response == Commands[indx++])
      {
         cout << "Entropy = " << setw(W) << Entropy() << "\n";
      }
      else if (response == Commands[indx++])
      {
         int no = SpeciesPotential_[0][0]->GetNoParameters();
         double* vals;
         int sze = no * ((NumberofSpecies_ + 1) * (NumberofSpecies_) / 2);
         vals = new double[sze];
         cout << "\tEnter new Parameter values > ";
         for (int i = 0; i < sze; ++i)
         {
            cin >> vals[i];
         }
         SetParameters(vals);
      }
      else if ((response == "?") || (response == "help"))
      {
         cout << setiosflags(ios::left);
         for (int i = 0; i < NOcommands / 2 + NOcommands % 2; ++i)
         {
            cout << "  " << setw(30) << Commands[i];
            if ((i == NOcommands / 2) && !NOcommands % 2)
            {
               cout << "\n";
            }
            else
            {
               cout << setw(30) << Commands[NOcommands / 2 + i] << "\n";
            }

            if (!((i + 1) % 30))
            {
               cout << "more...." << "\n";
               char ans;
               cin.sync(); // clear input
               ans = kbhitWait();
               if (ans == 'q')
               {
                  break;
               }
            }
         }
         cout << resetiosflags(ios::left) << "\n";
      }
      else if ((response == "\n") || (response == ""))
      {
      }
      else
      {
         cout << "!--- Error - Unknown command ---!" << "\n" << "\n";
      }

      cout << "\n" << prompt;
      getline(cin, response);
   }
}


void SCLDQMultiChainTPP::RefineEqbm(double const& Tol, int const& MaxItr,
                                    ostream* const out)
{
   Vector dx(DOFS, 0.0);
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

      SetDOF(DOF_ - dx);

      Stress = E1();

      if (out != 0)
      {
         *out << setw(20) << Stress;

         *out << itr << "\tdx " << dx.Norm() << "\tstress " << Stress.Norm() << "\n";
      }
   }
}
