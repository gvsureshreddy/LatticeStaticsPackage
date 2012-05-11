#include <fstream>
#include "FEAP.h"

using namespace std;

extern "C" void bfbfeap_main_(char const* const ffin, int* ffinlen,
                              int* const numDOFperNode, int* const numSpcDIM,
                              int* const numNodeInMesh, int* const numEqns);
extern "C" void bfbfeap_get_eqn_id_(int* bfb_id);
extern "C" void bfbfeap_get_eqn_bc_(int* bfb_bc);
extern "C" void bfbfeap_get_nodal_solution_(double* bfb_u);
extern "C" void bfbfeap_set_nodal_solution_(double* bfb_u);
extern "C" void bfbfeap_get_nodal_coords_(double* bfb_x);
extern "C" void bfbfeap_get_potential_energy_(double* bfb_epl);
extern "C" void bfbfeap_get_reduced_residual_(double* bfb_rd);
extern "C" void bfbfeap_get_reduced_tang_(double* bfb_tang);
extern "C" void bfbfeap_call_ener_();
extern "C" void bfbfeap_call_form_();
extern "C" void bfbfeap_call_tang_();
extern "C" void plstop_();
#define FORTRANSTRINGLEN 129

FEAP::~FEAP()
{
   plstop_(); // Shut down FEAP

   cout << "FEAP Function Calls:\n"
        << "\tEvaluation - w/o stiffness - " << EvaluationCount_[0] << "\n"
        << "\tEvaluation - w   stiffness - " << EvaluationCount_[1] << "\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";

   delete[] eqnID_;
   delete[] bcID_;
}

FEAP::FEAP(PerlInput const& Input, int const& Echo, int const& Width) :
   Lattice(Input, Echo),
   Lambda_(0.0),
   Width_(Width)
{
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash, "FEAP");

   ffin_ = Input.getString(Hash, "FEAPInputFileName");
   if (Input.ParameterOK(Hash, "Tolerance"))
   {
      Tolerance_ = Input.getDouble(Hash, "Tolerance");
   }
   else
   {
      Tolerance_ = Input.useDouble(1.0e-6, Hash, "Tolerance");  // Default Value
   }

   // Initialize FEAP, send input file name and get ndf, ndm, etc. back.
   int ffinlen = strlen(ffin_);
   bfbfeap_main_(ffin_, &ffinlen, &ndf_, &ndm_, &numnp_, &neq_);

   // get and store equation id's from FEAP
   eqnID_ = new int[ndf_ * numnp_];
   bfbfeap_get_eqn_id_(eqnID_);

   // get and store boundary condition id's from FEAP
   bcID_ = new int[ndf_ * numnp_];
   bfbfeap_get_eqn_bc_(bcID_);
   for (int i = 0; i < ndf_ * numnp_; ++i)
   {
      if (bcID_[i] != 0)
      {
         cerr << "*WARNING* FEAP::FEAP() -- Found displacement BC, but expecting none...\n";
      }
   }

   // set DOFS_
   DOFS_ = ndf_ * ( numnp_ - nbn_ ) + ndm_ * ndm_ + nbn_ / 2 ;
   DOFS_V_ = ndf_ * numnp_ ;
   // set DOF_ to initial value
   DOF_.Resize(DOFS_, 1.0);
   bfbfeap_get_nodal_solution_(&(DOF_[0]));

   // get and store reference coordinates
   X_.Resize(DOFS_);
   bfbfeap_get_nodal_coords_(&(X_[0]));

   // setup remaining variables
   E1CachedValue_.Resize(DOFS_);
   E1DLoadCachedValue_.Resize(DOFS_);
   E2CachedValue_.Resize(DOFS_, DOFS_);
   stiffdl_static.Resize(DOFS_, DOFS_);
   EmptyV_.Resize(DOFS_, 0.0);
   EmptyM_.Resize(DOFS_, DOFS_, 0.0);

   // set loadparameter to load (not temperature)
   LoadParameter_ = Load;

   for (int i = 0; i < cachesize; ++i)
   {
      // indicate that none of the current values are valid
      Cached_[i] = 0;
      // initialize call count (number of times member function is called) to zero
      CallCount_[i] = 0;
   }

   // initialize evaluation count (number of times UpdateValues is called) to zero
   EvaluationCount_[0] = 0;
   EvaluationCount_[1] = 0;
}

void FEAP::UpdateValues(UpdateFlag flag) const
{
   // Update FEAP solution vector
   bfbfeap_set_nodal_solution_(&(DOF_[0]));
   double eps = 2.0;
   double DispSum [ndm_];
   for (int i = 0; i < ndm_; ++i)
   {
      DispSum[i]=0.0;
      for (int j=0; j < numnp_; ++j)
      {
	 DispSum[i] += DOF_[j*ndf_+i];     
      }
    }	

   if (NoStiffness == flag)
   {
      int mode = 0;
      bfbfeap_call_ener_(); // needs a "TPLOt" and "ENER" command in FEAP input file to work
      bfbfeap_get_potential_energy_(&(E0CachedValue_));
      // Phantom energy term for E0
      for (int i=0; i < ndm_; ++i)
      {
         E0CachedValue_ += 1.0/eps*DispSum[i]*DispSum[i];
      }
      bfbfeap_call_form_();
      bfbfeap_get_reduced_residual_(&(E1CachedValue_[0]));
      // FEAP returns -E1, so fix it.
      for (int i = 0; i < E1CachedValue_.Dim(); ++i)
      {
         E1CachedValue_[i] = -E1CachedValue_[i];
      }
      // Phantom energy term for E1
      for (int i = 0; i < ndm_; ++i)
      {
         for (int j = 0; j < numnp_; ++j)
         {
           E1CachedValue_[j* ndf_ + i] += 2.0/eps*DispSum[i];
         }
      }
      Cached_[0] = 1;
      Cached_[1] = 1;
      EvaluationCount_[0]++;
   }
   else if (NeedStiffness == flag)
   {
      int mode = 1;
      bfbfeap_call_ener_();
      bfbfeap_get_potential_energy_(&(E0CachedValue_));
      // Phantom energy term for E0
      for (int i=0; i < ndm_; ++i)
      {
         E0CachedValue_ += 1.0/eps*DispSum[i]*DispSum[i];
      }

      bfbfeap_call_form_();
      bfbfeap_get_reduced_residual_(&(E1CachedValue_[0]));
      // FEAP returns -E1, so fix it.
      for (int i = 0; i < E1CachedValue_.Dim(); ++i)
      {
         E1CachedValue_[i] = -E1CachedValue_[i];
      }
      // Phantom energy term for E1
      for (int i = 0; i < ndm_; ++i)
      {
         for (int j = 0; j < numnp_; ++j)
         {
             E1CachedValue_[j * ndf_ + i] += 2.0/eps*DispSum[i];
         }
      }

      bfbfeap_call_tang_();
      bfbfeap_get_reduced_tang_(&(E2CachedValue_[0][0]));
      
      // Phantom energy term for E2
      for (int i=0; i < DOFS_; ++i)
      {
         for (int j=0; j<numnp_; ++j)
         {
            if (!((i%ndf_)>(ndm_ - 1)))
            E2CachedValue_[i][ndf_*j+(i%ndf_)] += 2.0/eps;
         }
      }
      // Need an E1DLoad
      Cached_[0] = 1;
      Cached_[1] = 1;
      Cached_[2] = 1;
      Cached_[3] = 1;
      EvaluationCount_[1]++;
   }
   else
   {
      cerr << "Error in FEAP::UpdateValues(), unknown UpdateFlag.\n";
      exit(-45);
   }
}

double FEAP::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues(NoStiffness);
   }
   CallCount_[0]++;

   return E0CachedValue_;
}

Vector const& FEAP::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
   }
   CallCount_[1]++;

   return E1CachedValue_;
}

Vector const& FEAP::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[2]++;

   return E1DLoadCachedValue_;
}

Matrix const& FEAP::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[3]++;

   return E2CachedValue_;
}

Matrix const& FEAP::StiffnessDL() const
{
   double load = Lambda_;

   Lambda_ = load + 10.0 * Tolerance_; for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   stiffdl_static = E2();
   Lambda_ = load - 10.0 * Tolerance_; for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   stiffdl_static -= E2();
   stiffdl_static /= 2.0 * Tolerance_;

   Lambda_ = load; for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   return stiffdl_static;
}

void FEAP::Print(ostream& out, PrintDetail const& flag,
                 PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions = 0;
   double engy;
   double E1norm;
   double mintestfunct;
   Vector TestFunctVals(NumTestFunctions());

   W = out.width();

   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }

   engy = E0();
   E1norm = E1().Norm();

   TestFunctions(TestFunctVals, LHS);
   mintestfunct = TestFunctVals[0];
   // check only the EigenValTFs
   for (int i = 0; i < DOFS_; ++i)
   {
      if ((UseEigenValTFs() == 1) && (TestFunctVals[i] < 0.0))
      {
         ++NoNegTestFunctions;
      }
      if (mintestfunct > TestFunctVals[i])
      {
         mintestfunct = TestFunctVals[i];
      }
   }

   switch (flag)
   {
      case PrintLong:
         out << "FEAP:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "FEAP:" << "\n" << "\n";
         }
      // passthrough to short
      case PrintShort:
         out << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "DOF: " << setw(W) << DOF_ << "\n"
             << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
             << "Potential Value: " << setw(W) << engy << "\n"
             << "Force Norm: " << setw(W) << E1norm << "\n";

         out << "Bifurcation Info: " << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "DOF: " << setw(W) << DOF_ << "\n"
                 << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
                 << "Potential Value: " << setw(W) << engy << "\n"
                 << "Force Norm: " << setw(W) << E1norm << "\n";

            cout << "Bifurcation Info: " << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out, FEAP& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}
