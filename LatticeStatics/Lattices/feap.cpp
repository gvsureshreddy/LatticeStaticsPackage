#include <fstream>
#include "feap.h"

using namespace std;

extern "C" void bfbfeap_main_(char const* const ffin, int* ffinlen,
                              int* const numDOFperNode, int* const numSpcDIM,
                              int* const numNodeInMesh, int* const numEqns);
extern "C" void bfbfeap_get_eqn_id_(int* bfb_id);
extern "C" void bfbfeap_get_eqn_bc_(int* bfb_bc);
extern "C" void bfbfeap_get_nodal_solution_(double* bfb_u);
extern "C" void bfbfeap_set_nodal_solution_(double* bfb_u);
extern "C" void bfbfeap_get_nodal_coords_(double* bfb_x);
extern "C" void bfbfeap_get_reduced_residual_(double* bfb_rd);
extern "C" void bfbfeap_get_reduced_tang_(double* bfb_tang);
extern "C" void bfbfeap_call_form_();
extern "C" void bfbfeap_call_tang_();

extern "C" void plstop_();
#define FORTRANSTRINGLEN 129

feap::~feap()
{
   char ret[10];
   cout << "calling plstop: press enter";
   cin >> ret;
   plstop_();
   cout << "feap Function Calls:\n"
        << "\tEvaluation - w/o stiffness - " << EvaluationCount_[0] << "\n"
        << "\tEvaluation - w   stiffness - " << EvaluationCount_[1] << "\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

feap::feap(PerlInput const& Input, int const& Echo, int const& Width) :
   Lattice(Input, Echo),
   Lambda_(0.0),
   Width_(Width)
{
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash, "feap");

   if (Input.ParameterOK(Hash, "Tolerance"))
   {
      Tolerance_ = Input.getDouble(Hash, "Tolerance");
   }
   else
   {
      Tolerance_ = Input.useDouble(1.0e-6, Hash, "Tolerance");  // Default Value
   }

   char const* ffin = "Ipatch_4el";

   int ndf;
   int ndm;
   int numnp;
   int neq;
   int ffinlen = strlen(ffin);
   bfbfeap_main_(ffin, &ffinlen, &ndf, &ndm, &numnp, &neq);

   int* id = new int[ndf*numnp];
   bfbfeap_get_eqn_id_(id);
   for (int i = 0; i < ndf*numnp; ++i)
   {
      cout << setw(Width) << id[i];
   }
   cout << endl;
   delete [] id;

   int* bc = new int[ndf*numnp];
   bfbfeap_get_eqn_bc_(bc);
   for (int i = 0; i < ndf*numnp; ++i)
   {
      cout << setw(Width) << bc[i];
   }
   cout << endl;
   delete [] bc;
   
   DOFS_ = ndf*numnp;
   cout << "DOFS_ is " << DOFS_ <<endl;

   DOF_.Resize(DOFS_, 1.0);
   cout << "DOF is " << setw(Width) << DOF_ << endl;
   bfbfeap_get_nodal_solution_(&(DOF_[0]));
   cout << "DOF is " << setw(Width) << DOF_ << endl;

   Vector feap_residual(neq, 1.0);
   cout << "feap_residual is   " << setw(Width) << feap_residual << endl;
   bfbfeap_get_reduced_residual_(&(feap_residual[0]));
   cout << "feap_residual is   " << setw(Width) << feap_residual << endl;
   bfbfeap_call_form_();
   cout << "called form\n";
   bfbfeap_get_reduced_residual_(&(feap_residual[0]));
   cout << "feap_residual is   " << setw(Width) << feap_residual << endl;
   
   Matrix TANG(neq,neq,0.0);
   bfbfeap_get_reduced_tang_(&(TANG[0][0]));
   cout << "TANG is" << endl << setw(Width) << TANG;
   bfbfeap_call_tang_();
   cout << "called tang\n";
   bfbfeap_get_reduced_tang_(&(TANG[0][0]));
   cout << "TANG is" << endl << setw(Width) << TANG;

   Vector X(DOFS_);
   bfbfeap_get_nodal_coords_(&(X[0]));
   cout << "X is " << setw(Width) << X << endl;
   for (int i=0;i<numnp;++i)
   {
      cout << setw(5) << i;
      for (int j=0;j<ndm;++j)
         cout << setw(Width) << X[ndm*i + j];
      cout << endl;
   }

   DOF_[2*0] = -1.0;
   DOF_[2*3] = -1.0;
   DOF_[2*6] = -1.0;
   DOF_[2*2] =  1.0;
   DOF_[2*5] =  1.0;
   DOF_[2*8] =  1.0;

   bfbfeap_set_nodal_solution_(&(DOF_[0]));

   Vector Chk(ndf*numnp,-1.0);
   bfbfeap_get_nodal_solution_(&(Chk[0]));
   cout << "DOF is " << setw(Width) << DOF_ << endl;
   cout << "Chk is " << setw(Width) << Chk << endl;

   bfbfeap_call_form_();
   cout << "called form\n";
   bfbfeap_get_reduced_residual_(&(feap_residual[0]));
   cout << "feap_residual is   " << setw(Width) << feap_residual << endl;

   cout << "check this value   " << setw(Width) << -TANG*DOF_ << endl;

//   E1CachedValue_.Resize(DOFS_);
//   E1DLoadCachedValue_.Resize(DOFS_);
//   E2CachedValue_.Resize(DOFS_, DOFS_);
//   stiffdl_static.Resize(DOFS_, DOFS_);
//   EmptyV_.Resize(DOFS_, 0.0);
//   EmptyM_.Resize(DOFS_, DOFS_, 0.0);
//
//   LoadParameter_ = Load;
//   for (int i = 0; i < cachesize; ++i)
//   {
//      Cached_[i] = 0;
//      CallCount_[i] = 0;
//   }
//   EvaluationCount_[0] = 0;
//   EvaluationCount_[1] = 0;
}

void feap::UpdateValues(UpdateFlag flag) const
{
   if (NoStiffness == flag)
   {
      int mode = 0;
      //qcbfb_energy_(mode, DOFS_, &(DOF_[0]), Lambda_, E0CachedValue_, &(E1CachedValue_[0]), 0, 0);
      Cached_[0] = 1;
      Cached_[1] = 1;
      EvaluationCount_[0]++;
   }
   else if (NeedStiffness == flag)
   {
      int mode = 1;
      //qcbfb_energy_(mode, DOFS_, &(DOF_[0]), Lambda_, E0CachedValue_, &(E1CachedValue_[0]),
      //              &(E2CachedValue_[0][0]), &(E1DLoadCachedValue_[0]));
      Cached_[0] = 1;
      Cached_[1] = 1;
      Cached_[2] = 1;
      Cached_[3] = 1;
      EvaluationCount_[1]++;
   }
   else
   {
      cerr << "Error in feap::UpdateValues(), unknown UpdateFlag.\n";
      exit(-45);
   }
}

double feap::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues(NoStiffness);
   }
   CallCount_[0]++;

   return E0CachedValue_;
}

Vector const& feap::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
   }
   CallCount_[1]++;

   return E1CachedValue_;
}

Vector const& feap::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[2]++;

   return E1DLoadCachedValue_;
}

Matrix const& feap::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[3]++;

   return E2CachedValue_;
}

Matrix const& feap::StiffnessDL() const
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

void feap::Print(ostream& out, PrintDetail const& flag,
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
         out << "feap:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "feap:" << "\n" << "\n";
         }
      // passthrough to short
      case PrintShort:
         out << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
             << "Potential Value: " << setw(W) << engy << "\n"
             << "Force Norm: " << setw(W) << E1norm << "\n";

         out << "Bifurcation Info: " << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
                 << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
                 << "Potential Value: " << setw(W) << engy << "\n"
                 << "Force Norm: " << setw(W) << E1norm << "\n";

            cout << "Bifurcation Info: " << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out, feap& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}
