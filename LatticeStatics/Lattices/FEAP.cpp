#include <fstream>
#include "FEAP.h"

using namespace std;

extern "C" void bfbfeap_main_(char const* const ffin, int* ffinlen,
                              int* const numDOFperNode, int* const numSpcDIM,
                              int* const numNodeInMesh, int* const numElemInMesh,
                              int* const numNodesPerEl, int* const numEqns);
extern "C" void bfbfeap_get_eqn_id_(int* bfb_id);
extern "C" void bfbfeap_get_elem_conn_(int* bfb_ix);
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
   delete[] elmConn_;
   delete[] bcID_;

   config_out_.close();
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
   if (Input.ParameterOK(Hash, "Nbn"))
     {
       nbn_ = Input.getInt(Hash, "Nbn");
     }
   else
     {
       nbn_ = Input.useInt(0, Hash, "Nbn");  // Default Value
     }

   cout << "nbn_ = " << nbn_ << endl;

   // Initialize FEAP, send input file name and get ndf, ndm, etc. back.
   int ffinlen = strlen(ffin_);
   bfbfeap_main_(ffin_, &ffinlen, &ndf_, &ndm_, &numnp_, &nel_, &nen1_, &neq_);

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

   // get and store element connectivity
   elmConn_ = new int[nen1_*nel_];
   bfbfeap_get_elem_conn_(elmConn_);

   // set DOFS_

   DOFS_ = ndf_ * ( numnp_ - nbn_ ) + ndm_ * ndm_ + nbn_ / 2 ;
   DOFS_F_ = ndf_ * numnp_;


   // set DOF_ to initial value
   DOF_.Resize(DOFS_, 0.0);
   DOF_F_.Resize(DOFS_F_, 0.0);


   // get and store reference coordinates
   X_.Resize(ndm_ * numnp_);
   bfbfeap_get_nodal_coords_(&(X_[0]));

   // setup mapping matrix between FEAP and BFB

   if (ndm_ == 3)
     cout << "*WARNING* 3D mesh not taken care of for mapping" << endl;
 
   Map_.Resize(DOFS_F_, DOFS_, 0.0);
   for (int k = 0; k < nbn_; ++k)
     {
       if (ndm_ == 2)
	 {
	   Map_[k * ndf_][0]=X_[k * ndm_];
	   Map_[k * ndf_][1]=X_[k * ndm_+1];
	   Map_[k * ndf_+1][2]=X_[k * ndm_];
	   Map_[k * ndf_+1][3]=X_[k * ndm_+1];	
	   
	   if (ndf_ > ndm_) //1 extra dof assumed
	       Map_[k * ndf_ +2][ndm_ * ndm_ + k%(nbn_ / 2)]=1.0;
	 }
       }

   for (int i = 0 ; i < (numnp_ - nbn_)*ndf_; ++i)
     {
       Map_[nbn_ * ndf_ + i][ndm_ * ndm_ + nbn_ / 2 + i]=1.0;
     }

   // Setup 1st PiolaKirchoff load tensor

   if (Input.ParameterOK(Hash,"FirstPKstress"))
   {
      Load_.Resize(2,2,0.0);
      Input.getMatrix(Load_,Hash,"FirstPKstress");
   }
   else
   {
      cout << "*Error: FEAP::FirstPKstress not found in input file.\n";
      exit(-1);
   }


   // setup remaining variables
   E1CachedValue_.Resize(DOFS_,0.0);
   E1CachedValue_F_.Resize(DOFS_F_,0.0);

   E1DLoadCachedValue_.Resize(DOFS_,0.0);

   E2CachedValue_.Resize(DOFS_, DOFS_);
   E2CachedValue_F_.Resize(DOFS_F_, DOFS_F_);

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

   // open config_out file
   string flnm(Input.LastInputFileName());
   int pos = flnm.find(".bfb");
   if (string::npos != pos)
   {
      flnm.erase(pos,flnm.length());
   }
   flnm.append(".gpl");
   config_out_.open(flnm.c_str() ,ios::out);
   config_count_ = 0;
}

void FEAP::UpdateValues(UpdateFlag flag) const
{
   // Update FEAP solution vector
   DOF_F_ = Map_ * DOF_;
   bfbfeap_set_nodal_solution_(&(DOF_F_[0]));
   
   // Set disp gradient for energy
   U_.Resize(2,2,0.0);
   U_[0][0] = DOF_[0];
   U_[0][1] = DOF_[1];
   U_[1][0] = DOF_[2];
   U_[1][1] = DOF_[3];


   if (NoStiffness == flag)
   {
      bfbfeap_call_ener_(); // needs a "TPLOt" and "ENER" command in FEAP input file to work
      bfbfeap_get_potential_energy_(&(E0CachedValue_));
      E0CachedValue_ -= Lambda_ * (Load_ * U_).Trace();

      bfbfeap_call_form_();
      bfbfeap_get_reduced_residual_(&(E1CachedValue_F_[0]));

      // FEAP returns -E1, so fix it.
      for (int i = 0; i < E1CachedValue_F_.Dim(); ++i)
      {
         E1CachedValue_F_[i] = -E1CachedValue_F_[i];
      }
      E1CachedValue_ = Map_.Transpose() * E1CachedValue_F_;

      E1CachedValue_[0] -= Lambda_ * Load_[0][0];
      E1CachedValue_[1] -= Lambda_ * Load_[1][0];
      E1CachedValue_[2] -= Lambda_ * Load_[0][1];
      E1CachedValue_[3] -= Lambda_ * Load_[1][1];   

      Cached_[0] = 1;
      Cached_[1] = 1;
      EvaluationCount_[0]++;
   }
   else if (NeedStiffness == flag)
   {
      bfbfeap_call_ener_();
      bfbfeap_get_potential_energy_(&(E0CachedValue_));
      E0CachedValue_ -= Lambda_ * (Load_ * U_).Trace();

      bfbfeap_call_form_();
      bfbfeap_get_reduced_residual_(&(E1CachedValue_F_[0]));
      // FEAP returns -E1, so fix it.
      for (int i = 0; i < E1CachedValue_F_.Dim(); ++i)
      {
         E1CachedValue_F_[i] = -E1CachedValue_F_[i];
      }
      E1CachedValue_ = Map_.Transpose() * E1CachedValue_F_;

      E1CachedValue_[0] -= Lambda_ * Load_[0][0];
      E1CachedValue_[1] -= Lambda_ * Load_[1][0];
      E1CachedValue_[2] -= Lambda_ * Load_[0][1];
      E1CachedValue_[3] -= Lambda_ * Load_[1][1];

      //E1DLoad 

      E1DLoadCachedValue_[0] = -Load_[0][0];
      E1DLoadCachedValue_[1] = -Load_[1][0];
      E1DLoadCachedValue_[2] = -Load_[0][1];
      E1DLoadCachedValue_[3] = -Load_[1][1];

      //E2

      bfbfeap_call_tang_();
      bfbfeap_get_reduced_tang_(&(E2CachedValue_F_[0][0]));
      
      E2CachedValue_ = Map_.Transpose() * E2CachedValue_F_ * Map_;

      
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

void FEAP::print_gpl_config(fstream& out) const
{
   out << "# Index num: " << config_count_ << endl;
   out << "########## configuration for Lambda = " << Lambda_
       << " ##########" << endl;
   for (int i=0;i<nel_;++i)
   {
      // node numbers for element i (-1 for zero-based values)
      int n1 = elmConn_[i*nen1_+0] - 1;
      int n2 = elmConn_[i*nen1_+1] - 1;
      out << setw(Width_) << X_[n1*ndm_+0] + DOF_F_[n1*ndf_+0]
          << setw(Width_) << X_[n1*ndm_+1] + DOF_F_[n1*ndf_+1]
          << endl;
      out << setw(Width_) << X_[n2*ndm_+0] + DOF_F_[n2*ndf_+0]
          << setw(Width_) << X_[n2*ndm_+1] + DOF_F_[n2*ndf_+1]
          << endl << endl;
   }
   out << endl;

   ++config_count_;
}

void FEAP::Print(ostream& out, PrintDetail const& flag,
                 PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions = 0;
   double engy;
   double E1norm;
   double mintestfunct[3];
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
   for (int i=0;i<3;++i) mintestfunct[i] = TestFunctVals[0];
   // check only the EigenValTFs
   for (int i = 0; i < DOFS_; ++i)
   {
      if ((UseEigenValTFs() == 1) && (TestFunctVals[i] < 0.0))
      {
         ++NoNegTestFunctions;
      }
      if (mintestfunct[0] > TestFunctVals[i])
      {
         mintestfunct[2] = mintestfunct[1];
         mintestfunct[1] = mintestfunct[0];
         mintestfunct[0] = TestFunctVals[i];
      }
   }

   print_gpl_config(config_out_);

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



         out << "Bifurcation Info: " << setw(W) << mintestfunct[0]
             << setw(W) << mintestfunct[1]
             << setw(W) << mintestfunct[2]
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "DOF: " << setw(W) << DOF_ << "\n"
                 << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
                 << "Potential Value: " << setw(W) << engy << "\n"
                 << "Force Norm: " << setw(W) << E1norm << "\n";



            cout << "Bifurcation Info: " << setw(W) << mintestfunct[0]
                 << setw(W) << mintestfunct[1]
                 << setw(W) << mintestfunct[2]
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
