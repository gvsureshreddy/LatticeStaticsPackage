#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <ostream>
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

   delete[] BoundNodes_;
   delete[] PeriodicNodes_;
   delete[] InnerNodes_;
   delete[] Map_[0];
   delete[] Map_;
   delete[] FMap_[0];
   delete[] FMap_;
   if (0 != N_) delete[] N_[0];
   delete[] N_;

   delete[] eqnID_;
   delete[] elmConn_;
   delete[] bcID_;

   config_out_.close();
   plot_out_.close();
   bloch_wave_out_.close();
   critical_eig_out_.close();
}

FEAP::FEAP(PerlInput const& Input, int const& Echo, int const& Width) :
   Lattice(Input, Echo),
   Lambda_(0.0),
   Width_(Width),
   BoundNodes_(0),
   PeriodicNodes_(0),
   InnerNodes_(0),
   Map_(0),
   FMap_(0),
   N_(0)
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
      Tolerance_ = Input.useDouble(1.0e-10, Hash, "Tolerance");  // Default Value
   }
   if (Input.ParameterOK(Hash, "Nbn"))
     {
       nbn_ = Input.getInt(Hash, "Nbn");
     }
   else
     {
       nbn_ = Input.useInt(0, Hash, "Nbn");  // Default Value
     }
   if (Input.ParameterOK(Hash, "Nuc"))
     {
       nuc_ = Input.getInt(Hash, "Nuc");
     }
   else
     {
       nuc_ = Input.useInt(1, Hash, "Nuc");  // Default Value
     }
   if (Input.ParameterOK(Hash, "HexagonDiagonalLength"))
     {
       HexSize_ = Input.getDouble(Hash, "HexagonDiagonalLength");
     }
   else
     {
       HexSize_ = Input.useDouble(1.0, Hash, "HexagonDiagonalLength");  // Default Value
     }


   if (Input.ParameterOK(Hash, "PhantomEnergyEpsilon"))
     {
       eps_ = Input.getDouble(Hash, "PhantomEnergyEpsilon");
     }
   else
     {
       eps_ = Input.useDouble(1.0, Hash, "PhantomEnergyEpsilon");  // Default Value
     }

   if (Input.ParameterOK(Hash, "BoundNodes"))
   {
     BoundNodes_=new int[nbn_/2];
     Input.getIntVector(BoundNodes_,nbn_/2,Hash,"BoundNodes");
     for (int i = 0; i < nbn_/2; ++i)
     {
       BoundNodes_[i]=BoundNodes_[i]-1;
//        cout << "BoundNodes_[" << i << "]=" <<BoundNodes_[i] << "\n";
     }
   }
   else
   {
      cout << "*Error: FEAP::lists of boundary nodes not found in input file.\n";
      exit(-1);
   }

   if (Input.ParameterOK(Hash, "PeriodicNodes"))
   {
     PeriodicNodes_=new int[nbn_/2];
     Input.getIntVector(PeriodicNodes_,nbn_/2,Hash,"PeriodicNodes");
     for (int i = 0; i < nbn_/2; ++i)
     {
       PeriodicNodes_[i]=PeriodicNodes_[i]-1;
//        cout << "PeriodicNodes_[" << i << "]=" <<PeriodicNodes_[i] << "\n";
     }
   }
   else
   {
      cout << "*Error: FEAP::lists of periodic nodes not found in input file.\n";
      exit(-1);
   }

   // Setup loading
   if (Input.ParameterOK(Hash,"Loading"))
   {
     PerlInput::HashStruct LoadHash = Input.getHash(Hash,"Loading");
     if (Input.ParameterOK(LoadHash,"Type"))
     {
       const char* LoadT = Input.getString(LoadHash,"Type");
       if (!strcmp("DeadLoad",LoadT))
       {
         LoadingType_ = DEAD_LOAD;
         Load_.Resize(2,2,0.0);
         Input.getMatrix(Load_,LoadHash,"BioStress");
       }
       else if (!strcmp("PressureLoad",LoadT))
       {
         LoadingType_ = PRESSURE_LOAD;
       }
       else if (!strcmp("DisplacementControl",LoadT))
       {
         LoadingType_ = DISPLACEMENT_CONTROL;
         StretchRatio_ = Input.getDouble(LoadHash,"StretchRatio");
       }
       else
       {
         cerr << "Error (FEAP()): Unknown value for Loading:Type" << "\n";
         exit(-3);
       }
     }
     else
     {
       cerr << "Error (FEAP()): Missing Loading{Type}" << "\n";
       exit(-3);
     }
   }
   else
   {
     cerr << "Error (FEAP()): Missing Loading" << "\n";
   }

   // Initialize FEAP, send input file name and get ndf, ndm, etc. back.
   int ffinlen = strlen(ffin_);
   bfbfeap_main_(ffin_, &ffinlen, &ndf_, &ndm_, &numnp_, &nel_, &nen1_, &neq_);

   KDirection_.Resize(2,0.0);
   PerlInput::HashStruct TFHash = Input.getHash(Hash, "ExtraTestFunctions");
   const char* TFtyp = Input.getString(TFHash, "Type");
   if ((!strcmp("None", TFtyp)) || (!strcmp("none", TFtyp)))
   {
      TFType_ = 0;
      NumExtraTFs_ = 0;
   }
   else if(!strcmp("BlochWaveAnalysis", TFtyp))
   {
      K_.Resize(2,0.0);

      if (Input.ParameterOK(TFHash, "KSpaceResolution"))
      {
        KSpaceResolution_ = Input.getInt(TFHash,"KSpaceResolution");

      }
      else
      {
         KSpaceResolution_ = 6;
      }

      if (Input.ParameterOK(TFHash, "DynamicalStiffnessInfo"))
      {
         N_rows_ = Input.getArrayLength(TFHash,"DynamicalStiffnessInfo");
         N_cols_ = Input.getArrayLength(TFHash,"DynamicalStiffnessInfo",1,-1,-1);

         int* N;
         N = new int[N_rows_*N_cols_];
         Input.getIntMatrix(N,N_rows_,N_cols_, TFHash, "DynamicalStiffnessInfo");
         N_ = new int*[N_rows_];
         N_[0] = new int[N_cols_*N_rows_];
         for (int i = 1; i < N_rows_; ++i)
         {
           N_[i] = (N_[i-1] + N_cols_);
         }
         for (int i = 0; i < N_rows_*N_cols_; ++i)
         {
            N_[i/N_cols_][i%N_cols_] = N[i];
         }
         delete[] N;

      }
      else
      {
         cerr << "*ERROR* Unknown DynamicalStiffnessInfo \n";
         exit(-3);
      }



      if (Input.ParameterOK(TFHash, "AnalysisType"))
      {
          const char* AnalysisType = Input.getString(TFHash, "AnalysisType");
          if (!strcmp("Full", AnalysisType))
          {
             TFType_ = 1;
             NumExtraTFs_ = ndf_*(numnp_-nbn_/2) * pow(KSpaceResolution_+1.0,2);
             if (!(Input.ParameterOK(TFHash, "KSpaceResolution")))
             {
                cout << "*WARNING* KSpaceResolution not specified. Using default value: 6 \n";
             }

          }
          else if (!strcmp("KDirection",AnalysisType))
          {
             TFType_ = 2;
             NumExtraTFs_ = ndf_*(numnp_-nbn_/2)*(KSpaceResolution_ + 1);
             KDirection_.Resize(2,0.0);
             KRange_.Resize(2,0.0);

             if (!(Input.ParameterOK(TFHash, "KSpaceResolution")))
             {
                cout << "*WARNING* KSpaceResolution not specified. Using default value: KSpaceResolution = 6 \n";
             }

             if (Input.ParameterOK(TFHash, "KDirection"))
             {
                Input.getVector(KDirection_,TFHash, "KDirection");
             }
             else
             {
                cerr << "*ERROR* KDirection not specified. \n";
                exit(-3);
             }

             if (Input.ParameterOK(TFHash, "KRange"))
             {
                Input.getVector(KRange_,TFHash,"KRange");
             }
             else
             {
                KRange_[0]=0.0;
                KRange_[1]=0.5;
                cout << "*WARNING* KRange not specified. Using default value: KRange = [0.0,0.5] \n";
             }

         }
         else if (!strcmp("KVectors",AnalysisType))
         {
            TFType_ = 3;
            NumKVectors_ = Input.getArrayLength(TFHash,"KVectors");
            KVectorMatrix_.Resize(NumKVectors_,2,0.0);
            Input.getMatrix(KVectorMatrix_, TFHash, "KVectors");
            NumExtraTFs_ = NumKVectors_*ndf_*(numnp_-nbn_/2);

         }
         else
         {
            cerr << "Error FEAP(): Unknown Bloch wave analysis type \n";
            exit(-3);
         }
      }
   }
   else if(!strcmp("LoadingParameter", TFtyp))
   {
      TFType_ = 4;
      NumExtraTFs_ = 1;

      TFLoad_=  Input.getDouble(TFHash,"LoadingParameter");
   }
   else
   {
      cerr << "Error (FEAP()): Unknown TestFunctions{Type}" << "\n";
      exit(-3);
   }


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
   DOFS_ = ndf_ * ( numnp_ - nbn_ / 2); // just the internal dofs
   switch (LoadingType_) // add any uniform deformation dofs
   {
     case DEAD_LOAD:
     case PRESSURE_LOAD:
       DOFS_ += ndm_ * (ndm_ + 1) / 2;
       break;
     case DISPLACEMENT_CONTROL:
       // nothing to add
       break;
   }
   DOFS_F_ = ndf_ * numnp_;

   // new CellArea_
   CellArea_ = HexSize_ * HexSize_ * 3.0 * sqrt(3.0) / 8.0;
   // old CellArea_
//    CellArea_ = HexSize_ * sqrt(3.0) / 4.0;

   cout << "ndf = " << ndf_ << " numnp = " << numnp_ << " ndm = " << ndm_ << " nbn = " << nbn_ << "\n";
   cout << "DOFS = " << DOFS_ << " DOFS_F = " << DOFS_F_ << "\n";

   // set DOF_ to initial value
   U_.Resize(2,2,0.0);
   F_.Resize(2,2,0.0);
   Vector initDOF(DOFS_,0.0);
   initDOF[0] = initDOF[1] = 1.0;
   DOF_.Resize(DOFS_, 0.0);
   SetDOF(initDOF);
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     SetLambda(1.0);
   }
   else
   {
     SetLambda(0.0);
   }
   DOF_F_.Resize(DOFS_F_, 0.0);

   // get and store reference coordinates
   X_.Resize(ndm_ * numnp_);
   bfbfeap_get_nodal_coords_(&(X_[0]));

   cout << "\n";
   X_F_.Resize(DOFS_F_,0.0);
   for (int i = 0; i < numnp_; ++i)
   {
      X_F_[ndf_ * i] = X_[ndm_ * i];
      X_F_[ndf_ * i + 1] = X_[ndm_ * i + 1];
   }


   // setup remaining variables
   E1CachedValue_.Resize(DOFS_,0.0);
   W1CachedValue_.Resize(DOFS_+1,0.0);
   DispE1CachedValue_.Resize(3,0.0);
   E1CachedValue_F_.Resize(DOFS_F_,0.0);

   E1DLoadCachedValue_.Resize(DOFS_,0.0);

   E2CachedValue_.Resize(DOFS_, DOFS_,0.0);
   W2CachedValue_.Resize(DOFS_+1, DOFS_+1,0.0);
   E2CachedValue_F_.Resize(DOFS_F_, DOFS_F_,0.0);

   stiffdl_static.Resize(DOFS_, DOFS_,0.0);
   EmptyV_.Resize(DOFS_, 0.0);
   EmptyM_.Resize(DOFS_, DOFS_, 0.0);

   Jacobian_.Resize(DOFS_F_, DOFS_,0.0);
   FJacobian_.Resize(DOFS_F_, DOFS_+1,0.0);
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     DispJacobian_.Resize(DOFS_F_, 3,0.0);
   }

   int dim = ndf_*(numnp_ - nbn_/2);
   Dk_.Resize(dim,dim,0.0);


   // Setup Map_ matrix (2D representation of 3D array d2(V)/d(U)2
   // The first pair (i,j) correspond to d2(Vk)/d(Ui)d(Uj) = 1.0, the second to 1/sqrt(2)
   // the % symbol in a%b stands for the rest of the euclidian division of a by b
   // mapping from FEAP dof to Lattice ones

   // This term doesn't contribute anything for displacement control, so just leave it be
   Map_ = new int*[DOFS_F_];
   Map_[0] = new int[4*DOFS_F_];
   FMap_ = new int*[DOFS_F_];
   FMap_[0] = new int[4*DOFS_F_];
   for (int i = 1; i < DOFS_F_; ++i)
   {
     Map_[i] = Map_[i-1] + 4;
     FMap_[i] = FMap_[i-1] + 4;
   }

   for (int i = 0; i < nbn_/2; ++i)
   {
     MapHelper(i, 3+(i%(nbn_/2))*ndf_, BoundNodes_);
   }

   for (int i = 0; i < nbn_/2; ++i)
   {
     MapHelper(i, 3+(i%(nbn_/2))*ndf_, PeriodicNodes_);
   }

   int l=-1;
   InnerNodes_ = new int[numnp_-nbn_];
   for (int i=0; i < numnp_; ++i)
   {
     l=l+1;
     InnerNodes_[l]=i;
//      cout << "InnerNodes_[" << l << "]=" <<InnerNodes_[l] << "\n";
     for (int j=0; j < nbn_/2; ++j)
     {
       if ((i==BoundNodes_[j])||(i==PeriodicNodes_[j]))
       {
	 	l=l-1;
       }
     }
   }

   for (int i = 0; i < numnp_-nbn_; ++i)
   {
     MapHelper(i, 3+(nbn_/2)*ndf_+i*ndf_, InnerNodes_);
   }

//    for (int i=0; i< numnp_*ndf_; ++i)
//    {
//      for (int j=0; j<4; ++j)
//      {
//        cout << "Map_[" << i << "][" << j << "]=" << Map_[i][j] << " ";
//      }
//      cout << "\n";
//    }
//    cout << "\n";

   // set loadparameter to load (not temperature)
   // load also includes disp. control here...
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

   // open critical_eig_out_file
   string flnm1(Input.LastInputFileName());
   int pos1 = flnm1.find(".bfb");
   if (string::npos != pos1)
   {
      flnm1.erase(pos1,flnm1.length());
   }
   flnm1.append(".eigs.gpl");
   critical_eig_out_.open(flnm1.c_str() ,ios::out);
   eig_count_ = 0;

   //Open bloch analysis out file
   string flnm2(Input.LastInputFileName());
   int pos2 = flnm2.find(".bfb");
   if (string::npos != pos2)
   {
      flnm2.erase(pos2,flnm2.length());
   }
   flnm2.append(".bloch.gpl");
   bloch_wave_out_.open(flnm2.c_str() ,ios::out);
   bloch_count_ = 0;

   //Open out file for plotting
   flnm=(Input.LastInputFileName());
   pos = flnm.find(".bfb");
   if (string::npos != pos)
   {
      flnm.erase(pos,flnm.length());
   }
   flnm.append(".plot.gpl");
   plot_out_.open(flnm.c_str() ,ios::out);
   plot_count_ = 0;

   cout << "End of constructor \n";


}

void FEAP::MapHelper(int const i, int const offset, int const* const Nodes)
{
  Map_[Nodes[i]*ndf_][0] = 0;
  Map_[Nodes[i]*ndf_][1] = offset;
  Map_[Nodes[i]*ndf_][2] = 2;
  Map_[Nodes[i]*ndf_][3] = offset+1;
  FMap_[Nodes[i]*ndf_][0] = 0;
  FMap_[Nodes[i]*ndf_][1] = 1+offset;
  FMap_[Nodes[i]*ndf_][2] = 1;
  FMap_[Nodes[i]*ndf_][3] = 1+offset+1;

  Map_[Nodes[i]*ndf_+1][0] = 1;
  Map_[Nodes[i]*ndf_+1][1] = offset+1;
  Map_[Nodes[i]*ndf_+1][2] = 2;
  Map_[Nodes[i]*ndf_+1][3] = offset;
  FMap_[Nodes[i]*ndf_+1][0] = 3;
  FMap_[Nodes[i]*ndf_+1][1] = 1+offset+1;
  FMap_[Nodes[i]*ndf_+1][2] = 2;
  FMap_[Nodes[i]*ndf_+1][3] = 1+offset;
  if (ndf_>ndm_)
  {
    Map_[Nodes[i]*ndf_+2][0] = -1;
    Map_[Nodes[i]*ndf_+2][1] = -1;
    Map_[Nodes[i]*ndf_+2][2] = -1;
    Map_[Nodes[i]*ndf_+2][3] = -1;
    FMap_[Nodes[i]*ndf_+2][0] = -1;
    FMap_[Nodes[i]*ndf_+2][1] = -1;
    FMap_[Nodes[i]*ndf_+2][2] = -1;
    FMap_[Nodes[i]*ndf_+2][3] = -1;
    if (ndf_>(ndm_+1))
    {
      Map_[Nodes[i]*ndf_+3][0] = -1;
      Map_[Nodes[i]*ndf_+3][1] = -1;
      Map_[Nodes[i]*ndf_+3][2] = -1;
      Map_[Nodes[i]*ndf_+3][3] = -1;
      FMap_[Nodes[i]*ndf_+3][0] = -1;
      FMap_[Nodes[i]*ndf_+3][1] = -1;
      FMap_[Nodes[i]*ndf_+3][2] = -1;
      FMap_[Nodes[i]*ndf_+3][3] = -1;
    }
  }
}

void FEAP::UpdateValues(UpdateFlag flag) const
{
  Matrix Eye(2,2,0.0);
  Eye.SetIdentity(2);
   // disp gradient U is set with lambda and DOF
   //
   // Update FEAP solution vector
   UpdateDOF_F();
//    cout << "\n DOF_F_=" << setw(20) << DOF_F_ <<"\n";
   bfbfeap_set_nodal_solution_(&(DOF_F_[0]));

   UpdateJacobian();

   //Sum of S term for phantom energy term
   double S1 = 0.0;
   double S2 = 0.0;
   double S3 = 0.0;
   int shift;
   switch (LoadingType_)
   {
     case DEAD_LOAD:
     case PRESSURE_LOAD:
       shift = ndm_ * (ndm_ + 1) / 2;
       break;
     case DISPLACEMENT_CONTROL:
       shift = 0;
       break;
   }

   for (int i = 0; i < (numnp_-nbn_/2); ++i)
   {
      S1 += DOF_[shift+i*ndf_];
      S2 += DOF_[shift+1+i*ndf_];
      if(ndf_>(ndm_+1))
      S3 += DOF_[shift+2+i*ndf_];
   }
//    for (int i = nbn_; i < numnp_; ++i)
//    {
//       S1 += DOF_[3+(nbn_/2)*ndf_+(i-nbn_)*ndf_];
//       S2 += DOF_[4+(nbn_/2)*ndf_+(i-nbn_)*ndf_];
//    }


   bfbfeap_call_ener_(); // needs a "TPLOt" and "ENER" command in FEAP input file to work
   bfbfeap_get_potential_energy_(&(E0CachedValue_));
//    cout << "\n E0_=" << E0CachedValue_ <<"\n";
   switch (LoadingType_)
   {
     case PRESSURE_LOAD:
       E0CachedValue_ += Lambda_ * ((U_[0][0]*U_[1][1] - U_[0][1]*U_[1][0]) - 1.0) * nuc_ * CellArea_;
       break;
     case DEAD_LOAD:
       E0CachedValue_ += -Lambda_ * (Load_ * (U_ - Eye)).Trace() * nuc_ * CellArea_;
       break;
     case DISPLACEMENT_CONTROL:
       // nothing to do for disp. control loading energy
       break;
   }

   E0CachedValue_ += 1.0/eps_*(S1*S1 + S2*S2 + S3*S3);

   bfbfeap_call_form_();
   bfbfeap_get_reduced_residual_(&(E1CachedValue_F_[0]));

//      cout << "\n E1CachedValue_F_=" << setw(20) << E1CachedValue_F_ << "\n";

   // FEAP returns -E1, so fix it.
   for (int i = 0; i < E1CachedValue_F_.Dim(); ++i)
   {
         E1CachedValue_F_[i] = -E1CachedValue_F_[i];
   }
   E1CachedValue_ = Jacobian_.Transpose() * E1CachedValue_F_;
   W1CachedValue_ = FJacobian_.Transpose() * E1CachedValue_F_;
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     DispE1CachedValue_ = DispJacobian_.Transpose() * E1CachedValue_F_;
   }

   switch (LoadingType_)
   {
     case PRESSURE_LOAD:
       E1CachedValue_[0] += Lambda_ * U_[1][1] * nuc_ * CellArea_;
       E1CachedValue_[1] += Lambda_ * U_[0][0] * nuc_ * CellArea_;
       E1CachedValue_[2] -= Lambda_ * (U_[0][1]*sqrt(2.0)) * nuc_ * CellArea_;
       break;
     case DEAD_LOAD:
       E1CachedValue_[0] += -Lambda_ * Load_[0][0] * nuc_ * CellArea_;
       E1CachedValue_[1] += -Lambda_ * Load_[1][1] * nuc_ * CellArea_;
       E1CachedValue_[2] += -sqrt(2.0)/2.0*Lambda_ * (Load_[1][0]+Load_[0][1]) * nuc_ * CellArea_;
       break;
     case DISPLACEMENT_CONTROL:
       // nothing to do for disp. control loading term
       break;
   }

   for (int i = 0; i < (numnp_-nbn_/2); ++i)
   {
     E1CachedValue_[shift+i*ndf_] += 2.0/eps_*S1;
     W1CachedValue_[4+i*ndf_] += 2.0/eps_*S1;
     E1CachedValue_[shift+1+i*ndf_] += 2.0/eps_*S2;
     W1CachedValue_[4+1+i*ndf_] += 2.0/eps_*S2;
     if(ndf_>(ndm_+1))
     {
       E1CachedValue_[shift+2+i*ndf_] += 2.0/eps_*S3;
       W1CachedValue_[4+2+i*ndf_] += 2.0/eps_*S3;
     }
   }

   Cached_[0] = 1;
   Cached_[1] = 1;
   EvaluationCount_[0]++;

   if (flag == NeedStiffness)
   {
     //E2
     bfbfeap_call_tang_();
     bfbfeap_get_reduced_tang_(&(E2CachedValue_F_[0][0]));

//      cout << "\n E2CachedValue_F_=" << setw(20) << E2CachedValue_F_ << "\n";

     E2CachedValue_ = Jacobian_.Transpose() * E2CachedValue_F_ * Jacobian_;
     W2CachedValue_ = FJacobian_.Transpose() * E2CachedValue_F_ * FJacobian_;
//
//      cout << "\n jacobian=" << setw(20) << Jacobian_ << "\n";
//      cout << "\n E2CachedValue_=" << setw(20) << E2CachedValue_ << "\n";

     // This term doesn't contribute for displacement control
     if (LoadingType_ != DISPLACEMENT_CONTROL)
     {
       for (int i = 0; i < DOFS_; ++i)
       {
         for (int j = 0; j < DOFS_; ++j)
         {
           for (int k = 0; k < DOFS_F_; ++k)
           {
             if ((Map_[k][0]==i && Map_[k][1]==j) || (Map_[k][0]==j && Map_[k][1]==i))
               E2CachedValue_[i][j] += E1CachedValue_F_[k];
             if ((Map_[k][2]==i && Map_[k][3]==j) || (Map_[k][2]==j && Map_[k][3]==i))
               E2CachedValue_[i][j] += E1CachedValue_F_[k] / sqrt(2.0);
           }
         }
       }
       for (int i = 0; i < DOFS_+1; ++i)
       {
         for (int j = 0; j < DOFS_+1; ++j)
         {
           for (int k = 0; k < DOFS_F_; ++k)
           {
             if ((FMap_[k][0]==i && FMap_[k][1]==j) || (FMap_[k][0]==j && FMap_[k][1]==i))
               W2CachedValue_[i][j] += E1CachedValue_F_[k];
             if ((FMap_[k][2]==i && FMap_[k][3]==j) || (FMap_[k][2]==j && FMap_[k][3]==i))
               W2CachedValue_[i][j] += E1CachedValue_F_[k];
           }
         }
       }
     }
     // Phantom Energy Term for E2
     for (int i = 0; i < (numnp_-nbn_/2); ++i)
     {
       for (int j = 0; j < (numnp_-nbn_/2); ++j)
       {
         E2CachedValue_[shift+i*ndf_][shift+j*ndf_] += 2.0/eps_;
         W2CachedValue_[4+i*ndf_][4+j*ndf_] += 2.0/eps_;
         E2CachedValue_[shift+1+i*ndf_][shift+1+j*ndf_] += 2.0/eps_;
         W2CachedValue_[4+1+i*ndf_][4+1+j*ndf_] += 2.0/eps_;
         if(ndf_>(ndm_+1))
         {
           E2CachedValue_[shift+2+i*ndf_][shift+2+j*ndf_] += 2.0/eps_;
           W2CachedValue_[4+2+i*ndf_][4+2+j*ndf_] += 2.0/eps_;
         }
       }
     }
     // Loading term for E2 //
     if (LoadingType_ == PRESSURE_LOAD)
     {
       E2CachedValue_[0][1] += Lambda_ * nuc_ * CellArea_;
       E2CachedValue_[1][0] += Lambda_ * nuc_ * CellArea_;
       E2CachedValue_[2][2] -= Lambda_ * nuc_ * CellArea_;
     }

     //E1DLoad
     switch (LoadingType_)
     {
       case PRESSURE_LOAD:
         E1DLoadCachedValue_[0] = U_[1][1] * nuc_ * CellArea_;
         E1DLoadCachedValue_[1] = U_[0][0] * nuc_ * CellArea_;
         E1DLoadCachedValue_[2] = -U_[0][1] * nuc_ * CellArea_;
         break;
       case DEAD_LOAD:
         E1DLoadCachedValue_[0] = -Load_[0][0]* nuc_ * CellArea_;
         E1DLoadCachedValue_[1] = -Load_[1][1]* nuc_ * CellArea_;
         E1DLoadCachedValue_[2] = -sqrt(2.0)/2.0*(Load_[1][0]+Load_[0][1])* nuc_ * CellArea_;
         break;
       case DISPLACEMENT_CONTROL:
         int ii;
         for (int i = 0; i < DOFS_; ++i)
         {
           E1DLoadCachedValue_[i] = 0.0;
           for (int p = 0; p < DOFS_F_; ++p)
           {
             for (int q = 0; q < DOFS_F_; ++q)
             {
               E1DLoadCachedValue_[i] += Jacobian_[p][i] * E2CachedValue_F_[p][q] *
                   (DispJacobian_[q][0] + StretchRatio_ *DispJacobian_[q][1]);
             }
           }
         }
         break;
     }

     Cached_[2] = 1;
     Cached_[3] = 1;
     EvaluationCount_[1]++;

//     Vector E1DLoad_F = Jacobian_ * E1DLoadCachedValue_;
//      Matrix A(DOFS_F_, DOFS_F_+1);
//      for (int i = 0; i < DOFS_F_; ++i)
//       {
//          for (int j = 0; j < DOFS_F_+1; ++j)
//          {
//             A[i][j] = E2CachedValue_F_[i][j];
//          }
//          A[i][DOFS_F_] = E1DLoad_F[i];
//       }
//      Matrix Q(DOFS_F_+1, DOFS_F_+1);
//   	 Matrix R(DOFS_F_+1, DOFS_F_);
//      QR(A, Q, R, 1);
//       cout << "\n E1DLoad_F=" << setw(20) << E1DLoad_F << "\n";
//   	 cout << "\n Q=" << setw(20) << Q;

//      cout << "E1DLoadCachedValue_=" << setw(20) << E1DLoadCachedValue_ << "\n";
   }
}

void FEAP::UpdateDOF_F() const
{
  int shift;
  switch (LoadingType_)
  {
    case DEAD_LOAD:
    case PRESSURE_LOAD:
      shift = ndm_ * (ndm_ + 1) / 2;
      break;
    case DISPLACEMENT_CONTROL:
      shift = 0;
      break;
  }
   int ii,jj;
   for (int i=0; i< nbn_/2; ++i)
   {
      ii = BoundNodes_[i]* ndf_;
      jj = (i%(nbn_/2))*ndf_;
      // Recall: F_ = U_
      DOF_F_[ii]=F_[0][0]*(X_F_[ii] + DOF_[shift+jj]) + F_[0][1]*(X_F_[ii+1] + DOF_[shift+1+jj]);
      DOF_F_[ii+1]=F_[1][1]*(X_F_[ii+1] + DOF_[shift+1+jj]) + F_[1][0]*(X_F_[ii] + DOF_[shift+jj]);
      if (ndf_>ndm_)  //1 Extra dof : theta
      {
        DOF_F_[ii+2]=DOF_[shift+2+jj];
        if (ndf_>(ndm_+1))  //1 Extra dof : u',v'
        {
          DOF_F_[ii+3]=DOF_[shift+3+jj];
        }
      }
   }
   for (int i=nbn_/2; i< nbn_; ++i)
   {
      ii = PeriodicNodes_[i-nbn_/2]* ndf_;
      jj = (i%(nbn_/2))*ndf_;
      // Recall: F_ = U_
      DOF_F_[ii]=F_[0][0]*(X_F_[ii] + DOF_[shift+jj]) + F_[0][1]*(X_F_[ii+1] + DOF_[shift+1+jj]);
      DOF_F_[ii+1]=F_[1][1]*(X_F_[ii+1] + DOF_[shift+1+jj]) + F_[1][0]*(X_F_[ii] + DOF_[shift+jj]);
      if (ndf_>ndm_)  //1 Extra dof : theta
      {
        DOF_F_[ii+2]=DOF_[shift+2+jj];
        if (ndf_>(ndm_+1))  //1 Extra dof : u', v'
        {
          DOF_F_[ii+3]=DOF_[shift+3+jj];
        }
      }
   }
   int offst = shift + nbn_ / 2 * ndf_;
   for (int i = nbn_; i < numnp_; ++i)
   {
      ii = InnerNodes_[i-nbn_]* ndf_;
      jj = offst+(i-nbn_)*ndf_;
      // Recall: F_ = U_
      DOF_F_[ii]=F_[0][0]*(X_F_[ii]+DOF_[jj]) + F_[0][1]*(X_F_[ii+1]+DOF_[jj+1]);
      DOF_F_[ii+1]=F_[1][0]*(X_F_[ii]+DOF_[jj]) + F_[1][1]*(X_F_[ii+1]+DOF_[jj+1]);
      if (ndf_>ndm_)  //1 Extra dof : theta
      {
        DOF_F_[ii+2]=DOF_[jj+2];
        if (ndf_>(ndm_+1))  //1 Extra dof : u', v'
        {
          DOF_F_[ii+3]=DOF_[jj+3];
        }
      }
   }
   DOF_F_ -= X_F_;

}

void FEAP::UpdateJacobian() const
{
  int shift;
  switch (LoadingType_)
  {
    case DEAD_LOAD:
    case PRESSURE_LOAD:
      shift = ndm_ * (ndm_ + 1) / 2;
      break;
    case DISPLACEMENT_CONTROL:
      shift = 0;
      break;
  }
   int ii, jj;
   for (int i = 0; i < nbn_/2; ++i)
   {
      ii = BoundNodes_[i]*ndf_;
      jj = (i%(nbn_/2))*ndf_;
      JacobianHelper(ii, jj, shift);
   }

   for (int i = nbn_/2; i < nbn_; ++i)
   {
      ii = PeriodicNodes_[i-nbn_/2]*ndf_;
      jj = (i%(nbn_/2))*ndf_;
      JacobianHelper(ii, jj, shift);
   }

   for (int i = 0; i < numnp_-nbn_; ++i)
   {
      ii = InnerNodes_[i]*ndf_;
      jj = nbn_/2*ndf_+(i)*ndf_;
      JacobianHelper(ii, jj, shift);
   }
}


void FEAP::JacobianHelper(int const ii, int const jj, int const shift) const
{
  if (LoadingType_ != DISPLACEMENT_CONTROL)
  {
    Jacobian_[ii][0] = X_F_[ii] + DOF_[shift+jj];
    Jacobian_[ii][2] = 1.0/sqrt(2.0)*(X_F_[ii+1] + DOF_[shift+1+jj]);
    FJacobian_[ii][0] = X_F_[ii] + DOF_[shift+jj];
    FJacobian_[ii][1] = X_F_[ii+1] + DOF_[shift+1+jj];

    Jacobian_[ii+1][1] = X_F_[ii+1] + DOF_[shift+1+jj];
    Jacobian_[ii+1][2] = 1.0/sqrt(2.0)*(X_F_[ii] + DOF_[shift+jj]);
    FJacobian_[ii+1][2] = X_F_[ii] + DOF_[shift+jj];
    FJacobian_[ii+1][3] = X_F_[ii+1] + DOF_[shift+1+jj];
  }
  else
  {
    DispJacobian_[ii][0] = X_F_[ii] + DOF_[shift+jj];
    DispJacobian_[ii][2] = 1.0/sqrt(2.0)*(X_F_[ii+1] + DOF_[shift+1+jj]);

    DispJacobian_[ii+1][1] = X_F_[ii+1] + DOF_[shift+1+jj];
    DispJacobian_[ii+1][2] = 1.0/sqrt(2.0)*(X_F_[ii] + DOF_[shift+jj]);
  }
  Jacobian_[ii][shift+jj]=U_[0][0];
  Jacobian_[ii][shift+1+jj]=U_[0][1];
  FJacobian_[ii][1+shift+jj]=F_[0][0];
  FJacobian_[ii][1+shift+1+jj]=F_[0][1];

  Jacobian_[ii+1][shift+jj]=U_[1][0];
  Jacobian_[ii+1][shift+1+jj]=U_[1][1];
  FJacobian_[ii+1][1+shift+jj]=F_[1][0];
  FJacobian_[ii+1][1+shift+1+jj]=F_[1][1];

  if (ndf_>ndm_)
  {
    Jacobian_[ii+2][shift+2+jj]=1.0;
    FJacobian_[ii+2][1+shift+2+jj]=1.0;
    if (ndf_>(ndm_+1))
    {
      Jacobian_[ii+3][shift+3+jj]=1.0;
      FJacobian_[ii+3][1+shift+3+jj]=1.0;
    }
  }
}

int FEAP::CriticalPointInfo(int* const CPCrossingNum, int const& TFIndex, Vector const& DrDt,
                          int const& CPorBif, int const& NumZeroEigenVals, double const& Tolerance,
                          int const& Width, PerlInput const& Input, ostream& out)
{

   int Bif;
   double pi = 4.0*atan(1.0);
   // do standard CPInfo stuff and output bfb restart file
   Bif = Lattice::CriticalPointInfo(CPCrossingNum, TFIndex, DrDt, CPorBif, NumZeroEigenVals,
                                    Tolerance, Width, Input, out);
   if (CPorBif >= 0 )
   {
      Matrix
         D2 = E2(),
         EigVec,
         EigVal = SymEigVal(D2, &EigVec);
      int count = 0;

      // Find the modes
      int* Ind = new int[DOFS_];
      for (int i = 0; i < DOFS_; i++)
      {
         Ind[i] = 0;
         if (fabs(EigVal[0][i]) < Tolerance)
         {
            Ind[count++] = i;
         }
      }

      critical_eig_out_ << "# Critical Point Crossing Number: " << TFIndex << endl;

      for (int i = 0; i < count; ++i)
      {
         critical_eig_out_ << "####  Index = " << eig_count_ << "\n" <<" # Mode[" << i << "] = " << Ind[i] << endl;
         for (int k = 0; k < DOFS_; ++k)
         {
            DOF_[k] = EigVec[k][Ind[i]];
         }
         UpdateDOF_F();

         for (int i = 0; i < numnp_; ++i)
         {
            critical_eig_out_ << setw(Width_) << X_[i*ndm_+0]
                              << setw(Width_) << X_[i*ndm_+1]
                              << "\t" << " # Node " << i+1 << endl;
            critical_eig_out_ << setw(Width_) << X_[i*ndm_+0] + 0.1*DOF_F_[i*ndf_+0]
                              << setw(Width_) << X_[i*ndm_+1] + 0.1*DOF_F_[i*ndf_+1]
                              << endl << endl;
         }

         critical_eig_out_ << endl;
         ++eig_count_;
      }

      delete[] Ind;
   }
   else if ((CPorBif < 0 ) && (TFType_ == 1)) //Bloch Wave Analysis
   {

      for (int i = 0; i < KSpaceResolution_+1; ++i)
      {
         for (int j = 0; j < KSpaceResolution_+1; ++j)
         {
            K_[0] = pi*(((double) i)/((double) KSpaceResolution_));
            K_[1] = pi*(((double) j)/((double) KSpaceResolution_));

            DynamicalMatrixBis(K_);

                  CMatrix EigVec;
               Matrix DkEigVal = HermiteEigVal(Dk_, &EigVec);
               int NumZeroEig = 0; //For treating the k = [0,0] case

               //Dk has 2 zero eigen values when k = [0,0], so special treatment
               if ((i==0) && (j==0))
               {
                  double max = DkEigVal.MaxElement();

                  for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
                  {
                     if ((fabs(DkEigVal[0][l]) < Tolerance) && (NumZeroEig < 2))
                     {
                        DkEigVal[0][l] = max;
                        ++NumZeroEig;
                     }
                  }
               }
               int count = 0;
               int* Ind = new int[DOFS_];
               for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
               {
                  if (fabs(DkEigVal[0][l]) < Tolerance)
                     Ind[count++] = l;
               }
               double min = DkEigVal.MinElement();
               if (count>0)
               {
                  cout << "Found zero eigen value in Bloch wave analysis at k = [ " << K_[0] << ", " << K_[1] << " ] \n"
                       << " Eigenvalue : " << min << "\n";
                  for (int k = 0; k < count; ++k)
                  {
                     cout << "Tangent => [";
                     for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)

		    {
                         if ((fabs(sin(K_[0])) < 1.0e-5) && (fabs(sin(K_[1])) < 1.0e-5))
                         {
                            cout << EigVec[l][Ind[k]].real() << ", ";
                         }
                         else
                         {
                            cout << EigVec[l][Ind[k]] << ", ";
                         }

                     }
                     cout << "] \n";

              /*       if ((fabs(sin(K_[0])) < 1.0e-5) && (fabs(sin(K_[1])) < 1.0e-5))
                     {
                        cout << "Tangent for bif path : [";
                        for (int l = 0; l = 2*ndf_*(numnp_ - nbn_/2); ++l)
                        {
                            cout << pow(-1,l/(ndf_*(numnp_-nbn_/2))) * EigVec[l%(ndf_*(numnp_-nbn_/2))][Ind[k]].real() << ", ";
                        }
                        cout << "] \n";
                     } */

                  }



                  out << "============================================="
                      << "Found zero eigen value in Bloch wave analysis at k = [ " << K_[0] << ", " << K_[1] << " ] \n"
                      << " Eigenvalue : " << min << "\n";

               }

               delete[] Ind;
         }
      }


   NumExtraTFs_ = 0;
   }
   else if ((CPorBif < 0 ) && (TFType_ == 3))
   {
      for (int i = 0; i < KVectorMatrix_.Rows(); ++i)
      {
         K_[0] = 2*pi*KVectorMatrix_[i][0];
         K_[1] = 2*pi*KVectorMatrix_[i][1];

         DynamicalMatrixBis(K_);
         CMatrix EigVec;
         Matrix DkEigVal = HermiteEigVal(Dk_, &EigVec);

         int count = 0;
         int* Ind = new int[DOFS_];
         for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
         {
            if (fabs(DkEigVal[0][l]) < Tolerance)
               Ind[count++] = l;
         }
         double min = DkEigVal.MinElement();
         if (count>0)
         {
            cout << "Found zero eigen value in Bloch wave analysis at k = [ " << K_[0] << ", " << K_[1] << " ] \n"
                 << " Eigenvalue : " << min << "\n";
            for (int k = 0; k < count; ++k)
            {
               cout << "Tangent : RealPart :";
               for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
               {
                     cout << EigVec[l][Ind[k]].real() << ", ";
               }
               cout << "\n";
               cout << "Tangent : ImagPart :";
               for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
               {
                  cout << EigVec[l][Ind[k]].imag() << ", ";
               }

               cout << "\n";
            }
         }

         delete[] Ind;
      }
   }
   return -1;
}

void FEAP::ExtraTestFunctions(Vector& TF) const
{

   double pi = 4.0*atan(1.0);
   if (TFType_ != 4)
   {
     bloch_wave_out_ << "## Index Number : " << bloch_count_ << "\n"
                     << "#### Lambda = " << Lambda_;
     switch (LoadingType_)
     {
       case DEAD_LOAD:
       case PRESSURE_LOAD:
         bloch_wave_out_ << "  U11 = " << DOF_[0]
                         << "  U22 = " << DOF_[1]
                         << "  U12 = " << DOF_[2] / sqrt(2.0) << endl;
         break;
       case DISPLACEMENT_CONTROL:
         bloch_wave_out_ << "  U11 = " << Lambda_
                         << "  U22 = " << StretchRatio_ * Lambda_
                         << "  U12 = " << 0.0 << endl;
         break;
     }
   }

   if(TFType_ == 1) // Bloch Wave Analysis
   {
      int k = 0;
      for (int i = 0; i < KSpaceResolution_+1; ++i)
      {
         for (int j = 0; j < KSpaceResolution_+1; ++j)
         {
            K_[0] = pi*(((double) i)/((double) KSpaceResolution_));
            K_[1] = pi*(((double) j)/((double) KSpaceResolution_));

            DynamicalMatrixBis(K_);

               Matrix DkEigVal = HermiteEigVal(Dk_);
               int NumZeroEig = 0; //For treating the k = [0,0] case

               for (int l = 0; l < ndf_*(numnp_-nbn_/2); ++l)
               {
                  //Dk has 2 zero eigen values when k = [0,0], so special treatment
                  if ((i==0) && (j==0) && (NumZeroEig < 2) && (fabs(DkEigVal[0][l]) < Tolerance_))
                  {
                     double max = DkEigVal.MaxElement();
                     DkEigVal[0][l] = max;
                     TF[k] = max;
                     ++NumZeroEig;
                  }
                  else
                  {
                     TF[k] = DkEigVal[0][l];
                  }
                  ++k;

               }

            bloch_wave_out_ << setw(Width_) << K_[0]/(2*pi)
                            << setw(Width_) << K_[1]/(2*pi)
                            << setw(Width_) << DkEigVal.MinElement() << "\n";


         }
         bloch_wave_out_ << endl;
      }
      bloch_wave_out_ << "\n" << "\n";
      bloch_count_++;

   }
   else if (TFType_ == 2)
   {

      KDirection_ /= KDirection_.MaxElement();
      int k = 0;

      for (int i = 0; i < KSpaceResolution_+1; ++i)
      {
         K_[0] = 2*pi*(KRange_[0] + (KRange_[1] - KRange_[0])*(((double) i)/((double) KSpaceResolution_)))* KDirection_[0];
         K_[1] = 2*pi*(KRange_[0] + (KRange_[1] - KRange_[0])*(((double) i)/((double) KSpaceResolution_)))* KDirection_[1];

         DynamicalMatrixBis(K_);
         CMatrix DkEigVec(Dk_.Rows(),Dk_.Cols());
         Matrix DkEigVal = HermiteEigVal(Dk_, &DkEigVec);
         double min = DkEigVal.MinElement();
         int foundmin=0;
         int minIndex=0;

         for (int  l=0; l < ndf_*(numnp_-nbn_/2); ++l)
         {
            int NumZeroEig = 0;
            if ((K_[0]==0.0) && (K_[1]==0.0) && (NumZeroEig < 2) && (fabs(DkEigVal[0][l]) < Tolerance_))
            {
               double max = DkEigVal.MaxElement();
               DkEigVal[0][l] = max;
               TF[k] = max;
               ++NumZeroEig;
            }
            else
            {
               TF[k]=DkEigVal[0][l];
            }
            ++k;

            if (!foundmin && DkEigVal[0][l] == min)
            {
               DkEigVal[0][l] += DkEigVal.MaxElement();
               foundmin = 1;
            }
         }
         bloch_wave_out_ << setw(Width_) << K_[0]/(2*pi)
                         << setw(Width_) << K_[1]/(2*pi)
                         << setw(Width_) << min
                         << setw(Width_) << DkEigVal.MinElement() << "\n";

      }
      bloch_wave_out_ << "\n" << "\n";
      bloch_count_++;
   }
   else if (TFType_ == 3)
   {
      int k = 0;
      for (int i = 0; i < KVectorMatrix_.Rows(); ++i)
      {
         K_[0] = 2*pi*KVectorMatrix_[i][0];
         K_[1] = 2*pi*KVectorMatrix_[i][1];
         DynamicalMatrixBis(K_);
         Matrix DkEigVal = HermiteEigVal(Dk_);
         double min = DkEigVal.MinElement();
         int foundmin=0;

         for (int  l=0; l < ndf_*(numnp_-nbn_/2); ++l)
         {

            int NumZeroEig = 0;
            if ((K_[0]==0.0) && (K_[1]==0.0) && (NumZeroEig < 2) && (fabs(DkEigVal[0][l]) < Tolerance_))
            {
               double max = DkEigVal.MaxElement();
               DkEigVal[0][l] = max;
               TF[k] = max;
               ++NumZeroEig;
            }
            else
            {
               TF[k]=DkEigVal[0][l];
            }
            ++k;

            if (!foundmin && DkEigVal[0][l] == min)
            {
               DkEigVal[0][l] += DkEigVal.MaxElement();
               foundmin = 1;
            }
         }
         bloch_wave_out_ << setw(Width_) << K_[0]/(2*pi)
                         << setw(Width_) << K_[1]/(2*pi)
                         << setw(Width_) << min
                         << setw(Width_) << DkEigVal.MinElement() << "\n";
         bloch_count_++;

      }
   }
   else if(TFType_ == 4) // LoadingParameter
   {
      TF[0] = TFLoad_ - Lambda_;
   }
   return;
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

   Lambda_ = load + 10.0 * Tolerance_;
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     U_[0][0] = Lambda_;
     U_[1][1] = StretchRatio_ * Lambda_;
     U_[0][1] = 0.0;
     U_[1][0] = 0.0;
   }
   for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   stiffdl_static = E2();
   Lambda_ = load - 10.0 * Tolerance_;
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     U_[0][0] = Lambda_;
     U_[1][1] = StretchRatio_ * Lambda_;
     U_[0][1] = 0.0;
     U_[1][0] = 0.0;
   }
   for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   stiffdl_static -= E2();
   stiffdl_static /= 2.0 * Tolerance_;

   Lambda_ = load;
   if (LoadingType_ == DISPLACEMENT_CONTROL)
   {
     U_[0][0] = Lambda_;
     U_[1][1] = StretchRatio_ * Lambda_;
     U_[0][1] = 0.0;
     U_[1][0] = 0.0;
   }
   for (int i = 0; i < cachesize; ++i)
   {
      Cached_[i] = 0;
   }
   return stiffdl_static;
      // ajout ///////////////////////////////////////////////////////////////////
   //cout << "StiffnessDL= " << setw(W) << StiffnessDL();
      // ajout ///////////////////////////////////////////////////////////////////
}

void FEAP::print_gpl_config(fstream& out) const
{
   out << "# Index num: " << config_count_ << endl;
   out << "########## configuration for Lambda = " << Lambda_
       << " ##########" << endl;
   for (int i=0;i<nel_;++i)
   {
     if (ndf_>ndm_)
     {
      // node numbers for element i (-1 for zero-based values)
      int n1 = elmConn_[i*nen1_+0] - 1;
      int n2 = elmConn_[i*nen1_+1] - 1;

      out << setw(Width_) << X_[n1*ndm_+0] + DOF_F_[n1*ndf_+0]
          << setw(Width_) << X_[n1*ndm_+1] + DOF_F_[n1*ndf_+1]
          << setw(Width_) << DOF_F_[n1*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n1+1 << endl;
      out << setw(Width_) << X_[n2*ndm_+0] + DOF_F_[n2*ndf_+0]
          << setw(Width_) << X_[n2*ndm_+1] + DOF_F_[n2*ndf_+1]
          << setw(Width_) << DOF_F_[n2*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n2+1
          << endl << endl;
    }
    else
    {
      // node numbers for element i (-1 for zero-based values)
      int n1 = elmConn_[i*nen1_+0] - 1;
      int n2 = elmConn_[i*nen1_+1] - 1;
      int n3 = elmConn_[i*nen1_+2] - 1;
      int n4 = elmConn_[i*nen1_+3] - 1;
      out << setw(Width_) << X_[n1*ndm_+0] + DOF_F_[n1*ndf_+0]
          << setw(Width_) << X_[n1*ndm_+1] + DOF_F_[n1*ndf_+1]
          << setw(Width_) << DOF_F_[n1*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n1+1 << endl;
      out << setw(Width_) << X_[n2*ndm_+0] + DOF_F_[n2*ndf_+0]
          << setw(Width_) << X_[n2*ndm_+1] + DOF_F_[n2*ndf_+1]
          << setw(Width_) << DOF_F_[n2*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n2+1 << endl;
      out << setw(Width_) << X_[n3*ndm_+0] + DOF_F_[n3*ndf_+0]
          << setw(Width_) << X_[n3*ndm_+1] + DOF_F_[n3*ndf_+1]
          << setw(Width_) << DOF_F_[n3*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n3+1 << endl;
      out << setw(Width_) << X_[n4*ndm_+0] + DOF_F_[n4*ndf_+0]
          << setw(Width_) << X_[n4*ndm_+1] + DOF_F_[n4*ndf_+1]
          << setw(Width_) << DOF_F_[n4*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n4+1 << endl;
      out << setw(Width_) << X_[n1*ndm_+0] + DOF_F_[n1*ndf_+0]
          << setw(Width_) << X_[n1*ndm_+1] + DOF_F_[n1*ndf_+1]
          << setw(Width_) << DOF_F_[n1*ndf_+2]
          << "\t" << "# Element " << i+1 << ": Node " << n1+1
          << endl << endl;
    }
   }
   out << endl;

   ++config_count_;
}
/*
void FEAP::DynamicalMatrix(Vector const& K) const
{

   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[3]++;

   CMatrix A;
   Matrix B;
   Vector R;
   R.Resize(ndm_,0.0);
   A.Resize(ndf_,ndf_);
   B.Resize(ndf_,ndf_);
   int dim = ndf_*(numnp_ - nbn_/2);
   Dk_.Resize(dim,dim,0.0);


   double pi = 4.0 * atan(1.0);
   MyComplexDouble Ic(0.0, 1.0);

   for (int i = 0; i < numnp_; ++i)
   {
      for (int j = 0; j < numnp_; ++j)
      {
         if (!(((i>=nbn_/2) && (i < nbn_)) || ((j>=nbn_/2) && (j < nbn_))))
         {

            B = E2CachedValue_F_.Extract(ndf_*i,ndf_*j,ndf_);
            for (int k = 0; k < ndf_; ++k)
            {
               for (int l = 0; l < ndf_; ++l)
               {
                  A[k][l] = MyComplexDouble(B[k][l],0.0);
               }
            }
            int ii = i;
            int jj = j;
            if (i>=nbn_)
               ii = (i-nbn_/2);
            if (j>=nbn_)
               jj = (j-nbn_/2);
            Dk_.AddInsert(A,ii*ndf_,jj*ndf_);
         }
      }
   }


   for (int i = 0; i < N_rows_; ++i)
   {
      for (int j = 0; j < N_cols_/2-1; ++j)
      {
      if (N_[i][2+2*j]>0)
      {
         B = E2CachedValue_F_.Extract(ndf_*(N_[i][2+2*j]-1),ndf_*(N_[i][3+2*j]-1+nbn_/2),ndf_);

         for (int k = 0; k < ndf_; ++k)
         {
            for (int l = 0; l < ndf_; ++l)
            {
               A[k][l] = MyComplexDouble(B[k][l],0.0);
            }
         }
         R[0] = N_[i][0];
         R[1] = N_[i][1];

         A *= exp((K * R)*Ic );
         Dk_.AddInsert(A,(N_[i][2+2*j]-1-nbn_/2)*ndf_,(N_[i][3+2*j]-1)*ndf_);


         for (int k = 0; k < ndf_; ++k)
         {
            for (int l = 0; l < ndf_; ++l)
            {
               A[k][l] = MyComplexDouble(B[k][l],0.0);
            }
         }

         A *= exp(-(K*R)* Ic );
         Dk_.AddInsert(A.Transpose(),(N_[i][3+2*j]-1)*ndf_,(N_[i][2+2*j]-1-nbn_/2)*ndf_);

         B = E2CachedValue_F_.Extract(ndf_*(N_[i][3+2*j]-1+nbn_/2),ndf_*(N_[i][3+2*j]-1+nbn_/2),ndf_);
         for (int k = 0; k < ndf_; ++k)
         {
            for (int l = 0; l < ndf_; ++l)
            {
               A[k][l] = MyComplexDouble(B[k][l],0.0);
            }
         }
         Dk_.AddInsert(A,(N_[i][3+2*j]-1)*ndf_,(N_[i][3+2*j]-1)*ndf_);
      }
      }
   }

   return;
}*/

void FEAP::DynamicalMatrixBis(Vector const& K) const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[3]++;

   int W = 12;


   Vector R(ndm_,0.0);
   int dim = ndf_*(numnp_-nbn_/2);
   Matrix A(dim,dim,0.0);
   CMatrix B(ndf_,ndf_,0.0);
   Dk_.Resize(dim,dim,0.0);
   MyComplexDouble Ic(0.0,1.0);


   switch (LoadingType_)
   {
     case DEAD_LOAD:
     case PRESSURE_LOAD:
       A = E2CachedValue_.Extract(3,3,dim);
       break;
     case DISPLACEMENT_CONTROL:
       A = E2CachedValue_;
       break;
   }
   // Remove the Phantom Energy Term from E2
   // Note: pressure_load terms don't impact the part of E2 used here.
   for (int i = 0; i < numnp_-nbn_/2; ++i)
   {
      for (int j = 0; j < numnp_ - nbn_ / 2; ++j)
      {
          A[ndf_*i][ndf_*j] -= 2.0/eps_;
          A[1+ndf_*i][1+ndf_*j] -= 2.0/eps_;
          // Don't remove the 3rd DOF Phantom Energy term
          // it is not associated with a rigid-body mode
      }
   }

   for (int i = 0; i < dim; ++i)
   {
      for (int j = 0; j < dim; ++j)
      {
         Dk_[i][j] = MyComplexDouble(A[i][j],0.0);
      }
   }

   for (int i = 0; i < N_rows_; ++i)
   {
          R[0] = N_[i][0];
          R[1] = N_[i][1];

      for (int j = 0; j < N_cols_/2-1; ++j)
      {
          if (N_[i][2+2*j]>0)
          {
             Dk_.MultiplyBlock(exp( (K*R) * Ic ) , Map_[ndf_*(N_[i][2+2*j]-1)][1]-3, Map_[ndf_*(N_[i][3+2*j]-1)][1]-3,ndf_);
             Dk_.MultiplyBlock(exp(-(K*R) * Ic ) , Map_[ndf_*(N_[i][3+2*j]-1)][1]-3, Map_[ndf_*(N_[i][2+2*j]-1)][1]-3,ndf_);
          }
          if (j<N_cols_/2-2 && (N_cols_/2-1>1) && ndf_==2)
          {
             Dk_.MultiplyBlock(exp( (K*R) * Ic ) , Map_[ndf_*(N_[i][2+2*(j+1)]-1)][1]-3, Map_[ndf_*(N_[i][3+2*j]-1)][1]-3,ndf_);
             Dk_.MultiplyBlock(exp(-(K*R) * Ic ) , Map_[ndf_*(N_[i][3+2*j]-1)][1]-3, Map_[ndf_*(N_[i][2+2*(j+1)]-1)][1]-3,ndf_);
          }
          if (j>0 && (N_cols_/2-1>1) && N_[i][2+2*j]>0 && ndf_==2)
          {
             Dk_.MultiplyBlock(exp( (K*R) * Ic ) , Map_[ndf_*(N_[i][2+2*(j-1)]-1)][1]-3, Map_[ndf_*(N_[i][3+2*j]-1)][1]-3,ndf_);
             Dk_.MultiplyBlock(exp(-(K*R) * Ic ) , Map_[ndf_*(N_[i][3+2*j]-1)][1]-3, Map_[ndf_*(N_[i][2+2*(j-1)]-1)][1]-3,ndf_);
          }

       }
   }
   return;
}




int sortFunction(const void *a, const void *b)
{
   double doubleOne = *((double*)a);
   double doubleTwo = *((double*)b);
   if (doubleOne < doubleTwo)
      return -1;
   if (doubleOne == doubleTwo)
      return 0;
   return 1;
}

double FEAP::RankOneConvex(Matrix const& d2WdFdF) const
{
  double const twoPi = 8.0*atan(1.0);
  double const dt=twoPi/360.0;  // 1 degree
  double min = 1.0e100;

  double t = 0.0;
  while (t < twoPi)
  {
    Vector n(2);
    n[0] = cos(t);
    n[1] = sin(t);

    Matrix A(2,2, 0.0);
    for (int i=0; i<2; ++i)
      for (int j=0; j<2; ++j)
        for (int k=0; k<2; ++k)
          for (int l=0; l<2;++l)
            A[i][j] += d2WdFdF[i*2+k][j*2+l]*n[k]*n[l];

    Vector Eigs=SymEigVal(A);
    for (int i=0; i<2; ++i)
      if (Eigs[i] < min) min = Eigs[i];

    t += dt;
  }

  return min;
}


void FEAP::Print(ostream& out, PrintDetail const& flag,
                 PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions = 0;
   double engy;
   double E1norm;
   double ConjToLambdaScal;
   Vector ConjToLambda(3);
   Vector mintestfunct(NumTestFunctions());
   Vector TestFunctVals(NumTestFunctions());
   Matrix Eye(2,2,0.0);
   Eye.SetIdentity(2);

   Matrix Q(DOFS_, DOFS_);
   Matrix R(DOFS_, DOFS_-1);

   W = out.width();

   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }


   engy = E0();
   E1norm = E1().Norm();
   switch (LoadingType_)
   {
     case DEAD_LOAD:
       ConjToLambdaScal = (Load_*(U_-Eye)).Trace() * nuc_ * CellArea_;
       ConjToLambda[0] = DOF_[0];
       ConjToLambda[1] = DOF_[1];
       ConjToLambda[2] = DOF_[2];
       break;
     case PRESSURE_LOAD:
       ConjToLambdaScal = (U_[0][0]*U_[1][1] - U_[0][1]*U_[1][0] - 1.0)*nuc_*CellArea_;
       ConjToLambda[0] = DOF_[0];
       ConjToLambda[1] = DOF_[1];
       ConjToLambda[2] = DOF_[2];
       break;
     case DISPLACEMENT_CONTROL:
       ConjToLambda = DispE1CachedValue_;
       ConjToLambdaScal = ConjToLambda[0] + StretchRatio_ * ConjToLambda[1];
       break;
   }

   TestFunctions(TestFunctVals, LHS);
   mintestfunct = TestFunctVals;


   // check only the EigenValTFs
   for (int i = 0; i < DOFS_; ++i)
   {
      if ((UseEigenValTFs() == 1) && (TestFunctVals[i] < 0.0))
      {
         ++NoNegTestFunctions;
      }
   }
   // sort eigenvalues
   qsort(&(mintestfunct[0]), mintestfunct.Dim(), sizeof(double), sortFunction);
   int minprint = (NoNegTestFunctions+2 < mintestfunct.Dim()) ?
      NoNegTestFunctions+2 : mintestfunct.Dim();

   // Condense d2W/dFdF
   Matrix A=W2CachedValue_.Extract(0,0,4);
   Matrix C=W2CachedValue_.Extract(4,4,W2CachedValue_.Rows()-4);
   Matrix B(4,C.Cols());
   for (int i=0; i<4; ++i)
   {
     for (int j=0; j<B.Cols(); ++j)
     {
       B[i][j] = W2CachedValue_[i][4+j];
     }
   }
   Matrix Cond=A-B*C.Inverse()*B.Transpose();
   double RankOne = RankOneConvex(Cond);

   print_gpl_config(config_out_);
   switch (LoadingType_)
   {
     case DEAD_LOAD:
     case PRESSURE_LOAD:
       plot_out_ << setw(Width_) << Lambda_
                 << setw(Width_) << DOF_[0]
                 << setw(Width_) << DOF_[1]
                 << setw(Width_) << DOF_[2]/sqrt(2.0)
                 << setw(Width_) << E0CachedValue_ << endl;
      break;
     case DISPLACEMENT_CONTROL:
       plot_out_ << setw(Width_) << -ConjToLambda[0]/(nuc_*CellArea_)
                 << setw(Width_) << -ConjToLambda[1]/(nuc_*CellArea_)
                 << setw(Width_) << -ConjToLambdaScal/(2.0*nuc_*CellArea_)
                 << setw(Width_) << Lambda_
                 << setw(Width_) << StretchRatio_ * Lambda_
                 << setw(Width_) << 0.0
                 << setw(Width_) << -E0CachedValue_ << endl;
       break;
   }
   switch (flag)
   {
      case PrintLong:
         out << "FEAP:" << "\n" << "\n";
         switch (LoadingType_)
         {
           case PRESSURE_LOAD:
             out << "Using: In-plane Pressure loading.\n";
             break;
           case DEAD_LOAD:
             out << "Using: In-plane Dead loading.\n";
             break;
           case DISPLACEMENT_CONTROL:
             out << "Using: In-plane Displacement Control.\n";
             break;
         }
         out << "nuc_*CellArea: " << setw(W) << nuc_*CellArea_ << "\n";

         if (Echo_)
         {
            cout << "FEAP:" << "\n" << "\n";
            switch (LoadingType_)
            {
              case PRESSURE_LOAD:
                cout << "Using: In-plane Pressure loading.\n";
                break;
              case DEAD_LOAD:
                cout << "Using: In-plane Dead loading.\n";
                break;
              case DISPLACEMENT_CONTROL:
                cout << "Using: In-plane Displacement Control.\n";
                break;
            }
         }
         cout << "nuc_*CellArea: " << setw(W) << nuc_*CellArea_ << "\n";

         // passthrough to short
      case PrintShort:
         out << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "ConjToLambdaScal: " << setw(W) << ConjToLambdaScal << "\n"
             << "ConjToLambda: " << setw(W) << ConjToLambda << "\n"
             << "DOF: " << setw(W) << DOF_ << "\n"
             << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
             << "Potential Value: " << setw(W) << engy << "\n"
             << "Force Norm: " << setw(W) << E1norm << "\n";
      // ajout ///////////////////////////////////////////////////////////////////
      //  out << "StiffnessDL= " << setw(W) << StiffnessDL() << "\n";
      // ajout ///////////////////////////////////////////////////////////////////

         out << "Bifurcation Info: ";
         for (int i=0;i<minprint; ++i) out << setw(W) << mintestfunct[i];
         out << setw(W) << NoNegTestFunctions << "\n";
         out << "RankOneConvex Min Eig: " << setw(W) << RankOne << "\n";

         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
                 << "ConjToLambdaScal: " << setw(W) << ConjToLambdaScal << "\n"
                 << "ConjToLambda: " << setw(W) << ConjToLambda << "\n"
              //   << "DOF_F_: " << setw(W) << DOF_F_ << "\n"
                 << "DOF: " << setw(W) << DOF_ << "\n"
                 << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
                 << "Potential Value: " << setw(W) << engy << "\n"
                 << "Force Norm: " << setw(W) << E1norm << "\n";
       //ajout ///////////////////////////////////////////////////////////////////
//   cout << "StiffnessDL= " << setw(W) << E2() << "\n";
//   cout << "X_F_= " << setw(W) << X_F_ << "\n";
//      QR(E2(), Q, R, 1);
//   cout << "QR= " << setw(W) << Q << "\n";

       //   cout << "nodal_coord= " << setw(W) << X_ << "\n";
       //ajout ///////////////////////////////////////////////////////////////////
       //
            cout << "Bifurcation Info: ";
            for (int i=0;i<minprint; ++i) cout << setw(W) << mintestfunct[i];
            cout << setw(W) << NoNegTestFunctions << "\n";
            cout << "RankOneConvex Min Eig: " << setw(W) << RankOne << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out, FEAP& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}
