#include "MultiLatticeKIM.h"
#include "UtilityFunctions.h"
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

int const MultiLatticeKIM::DIM3 = 3;

double const RoEig_[3] = {1.0, 2.0, 3.0};
double const TrEig_[3] = {4.0, 5.0, 6.0};



MultiLatticeKIM::~MultiLatticeKIM()
{
   delete[] EulerAng_;
   delete[] BodyForce_;
   delete CBK_;
   int status;

   //KIM_API_model_destroy(pkim_, &status);
   status = KIM_API_model_destroy(pkim_);
   KIM_API_free(&pkim_, &status);
}

MultiLatticeKIM::MultiLatticeKIM(PerlInput const& Input, int const& Echo, int const& Width,
                                 int const& Debug) :
   Lattice(Input, Echo)
{
   dbg_ = Debug;
   // Get Lattice definition
   PerlInput::HashStruct Hash = Input.getHash("Lattice", "MultiLatticeKIM");
   // Set default values
   KillTranslations_ = 1; // 1-true, 0-false
   StiffnessYes_ = 0;
   int needKillRotations = 1;
   KillRotations_ = 0; // 0-do nothing, 1-kill one rotation, 2-kill three rotations
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
   if (!strcmp("SymLagrangeCB", CBKin))
   {
      CBK_ = new SymLagrangeCB(Input, &Hash);
      KillTranslations_ = 0;
      needKillRotations = 0;
   }
   else if (!strcmp("SymLagrangeWTransCB", CBKin))
   {
      CBK_ = new SymLagrangeWTransCB(Input, &Hash);
      needKillRotations = 0;
   }
   else if (!strcmp("LagrangeCB", CBKin))
   {
      CBK_ = new LagrangeCB(Input, &Hash);
   }
   else if (!strcmp("MixedCB", CBKin))
   {
      CBK_ = new MixedCB(Input, &Hash);
   }
   else if (!strcmp("EulerCB", CBKin))
   {
      CBK_ = new EulerCB(Input, &Hash);
   }
   else
   {
      cerr << "Error Unknown MultiLattice{CBKinematics}{Type} specified" << "\n";
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
      }
      else
      {
         cerr << "Error (MultiLatticeKIM()): Unknown RotationConstraint type" << "\n";
         exit(-2);
      }
   }
   // Setup Bodyforce_
   BodyForce_ = new Vector[InternalAtoms_];

   for (int i = 0; i < InternalAtoms_; ++i)
   {
      BodyForce_[i].Resize(DIM3, 0.0);
   }

   // Get Thermo parameters
   Tref_ = Input.getDouble(Hash, "Tref");
   // PhiRef_ = Input.getDouble(Hash,"PhiRef");
   // EntropyRef_ = Input.getDouble(Hash,"EntropyRef");
   // HeatCapacityRef_ = Input.getDouble(Hash,"HeatCapacityRef");

   //START OF KIM COMPLIANT ALTERATIONS
   int status;
   char* modelname = new char[100];
   if (Input.ParameterOK(Hash, "KIMModel"))
   {
      modelname = (char*) Input.getString(Hash, "KIMModel"); //Reads in the name of the model from Input file
   }
   else
   {
      cerr << "No KIM Model in input file" << "\n";
      exit(-1);
   }

   numberOfParticles_ = InternalAtoms_; //From CBKHash

   numberParticleTypes_ = Input.getInt(CBKHash, "NumberOfSpecies"); // Gets total number of species. From CBKHash

   const char** SpeciesList = new const char*[numberParticleTypes_];
   const char** AtomSpeciesList  = new const char*[InternalAtoms_] ;

   for (int i = 0;i < numberParticleTypes_;i++)
   {
      SpeciesList[i] = Input.getString(CBKHash, "SpeciesList", i); // Gets list of Species. From CBKHash
   }


   /* Create empty KIM object conforming to fields in the KIM descriptor files */
   /* of the Test and Model */
   Write_KIM_descriptor_file(SpeciesList, numberParticleTypes_);  //Calls function to create a compatible descriptor file. Will need to augment so that it can read in info from a model and automatically decides the appropriate tests.

   char* Test_Descriptor_file = new char[10000];
   strcpy(Test_Descriptor_file, descriptor_file_.str().c_str());

//      cout << "Test_Descriptor_file = " << endl;
//      cout << Test_Descriptor_file << endl;

   status =  KIM_API_string_init(&pkim_, Test_Descriptor_file, modelname);
   KIM_API_set_test_buffer(pkim_, this, &status);

   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, (char*) __FILE__,(char*) "Test-Model coupling failure (see kim.log file for details).", status);
      exit(1);
   }
   particleTypes_ = new int[InternalAtoms_];

   coords_ = new double[InternalAtoms_*3];
   forces_ = new double[InternalAtoms_*3];
   //Force_.Resize(InternalAtoms_*3);
   KIM_API_setm_data(pkim_, &status, 10*4,
                     "numberOfParticles",           1,   &numberOfParticles_,              1,
                     "numberParticleTypes",         1,   &numberParticleTypes_,            1,
                     "particleTypes",               InternalAtoms_,   &(particleTypes_[0]),          1,
                     "coordinates",                 3*InternalAtoms_, &(coords_[0]),                 1,
                     "forces",                      3*InternalAtoms_, &(forces_[0]),                 1,
                     "cutoff",                      1,   &cutoff_,                         1,
                     "get_neigh",                                  1, &MultiLatticeKIM::get_neigh,                                                 1,
                     "process_dEdr",                               1, &MultiLatticeKIM::process_dEdr,                      1,
                     "process_d2Edr2",                             1, &MultiLatticeKIM::process_d2Edr2,                      1,
                     "energy",                      1,   &energy_,                                           1);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, (char*) __FILE__,(char*) "KIM_API_setm_data",status);

   status = KIM_API_model_init(pkim_);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, (char*) __FILE__,(char*) "KIM_API_model_init", status);

   NTemp_ = 1.0;
   InfluenceDist_ = cutoff_;
   LatSum_(CBK_, InternalAtoms_, &InfluenceDist_, &NTemp_);

   KIM_API_setm_data(pkim_, &status, 1*4,
                     "neighObject",                                 1,   &LatSum_,                                                         1);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, (char*) __FILE__,(char*) "KIM_API_setm_data",status);

   int partcl_type_code;
   for (int i = 0; i < InternalAtoms_; i++)  //figures out which atoms corresponds to each element in the species list
   {
      AtomSpeciesList[i] = Input.getString(CBKHash, "AtomSpeciesKIM", i);
      for (int j = 0; j < numberParticleTypes_; j++)
      {
         if(!strcmp(AtomSpeciesList[i],SpeciesList[j]))
         {
            partcl_type_code = KIM_API_get_partcl_type_code(pkim_, (char*) SpeciesList[j], &status);
            if (KIM_STATUS_OK > status){
               KIM_API_report_error(__LINE__, (char*) __FILE__, (char*) "KIM_API_get_partcl_type_code", status);
               exit(1);}
            particleTypes_[i] = partcl_type_code;
            j=numberParticleTypes_+1;
         }
      }
   }

   // Get Lattice parameters
   if (Input.ParameterOK(Hash, "Density"))
   {
      Density_ = Input.getInt(Hash, "Density");
   }
   else
   {
      Density_ = Input.useInt(1, Hash, "Density"); // Default Value
   }
   NormModulus_ = Input.getDouble(Hash, "NormModulus");
   ConvexityDX_ = Input.getDouble(Hash, "ConvexityDX");
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
      cerr << "Unknown Loading Parameter" << "\n";
      exit(-1);
   }
   Lambda_ = 0.0;
   EulerAng_ = new double[DIM3];
   EulerAng_[0] = Input.getDouble(Hash, "EulerAngle_X");
   EulerAng_[1] = Input.getDouble(Hash, "EulerAngle_Y");
   EulerAng_[2] = Input.getDouble(Hash, "EulerAngle_Z");
   LoadingProportions_.Resize(DIM3);
   Input.getVector(LoadingProportions_, Hash, "LoadProportions");
   // Calculate Rotation and Loading
   // Euler angles transformation Rotation_ = Z*Y*X
   Rotation_.Resize(DIM3, DIM3, 0.0);
   Rotation_[0][0] = cos(EulerAng_[1]) * cos(EulerAng_[2]);
   Rotation_[0][1] = cos(EulerAng_[2]) * sin(EulerAng_[0]) * sin(EulerAng_[1])
      - cos(EulerAng_[0]) * sin(EulerAng_[2]);
   Rotation_[0][2] = cos(EulerAng_[0]) * cos(EulerAng_[2]) * sin(EulerAng_[1])
      + sin(EulerAng_[0]) * sin(EulerAng_[2]);
   Rotation_[1][0] = cos(EulerAng_[1]) * sin(EulerAng_[2]);
   Rotation_[1][1] = cos(EulerAng_[0]) * cos(EulerAng_[2]) + sin(EulerAng_[0])
      * sin(EulerAng_[1]) * sin(EulerAng_[2]);
   Rotation_[1][2] = -cos(EulerAng_[2]) * sin(EulerAng_[0])
      + cos(EulerAng_[0]) * sin(EulerAng_[1]) * sin(EulerAng_[2]);
   Rotation_[2][0] = -sin(EulerAng_[1]);
   Rotation_[2][1] = cos(EulerAng_[1]) * sin(EulerAng_[0]);
   Rotation_[2][2] = cos(EulerAng_[0]) * cos(EulerAng_[1]);
   //
   // Loading_ = R*Lambda*R^T
   Loading_.Resize(DIM3, DIM3, 0.0);
   for (int i = 0; i < DIM3; ++i)
   {
      for (int j = 0; j < DIM3; ++j)
      {
         for (int k = 0; k < DIM3; ++k)
         {
            Loading_[i][j] += Rotation_[i][k] * LoadingProportions_[k] * Rotation_[j][k];
         }
      }
   }

   // needed to initialize reference length
   int iter;
   iter = Input.getPosInt(Hash, "MaxIterations");
   GridSize_ = Input.getPosInt(Hash, "BlochWaveGridSize");

   // values to identifiy REFERENCE CONFIGURATION
   if (Input.ParameterOK(Hash, "ReferenceTemperature"))
   {
      REFTemp_ = Input.getDouble(Hash, "ReferenceTemperature");
   }
   else
   {
      REFTemp_ = Input.useDouble(1.0, Hash, "ReferenceTemperature"); // Default Value
   }
   if (Input.ParameterOK(Hash, "ReferenceLambda"))
   {
      REFLambda_ = Input.getDouble(Hash, "ReferenceLambda");
   }
   else
   {
      REFLambda_ = Input.useDouble(0.0, Hash, "ReferenceLambda"); // Default Value
   }

   NewCBCellFlag_ = 0; // make sure initialized
   PerlInput::HashStruct TFHash = Input.getHash(Hash, "ExtraTestFunctions");
   const char* TFtyp = Input.getString(TFHash, "Type");
   if ((!strcmp("None", TFtyp)) || (!strcmp("none", TFtyp)))
   {
      TFType_ = 0;
      NumExtraTFs_ = 0;
   }
   else if((!strcmp("KVectors", TFtyp)) || (!strcmp("kvectors", TFtyp)))
   {
      //KVector is input as [h,k,l, c, d] -> (c/d)(h,k,l)
      DynMatrixDim_ = DIM3*InternalAtoms_;
      NumKVectors_ = Input.getArrayLength(TFHash,"KVectors");
      KVectorMatrix_.Resize(NumKVectors_, 5);
      Input.getMatrix(KVectorMatrix_,TFHash,"KVectors");

      TFType_ = 1;
      NumExtraTFs_ = DynMatrixDim_*NumKVectors_;

      NewCBCellFlag_ = 1;  // default value
      if (Input.ParameterOK(TFHash, "PrintNewCBCell"))
      {
         const char* NewCell = Input.getString(TFHash, "PrintNewCBCell");
         if ((!strcmp("No", NewCell)) || (!strcmp("no", NewCell)))
         {
            NewCBCellFlag_ = 0;
         }
      }
   }
   else if((!strcmp("LoadingParameters", TFtyp)) || (!strcmp("loadingparameters", TFtyp)))
   {
      TFType_ = 2;
      NumExtraTFs_ = Input.getArrayLength(TFHash,"LoadingParameters");

      TFLoad_.Resize(NumExtraTFs_);
      Input.getVector(TFLoad_,TFHash,"LoadingParameters");
   }
   else
   {
      cerr << "Error (MultiLatticeKIM()): Unknown TestFunctions{Type}" << "\n";
      exit(-3);
   }

   // Initialize various data storage space
   ME1_static.Resize(CBK_->DOFS(), 0.0);
   ME2_static.Resize(CBK_->DOFS(), CBK_->DOFS(), 0.0);
   A_static.Resize(InternalAtoms_ * DIM3, InternalAtoms_ * DIM3);
   EigVals_static.Resize(1, InternalAtoms_ * DIM3);
   InverseLat_static.Resize(DIM3, DIM3);
   Z_static.Resize(DIM3);
   str_static.Resize(CBK_->DOFS());
   stiff_static.Resize(CBK_->DOFS(), CBK_->DOFS());
   CondEV_static.Resize(1, CBK_->Fsize());
   TE_static.Resize(CBK_->DOFS());
   CondModuli_static.Resize(CBK_->Fsize(), CBK_->Fsize());
   TestFunctVals_static.Resize(NumTestFunctions());
   if(TFType_ == 2) // only print stiffness eigenvalues
   {
      TestFunctVals_Print.Resize(CBK_->DOFS());
   }
   else // print everything
   {
      TestFunctVals_Print.Resize(TestFunctVals_static.Dim());
   }
   K_static.Resize(DIM3);

   if (Input.ParameterOK(Hash, "InitialEqbm"))
   {
      const char* init_equil = Input.getString(Hash, "InitialEqbm");
      if (!strcmp("Yes",init_equil) || !strcmp("yes",init_equil))
      {
         cout << "Constructor: HERE " << endl;
         int err = 0;
         err = FindLatticeSpacing(iter);
         if (err)
         {
            cerr << "unable to find initial lattice spacing!" << "\n";
            exit(-1);
         }
      }
   }
   else
   {
      cout << "Constructor: HERE2 " << endl;
      Input.useString("Yes", Hash, "InitialEqbm");
      int err = 0;
      err = FindLatticeSpacing(iter);
      if (err)
      {
         cerr << "unable to find initial lattice spacing!" << "\n";
         exit(-1);
      }
   }

   // Setup initial status for parameters
   NTemp_ = Input.getDouble(Hash, "NTemp");
   Lambda_ = Input.getDouble(Hash, "Lambda");

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   UCIter_(GridSize_);

   Input.EndofInputSection();

   Cached_[0] = 0;
   Cached_[1] = 0;
}

int MultiLatticeKIM::FindLatticeSpacing(int const& iter)
{
   cout << ":FindLattieSpacing::HERE " << endl;
   const double Tol = DOF().Dim()*1.0e-13;

   Lambda_ = REFLambda_;
   NTemp_ = REFTemp_;

   CBK_->SetReferenceDOFs();
   LatSum_.Recalc();

   if (Echo_)
   {
      RefineEqbm(Tol, iter, &cout);
   }
   else
   {
      RefineEqbm(Tol, iter, 0);
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

   LatSum_.Recalc();

   cout << "End of FindLatticeSpacing " << endl;
   return 0;
}

void MultiLatticeKIM::SetParameters(double const* const Vals, int const& ResetRef)
{
   cout << "ERROR IN MultiLatticeKIM::SetParameters. Exiting" << endl;
   exit(1);
}


void MultiLatticeKIM::UpdateValues() const
{
//      cout << "CHECKPOINT MultiLatticeKIM::UpdateValues() " << endl;
   int status;
   Vector coordsTemp = CBK_->CBKtoCoords();

   for(int i = 0;i < (3*InternalAtoms_); i++)
   {
      coords_[i] = coordsTemp[i];
   }

   LatSum_.Recalc();
//    cout << "UpdateValues:: StiffnessYes_ = " << StiffnessYes_ << endl;

//    KIM_API_print(pkim_, &status);

   status = KIM_API_model_compute(pkim_);

//    KIM_API_print(pkim_, &status);

   // Reset StiffnessYes_ in case we don't need to compute Stiffness
   StiffnessYes_ = 0;

   for(int i = 0;i < (InternalAtoms_); i++)
   {
      for (int j = 0; j < 3; j++)
      {
         BodyForce_[i][j] = forces_[i*3 + j];
      }
   }
}

// Lattice Routines
double MultiLatticeKIM::E0() const
{
   //cout <<"CHECKPOINT MultiLatticeKIM::E0()" << endl;
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

   E0CachedValue_ +=  0.5 * (TrEig_[0] * Tsq_static[0] + TrEig_[1] * Tsq_static[1] + TrEig_[2] * Tsq_static[2])    + 0.5 * (RoEig_[0] * Rsq_static[0] + RoEig_[1] * Rsq_static[1] + RoEig_[2] * Rsq_static[2]);

   return E0CachedValue_;
}


double MultiLatticeKIM::energy(LDeriv const& dl) const
{
   double phi = 0.0;

   double Vr = Density_ ? CBK_->RefVolume() : 1.0;

   if (dl == L0)
   {
      if (!Cached_[0])
      {
         UpdateValues();
      }

      phi = energy_ / (Vr * NormModulus_);

      for (int i = 0; i < 3; ++i)
      {
         for (int j = 0; j < 3; ++j)
         {
            phi -= Lambda_ * Loading_[i][j] * ((CBK_->DOF())[CBK_->INDF(j, i)] - Del(j, i));
         }
      }
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::energy()" << "\n";
      exit(-1);
   }

   return phi;
}

double MultiLatticeKIM::ConjugateToLambda() const
{
   cerr << "Error in MultiLatticeKIM::ConjugateToLambda() empty function \n ";
   exit(-45);
}

Vector const& MultiLatticeKIM::E1() const
{
//      cout <<"CHECKPOINT MultiLatticeKIM::E1()" << endl;
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

   return ME1_static;
}

//Vector const& MultiLatticeKIM::stress(TDeriv const& dt, LDeriv const& dl) const
Vector const& MultiLatticeKIM::stress(LDeriv const& dl) const
{
   for(int i = 0;i < CBK_->DOFS(); i++)
   {
      ME1_static[i] = 0.0;
   }

   if (dl == L0)
   {
      if (!Cached_[0])
      {
         UpdateValues();
      }
      double Vr = Density_ ? CBK_->RefVolume() : 1.0;
      ME1_static *= 1.0 / ((Vr * NormModulus_));

      //      cout << "MultiLatticeKIM -> ME1_static Pt.1 = " << setw(15) << ME1_static << endl;

      for (int i = 0; i < 3; ++i)
      {
         for (int j = 0; j < 3; ++j)
         {
            ME1_static[(CBK_->INDF(i,j))] -= Lambda_ * Loading_[j][i];
         }
      }
      //      cout << "MultiLatticeKIM -> ME1_static Pt.2 = " << setw(15) << ME1_static << endl;
   }
   else if (dl == DL)
   {
      for (int i = 0; i < 3; ++i)
      {
         for (int j = 0; j < 3; ++j)
         {
            ME1_static[(CBK_->INDF(i,j))] -= Loading_[j][i];
         }
      }
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stress()" << "\n";
      exit(-1);
   }

   //   cout << "MultiLatticeKIM -> ME1_static Pt.3 = " << setw(15) << ME1_static << endl;
   return ME1_static;
}

Matrix const& MultiLatticeKIM::E2() const
{
//      cout << "CHECKPOINT MultiLatticeKIM::E2()"  << endl;
   stiffness();

   if (KillTranslations_)
   {
      for (int i = 0; i < InternalAtoms_; ++i)
      {
         for (int j = 0; j < DIM3; ++j)
         {
            for (int k = 0; k < InternalAtoms_; ++k)
            {
               ME2_static[CBK_->INDS(i, j)][CBK_->INDS(k, j)] += TrEig_[j] / InternalAtoms_;
            }
         }
      }
   }

   return ME2_static;
}

Matrix const& MultiLatticeKIM::stiffness(LDeriv const& dl) const
{
   for(int i = 0;i < CBK_->DOFS(); i++)
   {
      for(int j = 0;j < CBK_->DOFS(); j++)
      {
         ME2_static[i][j] = 0.0;
      }
   }

   if (dl == L0)
   {
      if (!Cached_[0])
      {
         StiffnessYes_ = 1;
         UpdateValues();
      }
      double Vr = Density_ ? CBK_->RefVolume() : 1.0;
      ME2_static *= 1.0 / ((Vr * NormModulus_));

   }
   else if (dl == DL)
   {
      // Nothing to do: Phi2_static is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stiffness()" << "\n";
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

int MultiLatticeKIM::CriticalPointInfo(int* const CPCrossingNum, int const& TFIndex,
                                       Vector const& DrDt, int const& CPorBif,
                                       int const& NumZeroEigenVals, double const& Tolerance,
                                       int const& Width, PerlInput const& Input, ostream& out)
{
   cerr << "Error in MultiLatticeKIM::CriticalPointInfo() empty function \n ";
   exit(-45);
}

void MultiLatticeKIM::ExtraTestFunctions(Vector& TF) const
{
   if(TFType_ == 1) // KVectors
   {
      Vector KV1(DIM3,0.0);
      Vector KV2(DIM3,0.0);
      CMatrix DynMat(DynMatrixDim_,DynMatrixDim_,0.0);
      Matrix DynMatEigVal(1,DynMatrixDim_,0.0);
      int k=0;
      for (int i=0;i<NumKVectors_;++i)
      {
         for(int j=0;j<DIM3;++j)
         {
            KV1[j] = KVectorMatrix_[i][j]*KVectorMatrix_[i][3]/KVectorMatrix_[i][4];
         }
         KV2 = InverseLat_static*KV1;
         DynMat = ReferenceDynamicalStiffness(KV2);
         DynMatEigVal = HermiteEigVal(DynMat);

         for (int l=0; l< DynMatrixDim_;++l)
         {
            TF[k] = DynMatEigVal[0][l];
            ++k;
         }
      }
   }
   else if(TFType_ == 2) // LoadingParameters
   {
      for(int i = 0; i < NumExtraTFs_; ++i)
      {
         if (LoadParameter_ == Temperature)
         {
            TF[i] = (TFLoad_[i] - Temp());
         }
         else if (LoadParameter_ == Load)
         {
            TF[i] = (TFLoad_[i] - Lambda());
         }

         else
         {
            cerr << "Error MultiLatticeKIM::ExtraTestFunctions:  Unknown LoadParameter.\n";
            exit(-2);
         }
      }
   }
}

void MultiLatticeKIM::KPrint(int TFIndex, int Width, ostream& out) const
{
   int counter;
   int whichTF;
   int whichKV;
   Vector KTest(DIM3,0.0);
   Vector KVectorPrint(5,0.0);
   whichTF = TFIndex - (CBK_->DOFS());

   for(int i=0;i<NumKVectors_;++i)
   {
      counter = i*DynMatrixDim_;
      for (int j = counter; j<(counter + DynMatrixDim_);++j)
      {
         if (whichTF == j)
         {
            whichKV = i;
            break;
         }
      }
   }
   for (int i=0; i<5;i++)
   {
      KVectorPrint[i] = KVectorMatrix_[whichKV][i];
   }
   out << setw(Width) << KVectorPrint << endl;
}

void MultiLatticeKIM::TFCritPtInfo(int TFIndex, int Width, ostream& out) const
{
   Vector KVec1(DIM3,0.0);
   Vector KVec2(DIM3,0.0);
   Vector KVectorPrint(5,0.0);
   CMatrix DynMat(DynMatrixDim_,DynMatrixDim_,0.0);
   Matrix DynMatEigVal(1,DynMatrixDim_,0.0);
   Vector DynMatEigValPrint(DynMatrixDim_, 0.0);

   int counter;
   int whichTF;
   int whichKV;
   whichTF = TFIndex - (CBK_->DOFS());

   for(int i=0;i<NumKVectors_;++i)
   {
      counter = i*DynMatrixDim_;
      for (int j = counter; j<(counter + DynMatrixDim_);++j)
      {
         if (whichTF == j)
         {
            whichKV = i;
            break;
         }
      }
   }

   for(int j=0;j<DIM3;++j)
   {
      KVec1[j] = KVectorMatrix_[whichKV][j]*KVectorMatrix_[whichKV][3]/KVectorMatrix_[whichKV][4];
   }
   KVec2 = InverseLat_static*KVec1;
   DynMat = ReferenceDynamicalStiffness(KVec2);
   DynMatEigVal = HermiteEigVal(DynMat);

   for (int i=0; i<5;i++)
   {
      KVectorPrint[i] = KVectorMatrix_[whichKV][i];
   }

   // Print out info
   out << "\n";
   out <<  "$TestFunctions{KVector} = [" << KVectorPrint[0] << ", " << KVectorPrint[1]
       << ", " << KVectorPrint[2] << ", " << KVectorPrint[3] << ", " << KVectorPrint[4]
       <<  "];" << "\n";
   for (int i = 0; i < DynMatrixDim_;i++)
   {
      DynMatEigValPrint[i] = DynMatEigVal[0][i];
   }
   out << "$ExtraTF{DynMatEigVal} = " << setw(Width) << DynMatEigValPrint << "\n";
   out <<  "$ExtraTF{DynMat} = " << setw(Width) << DynMat << "\n";
}

// The following is the Pairwise lattice reduction routine as per Tadmor, Sorkin & Arndt [2009]
// It takes in a bx3 matrix of b rows of Lattice vectors to be reduced
// It returns a bx3 matrix consisting of b row vectors (in R^3) for use with NewCBCellSingleK
Matrix const MultiLatticeKIM::PairwiseReduction(Matrix const& RowLatVects) const
{
   int b, terminate, l, s, m;
   double Norm1, Norm2, DotLS, NormS_squared, value;

   b = RowLatVects.Rows();
   Matrix ReduceLatticeVectorsMatrix(b,3);
   for(int i = 0;i < b;i++)
   {
      for(int j = 0; j< DIM3;j++)
      {
         ReduceLatticeVectorsMatrix[i][j] = RowLatVects[i][j];
      }
   }
   terminate = 0;
   while (terminate == 0)
   {
      terminate = 1;
      for (int i = 0; i < (b-1); i++)
      {
         for (int j = (i+1); j < b; j++)
         {
            Norm1 = sqrt(ReduceLatticeVectorsMatrix[i][0]*ReduceLatticeVectorsMatrix[i][0]+
                         ReduceLatticeVectorsMatrix[i][1]*ReduceLatticeVectorsMatrix[i][1]+
                         ReduceLatticeVectorsMatrix[i][2]*ReduceLatticeVectorsMatrix[i][2]);
            Norm2 = sqrt(ReduceLatticeVectorsMatrix[j][0]*ReduceLatticeVectorsMatrix[j][0]+
                         ReduceLatticeVectorsMatrix[j][1]*ReduceLatticeVectorsMatrix[j][1]+
                         ReduceLatticeVectorsMatrix[j][2]*ReduceLatticeVectorsMatrix[j][2]);
            if (Norm1 >= Norm2)
            {
               l = i;
               s = j;
            }
            else
            {
               l = j;
               s = i;
            }
            DotLS = ReduceLatticeVectorsMatrix[l][0]*ReduceLatticeVectorsMatrix[s][0]+
               ReduceLatticeVectorsMatrix[l][1]*ReduceLatticeVectorsMatrix[s][1]+
               ReduceLatticeVectorsMatrix[l][2]*ReduceLatticeVectorsMatrix[s][2];
            NormS_squared= ReduceLatticeVectorsMatrix[s][0]*ReduceLatticeVectorsMatrix[s][0]+
               ReduceLatticeVectorsMatrix[s][1]*ReduceLatticeVectorsMatrix[s][1]+
               ReduceLatticeVectorsMatrix[s][2]*ReduceLatticeVectorsMatrix[s][2];
            value = DotLS / NormS_squared;
            m = floor( value + 0.5 );
            if (m != 0)
            {
               ReduceLatticeVectorsMatrix[l][0] = ReduceLatticeVectorsMatrix[l][0] - m*ReduceLatticeVectorsMatrix[s][0];
               ReduceLatticeVectorsMatrix[l][1] = ReduceLatticeVectorsMatrix[l][1] - m*ReduceLatticeVectorsMatrix[s][1];
               ReduceLatticeVectorsMatrix[l][2] = ReduceLatticeVectorsMatrix[l][2] - m*ReduceLatticeVectorsMatrix[s][2];
               terminate = 0;
            }
         }
      }
   }
   return ReduceLatticeVectorsMatrix;
}

void MultiLatticeKIM::NewCBCellSingleK(int TFIndex, int Width, ostream& out) const
{
   //Note that MillerIndexCD is a vector in R^5 of the form [h,k,l,c,d]
   int counter;

   Matrix InitLatVects(2,DIM3);
   Matrix ReducedLatVects(2,DIM3);
   Matrix RefLat(DIM3,DIM3);
   Matrix RefLatTemp(DIM3,DIM3);
   Vector G1(DIM3, 0.0);
   Vector G2(DIM3, 0.0);
   Vector G3(DIM3, 0.0);
   Vector G1Star(DIM3, 0.0);
   Vector G2Star(DIM3, 0.0);
   Vector G3Star(DIM3, 0.0);
   Vector G1Plus(DIM3, 0.0);
   Vector G2Plus(DIM3, 0.0);
   Vector G3Plus(DIM3, 0.0);
   Vector GVector(DIM3,0.0);
   Vector K(5, 0.0);
   Vector MinValueVector(DIM3,0.0);

   RefLatTemp = CBK_->RefLattice();

   counter = 0;
   double MinValue;

   int num;
   for(int i = 0;i < DIM3; i++)
   {
      num = 0;
      for (int j = 0;j< DIM3;j++)
      {
         if (abs(RefLatTemp[i][j]) != 0)
         {
            MinValueVector[num] = abs(RefLatTemp[i][j]);
            num++;
         }
      }
      MinValue = MinValueVector[0];
      for (int j = 0;j< num;j++)
      {
         if (MinValueVector[j] < MinValue)
         {
            MinValue = MinValueVector[j];
         }
      }
      for (int j = 0; j<DIM3; j++)
      {
         RefLat[i][j] = RefLatTemp[i][j]/MinValue;
      }
   }

   for (int i = 0; i< DIM3; i++)
   {
      G1[i] = RefLat[0][i];
      G2[i] = RefLat[1][i];
      G3[i] = RefLat[2][i];
   }

   int k1, k2, k3;
   int CommonDivisor;

   int whichTF;
   int whichKV;
   whichTF = TFIndex - (CBK_->DOFS());

   for(int i=0;i<NumKVectors_;++i)
   {
      counter = i*DynMatrixDim_;
      for (int j = counter; j<(counter + DynMatrixDim_);++j)
      {
         if (whichTF == j)
         {
            whichKV = i;
            break;
         }
      }
   }

   for(int j=0;j<5;++j)
   {
      K[j] = KVectorMatrix_[whichKV][j];
   }

   k1 = floor(K[0]+0.5);
   k2 = floor(K[1]+0.5);
   k3 = floor(K[2]+0.5);

   for(int j = 0; j < DIM3; j++)
   {
      InitLatVects[0][j] = 0.0;
      InitLatVects[1][j] = 0.0;
      ReducedLatVects[0][j] = 0.0;
      ReducedLatVects[1][j] = 0.0;
   }

   if (K[2] != 0)
   {
      CommonDivisor = GCD(k2, k3);
      G1Star = (-k3*G2 + k2*G3)/CommonDivisor;
      CommonDivisor = GCD(k1, k3);
      G2Star = (-k3*G1 + k1*G3)/CommonDivisor;
   }
   else if (K[1] != 0)
   {
      CommonDivisor = GCD(k1, k2);
      G1Star = (-k2*G1 + k1*G2)/CommonDivisor;
      G2Star = G3;
   }
   else
   {
      G1Star = G2;
      G2Star = G3;
   }

   for (int j = 0; j < 3; j++)
   {
      InitLatVects[0][j] = G1Star[j];
      InitLatVects[1][j] = G2Star[j];
   }

   ReducedLatVects = PairwiseReduction(InitLatVects);

   for (int j = 0; j < 3; j++)
   {
      G1Plus[j] = ReducedLatVects[0][j];
      G2Plus[j] = ReducedLatVects[1][j];
   }


   //MINIMIZATION ROUTINE
   counter = 0;
   int dValue = 0;
   MinValue = 0;
   for (int i = -K[4]; i <= K[4]; i++)
   {
      for (int j = -K[4]; j <= K[4]; j++)
      {
         for (int k = -K[4]; k <= K[4]; k++)
         {
            dValue = i*K[0] + j*K[1] + k*K[2];
            GVector = i*G1 + j*G2 + k*G3;
            if(dValue == K[4])
            {
               if(counter == 1)
               {
                  if(MinValue > (GVector.Norm()))
                  {
                     MinValue = GVector.Norm();
                     G3Plus = GVector;
                  }
               }
               if(counter == 0)
               {
                  MinValue = GVector.Norm();
                  counter = 1;
                  G3Plus = GVector;
               }
            }
         }
      }
   }
   int MU[DIM3][DIM3];
   for(int j=0;j < DIM3; j++)
   {
      MU[0][j] = floor(G1Plus[j]);
      MU[1][j] = floor(G2Plus[j]);
      MU[2][j] = floor(G3Plus[j]);
   }


   Matrix SuperCellRefLattice(DIM3, DIM3);
   Vector* SuperCellIntPOS = NULL;
   int SuperCellIntAtoms;
   int SuperCellAtomSpecies[CBK_MAX_ATOMS];


   CBK_->SuperCellInfo(MU, SuperCellRefLattice, SuperCellIntAtoms, SuperCellIntPOS,
                       SuperCellAtomSpecies);

   out << "$Lattice{MultiLatticeKIM}{CBKinematics}{LatticeBasis} = [["
       << SuperCellRefLattice[0][0] << ",  " << SuperCellRefLattice[0][1] << ", "
       << SuperCellRefLattice[0][2] << "]," << "\n"
       << "                                                      ["
       << SuperCellRefLattice[1][0] << ",  " << SuperCellRefLattice[1][1] << ", "
       << SuperCellRefLattice[1][2] << "]," << "\n"
       << "                                                      ["
       << SuperCellRefLattice[2][0] << ",  " << SuperCellRefLattice[2][1] << ", "
       << SuperCellRefLattice[2][2] << "]];" << "\n";

   out << "$Lattice{MultiLatticeKIM}{CBKinematics}{InternalAtoms} = "
       << SuperCellIntAtoms << ";" << "\n";

   for (int i = 0; i < SuperCellIntAtoms; i++)
   {

      out << "$Lattice{MultiLatticeKIM}{CBKinematics}{AtomPositions}[" << i << "] = ["
          << SuperCellIntPOS[i][0] << ", " << SuperCellIntPOS[i][1] << ", "
          << SuperCellIntPOS[i][2] << "];" << "\n";
   }


   out << "$Lattice{MultiLatticeKIM}{CBKinematics}{AtomSpecies} = [";
   for (int i = 0; i < SuperCellIntAtoms; i++)
   {
      if (i <(SuperCellIntAtoms - 1))
      {
         out << SuperCellAtomSpecies[i] << ", ";
      }
      else
      {
         out << SuperCellAtomSpecies[i] << "];" << "\n";
      }
   }

   delete[] SuperCellIntPOS;
}

//The following is a function to find the greatest common divisor of two integeres x and y
//If none exists, it might be best to move this to a math object
int MultiLatticeKIM::GCD(int x, int y) const
{
   int t;
   while (y!=0)
   {
      t=y;
      y=x%y;
      x=t;
   }
   return x;
}

Matrix const& MultiLatticeKIM::CondensedModuli() const
{

   Matrix const& stiff = stiffness();
   int intrn = CBK_->Ssize();
   double factor = 1.0 / (intrn / DIM3);
   int fsz = CBK_->Fsize();
   Matrix IM(intrn, intrn);
   CM_static.Resize(fsz, fsz);

   for (int i = 0; i < fsz; i++)
   {
      for (int j = 0; j < fsz; j++)
      {
         CM_static[i][j] = stiff[i][j];
      }
   }

   // Make sure there are internal DOF's
   if (intrn)
   {
      for (int i = 0; i < intrn; i++)
      {
         for (int j = 0; j < intrn; j++)
         {
            IM[i][j] = stiff[fsz + i][fsz + j];

            // add translational stiffness to regularize IM, if needed
            if ((!CBK_->NoTrans()) && (i % DIM3 == j % DIM3))
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
                  CM_static[i][j] -= stiff[i][fsz + m] * IM[m][n] * stiff[fsz + n][j];
               }
            }
         }
      }
   }

   // If using symmetrized F, assume standard Voigt notation
   if (fsz == 6)
   {
      // Remove 2's and 4's
      for (int i = 3; i < 6; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            CM_static[i][j] /= 2.0;
            CM_static[j][i] /= 2.0;
         }

         for (int j = 3; j < 6; j++)
         {
            CM_static[i][j] /= 4.0;
         }
      }
   }

   return CM_static;
}

Vector const& MultiLatticeKIM::ThermalExpansion() const
{
   cerr << "Error in MultiLatticeKIM::ThermalExpansion() empty function \n ";
   exit(-45);
}

int MultiLatticeKIM::comp(void const* const a, void const* const b)
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

int MultiLatticeKIM::abscomp(void const* const a, void const* const b)
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

void MultiLatticeKIM::interpolate(Matrix* const EigVals, int const& zero, int const& one,
                                  int const& two)
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

CMatrix const& MultiLatticeKIM::ReferenceDynamicalStiffness(Vector const& K) const
{
   double pi = 4.0 * atan(1.0);
   MyComplexDouble Ic(0, 1);
   MyComplexDouble A = 2.0 * pi * Ic;
   int i, j;

   Dk_static.Resize(InternalAtoms_ * DIM3, InternalAtoms_ * DIM3, 0.0);

   return Dk_static;
}

void MultiLatticeKIM::ReferenceDispersionCurves(Vector const& K, int const& NoPTS,
                                                char const* const prefix, ostream& out) const
{
   int w = out.width();
   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }

   Matrix InverseLat(DIM3, DIM3);
   InverseLat = (CBK_->RefLattice()).Inverse();

   Matrix EigVal[DIM3];
   for (int i = 0; i < DIM3; ++i)
   {
      EigVal[i].Resize(1, InternalAtoms_ * DIM3);
   }

   Vector Z1(DIM3), Z2(DIM3);
   for (int k = 0; k < DIM3; ++k)
   {
      Z1[k] = K[k];
      Z2[k] = K[DIM3 + k];
   }
   Z1 = InverseLat * Z1;
   Z2 = InverseLat * Z2;

   Vector Z(DIM3),
      DZ = Z2 - Z1;
   double dz = 1.0 / (NoPTS - 1);
   for (int k = 0; k < 2; ++k)
   {
      Z = Z1 + (k * dz) * DZ;
      EigVal[k] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[k][0], InternalAtoms_ * DIM3, sizeof(double), &comp);

      out << prefix << setw(w) << k * dz;
      if (Echo_)
      {
         cout << prefix << setw(w) << k * dz;
      }
      for (int i = 0; i < InternalAtoms_ * DIM3; ++i)
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
      Z = Z1 + (k * dz) * DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[two][0], InternalAtoms_ * DIM3, sizeof(double), &comp);
      interpolate(EigVal, zero, one, two);

      out << prefix << setw(w) << k * dz;
      if (Echo_)
      {
         cout << prefix << setw(w) << k * dz;
      }
      for (int i = 0; i < InternalAtoms_ * DIM3; ++i)
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

int MultiLatticeKIM::ReferenceBlochWave(Vector& K) const
{
   InverseLat_static = (CBK_->RefLattice()).Inverse();

   // Iterate over points in cubic unit cell
   for (UCIter_.Reset(); !UCIter_.Done(); ++UCIter_)
   {
      for (int i = 0; i < DIM3; ++i)
      {
         K[i] = UCIter_[i];
      }

      Z_static = InverseLat_static * K;
      A_static = ReferenceDynamicalStiffness(Z_static);

      EigVals_static = HermiteEigVal(A_static);

      for (int i = 0; i < InternalAtoms_ * DIM3; ++i)
      {
         // if w^2 <= 0.0 --> Re(i*w*x) > 0 --> growing solutions --> unstable
         if (EigVals_static[0][i] <= 0.0)
         {
            return 0;
         }
      }
   }
   return 1;
}

// ---- needs to be updated---- //
void MultiLatticeKIM::LongWavelengthModuli(double const& dk, int const& gridsize,
                                           char const* const prefix, ostream& out) const
{
   cout << "ERROR IN MultiLatticeKIM::LongWavelengthModuli. Exiting" << endl;
   exit(1);
}

void MultiLatticeKIM::NeighborDistances(int const& cutoff, ostream& out) const
{
   cout << "ERROR IN MultiLatticeKIM::NeighborDistances. Exiting" << endl;
   exit(1);
}

void MultiLatticeKIM::Print(ostream& out, PrintDetail const& flag,
                            PrintPathSolutionType const& SolType)
{

   int W;
   int NoNegTestFunctions = 0;
   double engy, entropy, heatcapacity;
   int RankOneConvex;
   int BlochWaveStable;
   double mintestfunct;
   double conj;
   int NoFP = !FastPrint_;

   cout << "MultiLatticeKIM::Print::FastPrint_ = " << FastPrint_ << endl;
   cout << "MultiLatticeKIM::Print::NoFP =" << NoFP << endl;

   W = out.width();

   out.width(0);
   if (Echo_)
   {
      cout.width(0);
   }
   LatSum_.Recalc();
   LatSum_.Reset();

   engy = E0();
   //conj = ConjugateToLambda();
   str_static = stress();
   cout << "::Print::str_static = " << setw(15) << str_static << endl;
   cout << "::Print::E1() = " << setw(15) << E1() << endl;

   if (NoFP)
   {
      stiff_static = stiffness();
//      TE_static = ThermalExpansion();
//      entropy = Entropy();
//      heatcapacity = HeatCapacity();

      TestFunctions(TestFunctVals_static, LHS);
      mintestfunct = TestFunctVals_static[0];

      if(TFType_ == 2) // LoadingParameters
      {
         for(int i = 0; i < CBK_->DOFS(); ++i)
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
      else // KVectors or None
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

      CondModuli_static = CondensedModuli();
      CondEV_static = SymEigVal(CondModuli_static);
      RankOneConvex = FullScanRank1Convex3D(CBK_, CondModuli_static, ConvexityDX_);

      K_static.Resize(DIM3, 0.0);
      if (RankOneConvex)
      {
         BlochWaveStable = BlochWave(K_static);
      }
      else
      {
         BlochWaveStable = -1;
      }
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
                << "          Position : " << setw(W) << CBK_->AtomPositions(i) << "\n";
         }
         out << "REFTemp_ : " << setw(W) << REFTemp_ << "\n";
         out << "REFLambda_ : " << setw(W) << REFLambda_ << "\n";
         out << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
         out << "Tref = " << setw(W) << Tref_ << "\n";
         out << "Normalization Modulus : " << setw(W) << NormModulus_ << "\n";
         out << "EulerAngles : " << setw(W) << EulerAng_[0]
             << setw(W) << EulerAng_[1] << setw(W) << EulerAng_[2] << "\n";
         out << "Loading Proportions : " << setw(W) << LoadingProportions_ << "\n";

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
                    << "          Position : " << setw(W) << CBK_->AtomPositions(i) << "\n";
            }
            cout << "REFTemp_ : " << setw(W) << REFTemp_ << "\n";
            cout << "REFLambda_ : " << setw(W) << REFLambda_ << "\n";
            cout << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
            cout << "Tref = " << setw(W) << Tref_ << "\n";
            cout << "Normalization Modulus : " << setw(W) << NormModulus_ << "\n";
            cout << "EulerAngles : " << setw(W) << EulerAng_[0]
                 << setw(W) << EulerAng_[1] << setw(W) << EulerAng_[2] << "\n";
            cout << "Loading Proportions : " << setw(W) << LoadingProportions_ << "\n";
         }
         // passthrough to short

      case PrintShort:
         out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << "\n"
             << "Lambda (Normalized): " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << CBK_->DOF() << "\n"
             << "Potential Value (Normalized):" << setw(W) << engy << "\n";
         for (int i = 0; i < InternalAtoms_; ++i)
         {
            out << "BodyForce Value (KIM Output) " << i << setw(W) << BodyForce_[i] << "\n";
         }
         out << "Stress (Normalized):" << "\n" << setw(W) << str_static << "\n";
         out << "\nStiffness (Normalized):" << setw(W) << stiff_static
             << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n" << setw(W)
             << TestFunctVals_Print << "\n"
             << "Bifurcation Info:" << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n"
             << "Condensed Moduli (Normalized):" << setw(W) << CondModuli_static
             << "CondEV Info:" << setw(W) << CondEV_static
             << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
             << "BlochWave Stability (GridSize=" << GridSize_ << "):"
             << setw(W) << BlochWaveStable << ", "
             << setw(W) << K_static << endl;


         if (Echo_)
         {
            cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << "\n"
                 << "Lambda (Normalized): " << setw(W) << Lambda_ << "\n"
                 << "ConjugateToLambda: " << setw(W) << conj << "\n"
                 << "DOF's :" << "\n" << setw(W) << CBK_->DOF() << "\n"
                 << "Potential Value (Normalized):" << setw(W) << engy << "\n";
            for (int i = 0; i < InternalAtoms_; ++i)
            {
               cout << "BodyForce Value (KIM Output) " << i << setw(W) << BodyForce_[i] << "\n";
            }
            cout << "Stress (Normalized):" << "\n" << setw(W) << str_static << "\n";
            cout << "\nStiffness (Normalized):" << setw(W) << stiff_static
                 << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n" << setw(W)
                 << TestFunctVals_static << "\n"
                 << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n"
                 << "Condensed Moduli (Normalized):" << setw(W) << CondModuli_static
                 << "CondEV Info:" << setw(W) << CondEV_static
                 << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
                 << "BlochWave Stability (GridSize=" << GridSize_ << "):"
                 << setw(W) << BlochWaveStable << ", "
                 << setw(W) << K_static << endl;
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

ostream& operator<<(ostream& out, MultiLatticeKIM& A)
{
   A.Print(out, Lattice::PrintShort);
   return out;
}


// ---------------------- Debug Mode Handler --------------------------


void MultiLatticeKIM::DebugMode()
{
   const char* Commands[] = {
      "InternalAtoms_",
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
//      "Potential_",
      "ConvexityDX_",
      "ConjugateToLambda",
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
      "EulerAng_",
      "Rotation_",
      "Loading_",
      "PrintCrystal",
      "ThermalExpansion",
      "Entropy",
      "HeatCapacity",
      "TranslationProjection1D",
      "TranslationProjection3D",
   };
   int NOcommands = 53;

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
         cout << "InternalAtoms_ = " << InternalAtoms_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "CBK_->DOFS() = " << CBK_->DOFS() << "\n";
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
         for (int i = 0; i < CBK_->DOFS(); ++i)
         {
            cout << "DOF_[" << i << "] = " << (CBK_->DOF())[i] << "\n";
         }
      }
      else if (response == Commands[indx++])
      {
         cout << "RefLattice_= " << setw(W) << CBK_->RefLattice();
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
         for (int i = 0; i < InternalAtoms_; ++i)
         {
            cout << "BodyForce_[" << i << "]= " << setw(W)
                 << BodyForce_[i] << "\n";
         }
      }
      else if (response == Commands[indx++])
      {
         for (int i = 0; i < InternalAtoms_; ++i)
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
         for (int i = 0; i < InternalAtoms_; ++i)
         {
            for (int j = i; j < InternalAtoms_; ++j)
            {
//               cout << "Potential_[" << i << "][" << j << "]= "
               //                  << setw(W) << Potential_[i][j] << "\n";
               cout << "ERROR. NO POTENTIAL NECESSARY" << endl;
            }
         }
      }
      else if (response == Commands[indx++])
      {
         cout << "ConvexityDX_= " << ConvexityDX_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "ConjugateToLambda = " << setw(W) << ConjugateToLambda();
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
         Vector K(DIM3, 0.0);
         int NoPTS;
         string prefix;
         int oldEcho_ = Echo_;
         cout << "\tK > ";
         cin >> K;
         cout << "\tNoPTS > ";
         cin >> NoPTS;
         cout << "\tprefix > ";
         cin >> prefix;
         Echo_ = 0;
         cout << "ReferenceDispersionCurves= ";
         ReferenceDispersionCurves(K, NoPTS, prefix.c_str(), cout);
         Echo_ = oldEcho_;
      }
      else if (response == Commands[indx++])
      {
         Vector K(DIM3, 0.0);
         cout << "ReferenceBlochWave= " << ReferenceBlochWave(K) << "\t" << K << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "\tK > ";
         Vector K(DIM3, 0.0);
         cin >> K;
         cout << "ReferenceDynamicalStiffness= "
              << setw(W) << ReferenceDynamicalStiffness(K) << "\n";
      }
      else if (response == Commands[indx++])
      {
         Vector DOF(CBK_->DOFS(), 0.0);
         cout << "\tDOF > ";
         cin >> DOF;
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
         SetTemp(Temp);
      }
      else if (response == Commands[indx++])
      {
         double dist;
         cout << "\tInfluenceDist > ";
         cin >> dist;
         SetInfluenceDist(dist);
      }
      else if (response == Commands[indx++])
      {
         cout << "energy= " << E0() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "E0= " << E0() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "E1= " << setw(W) << E1() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "E2= " << setw(W) << E2() << "\n";
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
         SetGridSize(GridSize);
      }
      else if (response == Commands[indx++])
      {
         int oldEcho_ = Echo_;
         int cutoff;
         cout << "\tcutoff > ";
         cin >> cutoff;
         Echo_ = 0;
         cout << setw(W);
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
         cout << "EulerAng_ = "
              << setw(W) << EulerAng_[0]
              << setw(W) << EulerAng_[1]
              << setw(W) << EulerAng_[2] << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "Rotation_ = "
              << setw(W) << Rotation_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "Loading_ = "
              << setw(W) << Loading_ << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << setw(W);
         PrintCurrentCrystalParamaters(cout);
      }
      else if (response == Commands[indx++])
      {
         cout << "ThermalExpansion = " << setw(W) << ThermalExpansion() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "Entropy = " << setw(W) << Entropy() << "\n";
      }
      else if (response == Commands[indx++])
      {
         cout << "HeatCapacity = " << setw(W) << HeatCapacity() << "\n";
      }
      else if (response == Commands[indx++])
      {
         int n;
         cout << "\tNoAtoms > ";
         cin >> n;
         cout << "P.Transpose()\n"
              << setw(20) << TranslationProjection1D(n).Transpose() << "\n";
      }
      else if (response == Commands[indx++])
      {
         int n, f;
         cout << "\tFsize > ";
         cin >> f;
         cout << "\tNoAtoms > ";
         cin >> n;
         cout << "P.Transpose()\n"
              << setw(20) << TranslationProjection3D(f, n).Transpose() << "\n";
      }
      else if (response == Commands[indx++])
      {
         /*
           int no = SpeciesPotential_[0][0]->GetNoParameters();
           double* vals;
           int sze = no * ((CBK_->NumberofSpecies() + 1) * (CBK_->NumberofSpecies()) / 2);
           vals = new double[sze];
           cout << "\tEnter new Parameter values > ";
           for (int i = 0; i < sze; ++i)
           {
           cin >> vals[i];
           }
           SetParameters(vals);
         */
         cout << "ERROR.  POTENTIAL NOT APPLICABLE IN KIM" << endl;
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


void MultiLatticeKIM::RefineEqbm(double const& Tol, int const& MaxItr, ostream* const out)
{
   Vector dx(CBK_->DOFS(), 0.0);
   cout << "RefEqbm::CBK_->DOF()::1 = " << setw(15)<< CBK_->DOF() << endl;

   Vector Stress = E1();
   cout << "RefEqbm::E1() = " << setw(15) << Stress << endl;
   int itr = 0;
   while ((itr < MaxItr) && Stress.Norm() > Tol)
   {
      ++itr;
      cout << "E2().Det() = " << setw(15) << E2().Det() << endl;
      cout << "dx = " << setw(15) << dx << endl;

#ifdef SOLVE_SVD
      dx = SolveSVD(E2(), Stress, MAXCONDITION, Echo_);
#else
      dx = SolvePLU(E2(), Stress);
#endif

      break;

      SetDOF(CBK_->DOF() - dx);
      Stress = E1();
      if (out != 0)
      {
         *out << setw(20) << Stress;

         *out << itr << "\tdx " << dx.Norm() << "\tstress " << Stress.Norm() << "\n";
      }
   }
   cout << "RefEqbm::CBK_->DOF()::2 = " << setw(15) << CBK_->DOF() << endl;
   cout << "End of RefineEqbm " << endl;
}

void MultiLatticeKIM::PrintCurrentCrystalParamaters(ostream& out) const
{
   Matrix F(DIM3, DIM3, 0.0);
   Vector CurrentLattice[DIM3];
   int W = out.width();
   out.width(0);

   for (int i = 0; i < DIM3; ++i)
   {
      CurrentLattice[i].Resize(DIM3);
      CurrentLattice[i] = CBK_->CurrentLatticeVec(i);
   }


   out << "TITLE LatticeStatics crystal structure scaled by 10.0" << "\n";
   out << "DIMENSION 3" << "\n";
   out << "CELL" << setw(W) << 10.0 * CurrentLattice[0].Norm()
       << setw(W) << 10.0 * CurrentLattice[1].Norm()
       << setw(W) << 10.0 * CurrentLattice[2].Norm();

   double alpha, beta, gamma;
   double pi = 4.0 * atan(1.0);

   alpha = acos(CurrentLattice[1] * CurrentLattice[2]
                / (CurrentLattice[1].Norm() * CurrentLattice[2].Norm()));
   beta = acos(CurrentLattice[0] * CurrentLattice[2]
               / (CurrentLattice[0].Norm() * CurrentLattice[2].Norm()));
   gamma = acos(CurrentLattice[0] * CurrentLattice[1]
                / (CurrentLattice[0].Norm() * CurrentLattice[1].Norm()));

   out << setw(W) << alpha * 180.0 / pi
       << setw(W) << beta * 180.0 / pi
       << setw(W) << gamma * 180.0 / pi << "\n";

   out << "SYMMETRY  NUMBER 1  LABEL P1  " << "\n";
   out << "SYM MAT  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0 0.0000 0.0000 0.0000" << "\n";
   out << "\n" << "ATOMS" << "\n"
       << "NAME" << setw(W) << "X" << setw(W) << "Y" << setw(W) << "Z" << "\n";
   char const* species[] = {"Ni", "Ti", "C"};
   out << setw(4) << species[(CBK_->AtomSpecies(0) > 3) ? 3 : CBK_->AtomSpecies(0)]
       << setw(W) << CBK_->FractionalPosVec(0) << "\n";
   for (int i = 1; i < InternalAtoms_; ++i)
   {
      out << setw(4) << species[(CBK_->AtomSpecies(i) > 3) ? 3 : CBK_->AtomSpecies(i)];
      out << setw(W) << CBK_->FractionalPosVec(i) << "\n";
   }

   out << "EOF" << "\n";

   out << "\n"
       << "Temperature : " << setw(W) << NTemp_ << "\n"
       << "Lambda : " << setw(W) << Lambda_ << "\n"
       << "DOFs : " << setw(W) << CBK_->DOF() << "\n"
       << "Lattice Vectors : " << "\n"
       << setw(W) << CurrentLattice[0] << "\n"
       << setw(W) << CurrentLattice[1] << "\n"
       << setw(W) << CurrentLattice[2] << "\n";
}

//CHECK OUT SSTREAM STRING INCLUDES IN QC. Line 527 in QC.cpp, Line 598 to return pointer to a char of string object
void MultiLatticeKIM::Write_KIM_descriptor_file(const char** SpeciesList, int numberParticleTypes_)
{
   descriptor_file_ << "########################################################################################################" << endl;
   descriptor_file_ << "#" << endl;
   descriptor_file_ << "# Copyright 2011 Ellad B. Tadmor, Ryan S. Elliott, and James P. Sethna" << endl;
   descriptor_file_ << "# All rights reserved. " << endl;
   descriptor_file_ << "# " << endl;
   descriptor_file_ << "# Author: Automatically generated by BFB" << endl;
   descriptor_file_ << "#" << endl;
   descriptor_file_ << "#" << endl;
   descriptor_file_ << "#######################################################################################################" << endl;
   descriptor_file_ << "TEST_NAME := BFB" << endl;
   descriptor_file_ << "Unit_length      := A  " << endl;
   descriptor_file_ << "Unit_energy      := eV " << endl;
   descriptor_file_ << "Unit_charge      := e " << endl;
   descriptor_file_ << "Unit_temperature := K " << endl;
   descriptor_file_ << "Unit_time        := ps     " << endl;
   descriptor_file_ << "#######################################################################################################" << endl;
   descriptor_file_ << "SUPPORTED_ATOM/PARTICLES_TYPES:                                                              " << endl;

   for (int i = 0; i <numberParticleTypes_;i++)
   {
      descriptor_file_ << SpeciesList[i] << "                          spec                    1" << endl;
   }

   descriptor_file_ << "#######################################################################################################" << endl;
   descriptor_file_ << "CONVENTIONS:" << endl;
   descriptor_file_ << "# Name                      Type" << endl;
   descriptor_file_ << "ZeroBasedLists              flag" << endl;
   descriptor_file_ << "Neigh_BothAccess            flag" << endl;
//      descriptor_file_ << "CLUSTER                    flag" << endl;
   descriptor_file_ << "NEIGH_RVEC_F                           flag" << endl;
   descriptor_file_ << "#######################################################################################################" << endl;
   descriptor_file_ << "MODEL_INPUT:" << endl;
   descriptor_file_ << "# Name                      Type         Unit                Shape              Requirements" << endl;
   descriptor_file_ << "numberOfParticles           integer      none                []" << endl;
   descriptor_file_ << "numberParticleTypes         integer      none                []" << endl;
   descriptor_file_ << "particleTypes               integer      none                [numberOfParticles]" << endl;
   descriptor_file_ << "coordinates                 real*8       length              [numberOfParticles,3]" << endl;
   descriptor_file_ << "get_neigh                   method       none                []" << endl;
   descriptor_file_ << "neighObject                 pointer      none                []" << endl;
   descriptor_file_ << "process_dEdr                method       none                []" << endl;
   descriptor_file_ << "process_d2Edr2              method       none                []" << endl;
   descriptor_file_ << "#######################################################################################################" << endl;
   descriptor_file_ << "MODEL_OUTPUT:" << endl;
   descriptor_file_ << "# Name                      Type         Unit                Shape              Requirements" << endl;
   descriptor_file_ << "destroy                     method       none                []" << endl;
   descriptor_file_ << "compute                     method       none                []" << endl;
   descriptor_file_ << "cutoff                      real*8       length              []" << endl;
   descriptor_file_ << "energy                      real*8       energy              []" << endl;
   descriptor_file_ << "forces                      real*8       force               [numberOfParticles,3]" << endl;
}


int MultiLatticeKIM::get_neigh(void* kimmdl, int *mode, int *request, int* atom,
                               int* numnei, int** nei1atom, double** Rij)
{
   intptr_t* pkim = *((intptr_t**) kimmdl);
   int atomToReturn;
   int status;

   PPSumKIM* LatSum;
   LatSum = (PPSumKIM*) KIM_API_get_data(pkim, "neighObject", &status);

   int* numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);

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
         KIM_API_report_error(__LINE__, (char*) __FILE__, (char*) "Invalid request in get_periodic_neigh", KIM_STATUS_NEIGH_INVALID_REQUEST);
         return KIM_STATUS_NEIGH_INVALID_REQUEST;
      }

   }
   else if (1 == *mode) // locator mode
   {
      if ((*request >= *numberOfParticles) || (*request < 0)) // invalid id
      {
         KIM_API_report_error(__LINE__,  (char*) __FILE__, (char*) "Invalid atom ID in get_neigh", KIM_STATUS_PARTICLE_INVALID_ID);
         return KIM_STATUS_PARTICLE_INVALID_ID;
      }
      else
      {
         LatSum->Reset();
         for (int i = 0; i < *request;i++)
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
      KIM_API_report_error(__LINE__, (char*)__FILE__, (char*)"Invalid mode in get_periodic_neigh", KIM_STATUS_NEIGH_INVALID_MODE);
      return KIM_STATUS_NEIGH_INVALID_MODE;
   }

}


int MultiLatticeKIM::process_dEdr(void *kimmdl, double *dEdr, double *r, double **dx, int *i, int *j)
{

   int status;
   intptr_t* pkim = *((intptr_t**) kimmdl);
   MultiLatticeKIM* obj;
   obj = (MultiLatticeKIM*) KIM_API_get_test_buffer(pkim, &status);

   double DX[3];

   Matrix InverseF = (obj->CBK_->F()).Inverse();

   double temp1, temp2;
   for(int rows = 0; rows < 3; rows++)
   {
      DX[rows] = 0.0;
      for(int cols = 0; cols < 3; cols++)
      {
         DX[rows] += InverseF[rows][cols] * (*dx)[cols];
      }
   }

   //dEdUij
   for(int rows = 0; rows < 3; rows++)
   {
      for(int cols = 0; cols < 3; cols++)
      {
         obj->ME1_static[(obj->CBK_->INDF(rows,cols))] += *dEdr * (obj->CBK_->DyDF(*dx, DX, rows, cols))/(2.0*(*r));
      }
   }

   //dEdSij
   for(int atom = 0; atom < (obj->CBK_->InternalAtoms()); atom++)
   {
      for(int k = 0; k < 3; k++)
      {
         obj->ME1_static[(obj->CBK_->INDS(atom,k))] += *dEdr * (obj->CBK_->DyDS(*dx , *i,*j, atom, k))/(2.0*(*r));
      }
   }

   if ((obj->StiffnessYes_) == 1)
   {
      double DyDF[3][3];
      int i1, j1, k1, l1;

      // Upper Diag Block (CBK_->Fsize(),CBK_->Fsize())
      for (i1 = 0; i1 < 3; i1++)
      {
         for (j1 = 0; j1 < 3; j1++)
         {
            for (k1 = 0; k1 < 3; k1++)
            {
               for (l1 = 0; l1 < 3; l1++)
               {
                  obj->ME2_static[obj->CBK_->INDF(i1,j1)][obj->CBK_->INDF(k1,l1)] +=  (*dEdr) *((0.5/ (*r)) * (obj->CBK_->D2yDFF(DX, i1, j1, k1, l1)) - (0.25/ ((*r)*(*r)*(*r))) * (obj->CBK_->DyDF((*dx), DX, i1, j1)) * (obj->CBK_->DyDF((*dx), DX, k1, l1)));
               }
            }
         }
      }

      //Lower Diagonal blocks
      for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
      {
         for (int atom1 = 0; atom1 < (obj->CBK_->InternalAtoms()); atom1++)
         {
            for (k1 = 0; k1 < 3; k1++)
            {
               for (l1 = 0; l1 < 3; l1++)
               {

                  obj->ME2_static[obj->CBK_->INDS(atom0,k1)][obj->CBK_->INDS(atom1,l1)] += (*dEdr) * ((0.5/ (*r)) * obj->CBK_->D2yDSS(*i, *j, atom0, k1, atom1, l1)  - (0.25/ ((*r)*(*r)*(*r))) * (obj->CBK_->DyDS(*dx , *i,*j, atom0, k1)) * (obj->CBK_->DyDS(*dx , *i,*j, atom1, l1)));

               }
            }
         }

      }

      //Off-diagonal blocks
      for (i1 = 0; i1 < 3; i1++)
      {
         for (j1 = 0; j1 < 3; j1++)
         {
            for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
            {
               for (k1 = 0; k1 < 3; k1++)
               {
                  obj->ME2_static[obj->CBK_->INDF(i1,j1)][obj->CBK_->INDS(atom0,k1)] += (*dEdr) * ((0.5/ (*r)) * obj->CBK_->D2yDFS(*dx, DX, *i, *j, i1, j1, atom0, k1) - (0.25/ ((*r)*(*r)*(*r))) * (obj->CBK_->DyDS(*dx , *i,*j, atom0, k1)) * (obj->CBK_->DyDF((*dx), DX, i1, j1)));

                  obj->ME2_static[obj->CBK_->INDS(atom0,k1)][obj->CBK_->INDF(i1,j1)] += (*dEdr) * ((0.5/ (*r)) * obj->CBK_->D2yDFS(*dx, DX, *i, *j, i1, j1, atom0, k1) - (0.25/ ((*r)*(*r)*(*r))) * (obj->CBK_->DyDS(*dx , *i,*j, atom0, k1)) * (obj->CBK_->DyDF((*dx), DX, i1, j1)));

               }

            }
         }
      }

   }
   return KIM_STATUS_OK;
}


int MultiLatticeKIM::process_d2Edr2(void *kimmdl, double *d2Edr2, double **r, double **dx, int **i, int **j)
{
   int status;
   intptr_t* pkim = *((intptr_t**) kimmdl);
   MultiLatticeKIM* obj;
   obj = (MultiLatticeKIM*) KIM_API_get_test_buffer(pkim, &status);

   double DX[2][3];
   double Dx[2][3];

   Matrix InverseF = (obj->CBK_->F()).Inverse();

   for (int k = 0; k < 3; k++)
   {
      Dx[0][k] = (*dx)[k];
      Dx[1][k] = (*dx)[k+3];

   }

   for (int atoms = 0; atoms < 2; atoms++)
   {
      for(int rows = 0; rows < 3; rows++)
      {
         DX[atoms][rows] = 0.0;
         for(int cols = 0; cols < 3; cols++)
         {
            DX[atoms][rows] += InverseF[rows][cols] * Dx[atoms][cols];
         }
      }
   }
   int i1, j1, k1, l1;

   //Upper Diagonal
   for (i1 = 0; i1 < 3; i1++)
   {
      for (j1 = 0; j1 < 3; j1++)
      {
         for (k1 = 0; k1 < 3; k1++)
         {
            for (l1 = 0; l1 < 3; l1++)
            {
               obj->ME2_static[obj->CBK_->INDF(i1,j1)][obj->CBK_->INDF(k1,l1)] += (0.5/((*r)[0])) * (0.5/((*r)[1])) * (*d2Edr2) * (obj->CBK_->DyDF(Dx[0], DX[0], i1, j1)) * (obj->CBK_->DyDF(Dx[1], DX[1], k1, l1));
            }
         }
      }
   }

   //Lower Diagonal blocks
   for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
   {
      if(((atom0 == (*i)[0]) || (atom0 == (*j)[0])))
      {
         for (int atom1 = 0; atom1 < (obj->CBK_->InternalAtoms()); atom1++)
         {
            if(((atom1 == (*i)[1]) || (atom1 == (*j)[1])))
            {
               for (k1 = 0; k1 < 3; k1++)
               {
                  for (l1 = 0; l1 < 3; l1++)
                  {
                     obj->ME2_static[obj->CBK_->INDS(atom0,k1)][obj->CBK_->INDS(atom1,l1)] += (0.5/((*r)[0])) * (0.5/((*r)[1])) * (*d2Edr2) * (obj->CBK_->DyDS(Dx[0] , (*i)[0],(*j)[0], atom0, k1)) * (obj->CBK_->DyDS(Dx[1] , (*i)[1], (*j)[1], atom1, l1));
                  }
               }
            }
         }
      }
   }

   //Off Diagonal
   int RowIndex;
   for (int atom0 = 0; atom0 < (obj->CBK_->InternalAtoms()); atom0++)
   {
      for (int k = 0; k < 2; k++)
      {
         if ((atom0 == (*i)[k] || atom0 == (*j)[k]))
         {
            RowIndex = k;
            k = 2;

            for (i1 = 0; i1 < 3; i1++)
            {
               for (j1 = 0; j1 < 3; j1++)
               {
                  for (k1 = 0; k1 < 3; k1++)
                  {
                     obj->ME2_static[obj->CBK_->INDF(i1,j1)][obj->CBK_->INDS(atom0,k1)] += (0.5/((*r)[RowIndex])) * (0.5/((*r)[RowIndex])) * (*d2Edr2) * obj->CBK_->DyDF(Dx[RowIndex], DX[RowIndex], i1, j1) * obj->CBK_->DyDS(Dx[RowIndex] , (*i)[RowIndex], (*j)[RowIndex], atom0, k1);

                     obj->ME2_static[obj->CBK_->INDS(atom0,k1)][obj->CBK_->INDF(i1,j1)] += (0.5/((*r)[RowIndex])) * (0.5/((*r)[RowIndex])) * (*d2Edr2) * obj->CBK_->DyDF(Dx[RowIndex], DX[RowIndex], i1, j1) * obj->CBK_->DyDS(Dx[RowIndex] , (*i)[RowIndex], (*j)[RowIndex], atom0, k1);
                  }
               }
            }
         }
      }
   }

   return KIM_STATUS_OK;
}


