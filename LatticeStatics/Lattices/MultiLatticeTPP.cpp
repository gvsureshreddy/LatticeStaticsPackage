#include "MultiLatticeTPP.h"
#include "UtilityFunctions.h"
#include <cmath>

using namespace std;

int const MultiLatticeTPP::DIM3 = 3;

double const RoEig_[3] = {1.0,2.0,3.0};
double const TrEig_[3] = {4.0,5.0,6.0};

MultiLatticeTPP::~MultiLatticeTPP()
{
   delete [] EulerAng_;
   delete [] BodyForce_;
   delete [] SpeciesMass_;
   delete [] AtomicMass_;
   for (int i=0;i<NumberofSpecies_;++i)
      for (int j=i;j<NumberofSpecies_;++j)
         delete SpeciesPotential_[i][j];
   delete [] SpeciesPotential_[0];
   delete [] SpeciesPotential_;
   delete [] Potential_[0];
   delete [] Potential_;
   delete CBK_;
}

MultiLatticeTPP::MultiLatticeTPP(PerlInput const& Input,int const& Echo,int const& Width,
                                 int const& Debug)
   : Lattice(Input)
{
   Echo_ = Echo;
   dbg_ = Debug;
   // Get Lattice definition
   PerlInput::HashStruct Hash = Input.getHash("Lattice","MultiLatticeTPP");
   
   // Set default values
   KillTranslations_ = 1; // 1-true, 0-false
   int needKillRotations = 1;
   KillRotations_ = 0; // 0-do nothing, 1-kill one rotation, 2-kill three rotations

   PerlInput::HashStruct CBKHash = Input.getHash(Hash,"CBKinematics");
   const char *CBKin = Input.getString(CBKHash,"Type");
   if (!strcmp("SymLagrangeCB",CBKin))
   {
      CBK_ = new SymLagrangeCB(Input,&Hash);
      KillTranslations_ = 0;
      needKillRotations = 0;
   }
   else if (!strcmp("SymLagrangeWTransCB",CBKin))
   {
      CBK_ = new SymLagrangeWTransCB(Input,&Hash);
      needKillRotations = 0;
   }
   else if (!strcmp("LagrangeCB",CBKin))
   {
      CBK_ = new LagrangeCB(Input,&Hash);
   }
   else if (!strcmp("MixedCB",CBKin))
   {
      CBK_ = new MixedCB(Input,&Hash);
   }
   else if (!strcmp("EulerCB",CBKin))
   {
      CBK_ = new EulerCB(Input,&Hash);
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
      Vector R(DIM3),r(DIM3);
      double norm;
      const char *KillRot = Input.getString(Hash,"RotationConstraint",0);
      if (!strcmp("FullRotationConstraint",KillRot))
      {
         KillRotations_=2;
      }
      else if (!strcmp("OneRotationConstraint",KillRot))
      {
         KillRotations_=1;
         Input.getVector(r,Hash,"RotationConstraint",1);
         Input.getVector(R,Hash,"RotationConstraint",2);
         KillOneRotation_.Resize(CBK_->DOFS(),0.0);
         for (int i=0;i<DIM3;++i)
            for (int j=0;j<DIM3;++j)
            {
               KillOneRotation_[CBK_->INDF(i,j)] = r[i]*R[j];
            }
         norm=KillOneRotation_.Norm();
         for (int i=0;i<KillOneRotation_.Dim();++i)
            KillOneRotation_[i] /= norm;
      }
      else if (!strcmp("NoRotationConstraint",KillRot))
      {
            KillRotations_=0;
      }
      else
      {
         cerr << "Error (MultiLatticeTPP()): Unknown RotationConstraint type" << "\n";
         exit(-2);
      }
   }
   
   if (DOFMAX < CBK_->DOFS())
   {
      cerr << "Error (MultiLatticeTPP()): DOFMAX < " << CBK_->DOFS() << " in Lattice.h" << "\n";
      exit(-5);
   }
   
   // Setup Bodyforce_
   BodyForce_ = new Vector[InternalAtoms_];
   for (int i=0;i<InternalAtoms_;++i)
      BodyForce_[i].Resize(DIM3,0.0);
   
   // Get Thermo parameters
   Tref_ = Input.getDouble(Hash,"Tref");
   //PhiRef_ = Input.getDouble(Hash,"PhiRef");
   //EntropyRef_ = Input.getDouble(Hash,"EntropyRef");
   //HeatCapacityRef_ = Input.getDouble(Hash,"HeatCapacityRef");

   Input.getPosIntVector(AtomSpecies_,InternalAtoms_,Hash,"AtomSpecies");
   NumberofSpecies_ = AtomSpecies_[0];
   for (int i=1;i<InternalAtoms_;++i)
      if (NumberofSpecies_ < AtomSpecies_[i])
         NumberofSpecies_ = AtomSpecies_[i];
   NumberofSpecies_++;
   
   // Get Potential Parameters
   SpeciesPotential_ = new PairPotentials**[NumberofSpecies_];
   SpeciesPotential_[0] = new PairPotentials*[NumberofSpecies_*NumberofSpecies_];
   for (int i=1;i<NumberofSpecies_;++i)
   {
      SpeciesPotential_[i] = SpeciesPotential_[i-1] + NumberofSpecies_;
   }
   Potential_ = new PairPotentials**[InternalAtoms_];
   Potential_[0] = new PairPotentials*[InternalAtoms_*InternalAtoms_];
   for (int i=1;i<InternalAtoms_;++i)
   {
      Potential_[i] = Potential_[i-1] + InternalAtoms_;
   }
   
   SpeciesMass_ = new double[NumberofSpecies_];
   AtomicMass_ = new double[InternalAtoms_];
   
   for (int i=0;i<NumberofSpecies_;++i)
   {
      for (int j=i;j<NumberofSpecies_;++j)
      {
         SpeciesPotential_[i][j] = SpeciesPotential_[j][i]
            = InitializePairPotential(Hash,Input,i,j);
      }
      SpeciesMass_[i] = Input.getDouble(Hash,"AtomicMasses",i);
   }
   
   for (int i=0;i<InternalAtoms_;++i)
   {
      for (int j=i;j<InternalAtoms_;++j)
      {
         Potential_[i][j] = Potential_[j][i]
            = SpeciesPotential_[AtomSpecies_[i]][AtomSpecies_[j]];
      }
      
      AtomicMass_[i] = SpeciesMass_[AtomSpecies_[i]];
   }
   
   // Get Lattice parameters
   NTemp_ = 1.0;
   InfluenceDist_ = Input.getDouble(Hash,"InfluenceDist");
   if (Input.ParameterOK(Hash,"Density"))
   {
      Density_ = Input.getInt(Hash,"Density");
   }
   else
   {
      Density_ = Input.useInt(1,Hash,"Density"); // Default Value
   }
   NormModulus_ = Input.getDouble(Hash,"NormModulus");
   ConvexityDX_ = Input.getDouble(Hash,"ConvexityDX");
   
   // Get Loading parameters
   const char *loadparam = Input.getString(Hash,"LoadingParameter");
   if (!strcmp("Temperature",loadparam))
   {
      LoadParameter_ = Temperature;
   }
   else if (!strcmp("Load",loadparam))
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
   EulerAng_[0] = Input.getDouble(Hash,"EulerAngle_X");
   EulerAng_[1] = Input.getDouble(Hash,"EulerAngle_Y");
   EulerAng_[2] = Input.getDouble(Hash,"EulerAngle_Z");
   LoadingProportions_.Resize(DIM3);
   Input.getVector(LoadingProportions_,Hash,"LoadProportions");
   // Calculate Rotation and Loading
   // Euler angles transformation Rotation_ = Z*Y*X
   Rotation_.Resize(DIM3,DIM3,0.0);
   Rotation_[0][0] = cos(EulerAng_[1])*cos(EulerAng_[2]);
   Rotation_[0][1] = cos(EulerAng_[2])*sin(EulerAng_[0])*sin(EulerAng_[1])
      - cos(EulerAng_[0])*sin(EulerAng_[2]);
   Rotation_[0][2] = cos(EulerAng_[0])*cos(EulerAng_[2])*sin(EulerAng_[1])
      + sin(EulerAng_[0])*sin(EulerAng_[2]);
   Rotation_[1][0] = cos(EulerAng_[1])*sin(EulerAng_[2]);
   Rotation_[1][1] = cos(EulerAng_[0])*cos(EulerAng_[2]) + sin(EulerAng_[0])
      *sin(EulerAng_[1])*sin(EulerAng_[2]);
   Rotation_[1][2] = -cos(EulerAng_[2])*sin(EulerAng_[0])
      + cos(EulerAng_[0])*sin(EulerAng_[1])*sin(EulerAng_[2]);
   Rotation_[2][0] = -sin(EulerAng_[1]);
   Rotation_[2][1] = cos(EulerAng_[1])*sin(EulerAng_[0]);
   Rotation_[2][2] = cos(EulerAng_[0])*cos(EulerAng_[1]);
   //
   // Loading_ = R*Lambda*R^T
   Loading_.Resize(DIM3,DIM3,0.0);
   for (int i=0;i<DIM3;++i)
      for (int j=0;j<DIM3;++j)
         for (int k=0;k<DIM3;++k)
            Loading_[i][j] += Rotation_[i][k]*LoadingProportions_[k]*Rotation_[j][k];
   
   
   // needed to initialize reference length
   int iter;
   iter = Input.getPosInt(Hash,"MaxIterations");
   GridSize_ = Input.getPosInt(Hash,"BlochWaveGridSize");

   // Initialize various data storage space
   ME1_static.Resize(CBK_->DOFS(),0.0);
   ME2_static.Resize(CBK_->DOFS(),CBK_->DOFS(),0.0);
   A_static.Resize(InternalAtoms_*DIM3,InternalAtoms_*DIM3);
   EigVals_static.Resize(1,InternalAtoms_*DIM3);
   InverseLat_static.Resize(DIM3,DIM3);
   Z_static.Resize(DIM3);
   str_static.Resize(CBK_->DOFS());
   stiff_static.Resize(CBK_->DOFS(),CBK_->DOFS());
   CondEV_static.Resize(1,CBK_->Fsize());
   TE_static.Resize(1,CBK_->DOFS());
   CondModuli_static.Resize(CBK_->Fsize(),CBK_->Fsize());
   TestFunctVals_static.Resize(CBK_->DOFS());
   K_static.Resize(DIM3);
   
   
   // Initiate the Lattice Sum object
   LatSum_(CBK_,InternalAtoms_,Potential_,&InfluenceDist_,&NTemp_);
   
   int err=0;
   err=FindLatticeSpacing(iter);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << "\n";
      exit(-1);
   }
   
   // Setup initial status for parameters
   NTemp_ = Input.getDouble(Hash,"NTemp");
   Lambda_ = Input.getDouble(Hash,"Lambda");
   // Make any changes to atomic potentials that might be required
   for (int i=0;i<InternalAtoms_;++i)
   {
      for (int j=i;j<InternalAtoms_;++j)
      {
         if (AtomSpecies_[i] < AtomSpecies_[j])
            UpdatePairPotential(Hash,Input,
                                AtomSpecies_[i],AtomSpecies_[j],Potential_[i][j]);
         else
            UpdatePairPotential(Hash,Input,
                                AtomSpecies_[j],AtomSpecies_[i],Potential_[j][i]);
      }
   }
   LatSum_.Recalc();
   
   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   UCIter_(GridSize_);

   Input.EndofInputSection();
}

int MultiLatticeTPP::FindLatticeSpacing(int const& iter)
{
   Lambda_=0.0;
   NTemp_=1.0;

   CBK_->SetReferenceDOFs();
   LatSum_.Recalc();

   if (Echo_)
      RefineEqbm(1.0e-13,iter,&cout);
   else
      RefineEqbm(1.0e-13,iter,0);

   // Clean up numerical round off (at least for zero values)
   Vector doftmp(CBK_->DOFS(),0.0);
   for (int i=0;i<CBK_->DOFS();++i)
   {
      if (fabs(CBK_->DOF()[i]) < 1.0e-13)
         doftmp[i] = 0.0;
      else
         doftmp[i] = CBK_->DOF()[i];
   }
   CBK_->SetDOF(doftmp);
   
   CBK_->SetReferenceToCurrent();
   
   LatSum_.Recalc();
   return 0;
}

void MultiLatticeTPP::SetParameters(double const* const Vals,int const& ResetRef)
{
   int no = SpeciesPotential_[0][0]->GetNoParameters();
   int cur = 0;
   for (int i=0;i<NumberofSpecies_;++i)
      for (int j=i;j<NumberofSpecies_;++j)
      {
         SpeciesPotential_[i][j]->SetParameters(&(Vals[cur]));
         cur += no;
      }

   LatSum_.Recalc();
   if (ResetRef)
   {
      FindLatticeSpacing(50);
   }
}

// Lattice Routines
double MultiLatticeTPP::E0() const
{
   Phi0_static = energy();
   
   if (KillTranslations_)
   {
      for (int j=0;j<DIM3;++j)
      {
         Tsq_static[j]=0.0;
         for (int i=0;i<InternalAtoms_;++i)
         {
            Tsq_static[j]+=CBK_->DOF()[CBK_->INDS(i,j)];
         }
         Tsq_static[j] = (Tsq_static[j]*Tsq_static[j])/InternalAtoms_;
      }
   }
   switch (KillRotations_)
   {
      case 2:
         // Kill three rotations
         Rsq_static[0] = (CBK_->DOF()[CBK_->INDF(0,1)] - CBK_->DOF()[CBK_->INDF(1,0)]);
         Rsq_static[1] = (CBK_->DOF()[CBK_->INDF(1,2)] - CBK_->DOF()[CBK_->INDF(2,1)]);
         Rsq_static[2] = (CBK_->DOF()[CBK_->INDF(2,0)] - CBK_->DOF()[CBK_->INDF(0,2)]);
         for (int i=0;i<DIM3;++i)
            Rsq_static[i] *= 0.5*Rsq_static[i];
         break;
      case 1:
         // Kill one rotation
         for (int i=0;i<DIM3;++i)
         {
            Rsq_static[i]=0.0;
            for (int j=0;j<DIM3;++j)
               Rsq_static[0] += KillOneRotation_[CBK_->INDF(i,j)]*CBK_->DOF()[CBK_->INDF(i,j)];
         }
         Rsq_static[0] *= Rsq_static[0];
         break;
   }
   
   return Phi0_static
      + 0.5*(TrEig_[0]*Tsq_static[0] + TrEig_[1]*Tsq_static[1] + TrEig_[2]*Tsq_static[2])
      + 0.5*(RoEig_[0]*Rsq_static[0] + RoEig_[1]*Rsq_static[1] + RoEig_[2]*Rsq_static[2]);
}

double MultiLatticeTPP::energy(PairPotentials::TDeriv const& dt) const
{
   double Phi = 0.0;
   double Vr;
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Phi
      Phi += Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
         NTemp_,LatSum_.r2(),PairPotentials::Y0,dt);
   }
   
   // Phi = Phi/(2*Vr*NormModulus)
   Vr = Density_ ? CBK_->RefVolume() : 1.0;
   Phi *= 1.0/(2.0*(Vr*NormModulus_));
   
   // Apply loading potential and Thermal term
   if (dt == PairPotentials::T0)
   {
      // Loading
      for (int i=0;i<DIM3;++i)
         for (int j=0;j<DIM3;++j)
            Phi -= Lambda_*Loading_[i][j]*((CBK_->DOF())[CBK_->INDF(j,i)] - Del(j,i));
      
      // Thermal term
      //Phi += (PhiRef_ -
      //              (NTemp_*Tref_)*EntropyRef_ -
      //              HeatCapacityRef_*(NTemp_*Tref_)*(log(NTemp_*Tref_) - 1.0)
      //         )/NormModulus_;
   }
   else if (dt == PairPotentials::DT)
   {
      // Loading
      
      // Thermal term
      //Phi += (-EntropyRef_ - HeatCapacityRef_*log(NTemp_*Tref_))/NormModulus_;
   }
   else if (dt == PairPotentials::D2T)
   {
      // Loading
      
      // Thermal term
      //Phi += (-HeatCapacityRef_/(NTemp_*Tref_))/NormModulus_;
   }
   else
   {
      cerr << "Error in MultiLatticeTPP::energy" << "\n";
      exit(-1);
   }
   
   return Phi;
}

Vector const& MultiLatticeTPP::E1() const
{
   ME1_static = stress();
   
   if (KillTranslations_)
   {
      for (int j=0;j<DIM3;++j)
      {
         T_static[j]=0.0;
         for (int i=0;i<InternalAtoms_;++i)
            T_static[j]+=CBK_->DOF()[CBK_->INDS(i,j)];
         T_static[j]/=InternalAtoms_;
      }
      for (int i=0;i<InternalAtoms_;++i)
         for (int j=0;j<DIM3;++j)
            ME1_static[CBK_->INDS(i,j)] += TrEig_[j]*T_static[j];
   }
   
   switch (KillRotations_)
   {
      case 2:
         // Kill three rotations
         R_static[0] = (CBK_->DOF()[CBK_->INDF(0,1)] - CBK_->DOF()[CBK_->INDF(1,0)])/2.0;
         ME1_static[CBK_->INDF(0,1)] += RoEig_[0]*R_static[0];
         ME1_static[CBK_->INDF(1,0)] -= RoEig_[0]*R_static[0];
         R_static[1] = (CBK_->DOF()[CBK_->INDF(1,2)] - CBK_->DOF()[CBK_->INDF(2,1)])/2.0;
         ME1_static[CBK_->INDF(1,2)] += RoEig_[1]*R_static[1];
         ME1_static[CBK_->INDF(2,1)] -= RoEig_[1]*R_static[1];
         R_static[2] = (CBK_->DOF()[CBK_->INDF(2,0)] - CBK_->DOF()[CBK_->INDF(0,2)])/2.0;
         ME1_static[CBK_->INDF(2,0)] += RoEig_[2]*R_static[2];
         ME1_static[CBK_->INDF(0,2)] -= RoEig_[2]*R_static[2];
         break;
      case 1:
         // Kill one rotation
         for (int i=0;i<DIM3;++i)
         {
            R_static[i]=0.0;
            for (int j=0;j<DIM3;++j)
               R_static[0] += KillOneRotation_[CBK_->INDF(i,j)]*CBK_->DOF()[CBK_->INDF(i,j)];
         }
         
         for (int i=0;i<ME1_static.Dim();++i)
            ME1_static[i] += RoEig_[0]*R_static[0]*KillOneRotation_[i];
         break;
   }
   
   return ME1_static;
}

Vector const& MultiLatticeTPP::stress(PairPotentials::TDeriv const& dt,LDeriv const& dl) const
{
   double ForceNorm = 0.0;
   double phi,Vr;
   int i,j;
   
   S_static.Resize(CBK_->DOFS(),0.0);
   
   Vr = Density_ ? CBK_->RefVolume() : 1.0;
   
   if (dl==L0)
   {
      for (i=0;i<InternalAtoms_;++i)
      {
         for (j=0;j<DIM3;++j)
         {
            BodyForce_[i][j] = 0.0;
         }
      }
      
      for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
      {
         // Calculate bodyforce
         // NOTE: phi1 = d(phi)/d(r2)
         // We need d(phi)/dr = 2*r*d(phi)/d(r2)
         phi = 2.0*sqrt(LatSum_.r2())*LatSum_.phi1();
         if (ForceNorm < fabs(-phi/2.0))
         {
            ForceNorm = fabs(-phi/2.0);
         }
         for (i=0;i<DIM3;i++)
         {
            BodyForce_[LatSum_.Atom(0)][i] += -phi*LatSum_.Dx(i)/(2.0*sqrt(LatSum_.r2()));
         }
         
         // Claculate the Stress
         if (dt == PairPotentials::T0)
            phi=LatSum_.phi1();
         else if (dt == PairPotentials::DT)
            phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
               NTemp_,LatSum_.r2(),PairPotentials::DY,dt);
         else
         {
            cerr << "Error in MultiLatticeTPP::stress" << "\n";
            exit(-1);
         }
         
         for (i=0;i<DIM3;i++)
         {
            for (j=0;j<DIM3;j++)
            {
               S_static[CBK_->INDF(i,j)] += phi*CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j);
            }
         }
         for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         {
            for (j=0;j<DIM3;j++)
            {
               S_static[CBK_->INDS(i,j)] += phi*CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                       LatSum_.Atom(1),i,j);
            }
         }
      }
      
      // BodyForce[i] = BodyForce[i] / ForceNorm
      for (i=0;i<InternalAtoms_;i++)
      {
         for (j=0;j<DIM3;j++)
         {
            BodyForce_[i][j] /= ForceNorm;
         }
      }
      
      // S = S/(2*Vr*NormModulus)
      S_static *= 1.0/(2.0*(Vr*NormModulus_));
      
      // Load terms
      if (dt == PairPotentials::T0)
      {
         for (i=0;i<DIM3;++i)
            for (j=0;j<DIM3;++j)
            {
               S_static[CBK_->INDF(i,j)] -= Lambda_*Loading_[j][i];
            }
      }
      
   }
   else if (dl==DL)
   {
      // dl=DL
      for (i=0;i<DIM3;++i)
         for (j=0;j<DIM3;++j)
            S_static[CBK_->INDF(i,j)] -= Loading_[j][i];
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stress()" << "\n";
      exit(-1);
   }
   
   return S_static;
}

Matrix const& MultiLatticeTPP::E2() const
{
   ME2_static = stiffness();
   
   if (KillTranslations_)
   {
      
      for (int i=0;i<InternalAtoms_;++i)
         for (int j=0;j<DIM3;++j)
            for (int k=0;k<InternalAtoms_;++k)
            {
               ME2_static[CBK_->INDS(i,j)][CBK_->INDS(k,j)] += TrEig_[j]/InternalAtoms_;
            }
   }
   
   switch (KillRotations_)
   {
      case 2:
         // Kill three rotations
         ME2_static[CBK_->INDF(0,1)][CBK_->INDF(0,1)] += RoEig_[0]/2.0;
         ME2_static[CBK_->INDF(0,1)][CBK_->INDF(1,0)] -= RoEig_[0]/2.0;
         ME2_static[CBK_->INDF(1,0)][CBK_->INDF(1,0)] += RoEig_[0]/2.0;
         ME2_static[CBK_->INDF(1,0)][CBK_->INDF(0,1)] -= RoEig_[0]/2.0;
         
         ME2_static[CBK_->INDF(1,2)][CBK_->INDF(1,2)] += RoEig_[1]/2.0;
         ME2_static[CBK_->INDF(1,2)][CBK_->INDF(2,1)] -= RoEig_[1]/2.0;
         ME2_static[CBK_->INDF(2,1)][CBK_->INDF(2,1)] += RoEig_[1]/2.0;
         ME2_static[CBK_->INDF(2,1)][CBK_->INDF(1,2)] -= RoEig_[1]/2.0;
         
         ME2_static[CBK_->INDF(2,0)][CBK_->INDF(2,0)] += RoEig_[2]/2.0;
         ME2_static[CBK_->INDF(2,0)][CBK_->INDF(0,2)] -= RoEig_[2]/2.0;
         ME2_static[CBK_->INDF(0,2)][CBK_->INDF(0,2)] += RoEig_[2]/2.0;
         ME2_static[CBK_->INDF(0,2)][CBK_->INDF(2,0)] -= RoEig_[2]/2.0;
         break;
      case 1:
         // Kill one rotation
         for (int i=0;i<ME2_static.Rows();++i)
            for (int j=0;j<ME2_static.Cols();++j)
               ME2_static[i][j] += RoEig_[0]*KillOneRotation_[i]*KillOneRotation_[j];
         break;
   }
   
   return ME2_static;
}

Matrix const& MultiLatticeTPP::stiffness(PairPotentials::TDeriv const& dt,LDeriv const& dl)
   const
{
   Matrix F(DIM3,DIM3);
   double phi,phi1;
   int i,j,k,l;
   
   Phi2_static.Resize(CBK_->DOFS(),CBK_->DOFS(),0.0);
   
   if (dl==L0)
   {
      for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
      {
         if (dt==PairPotentials::T0)
         {
            phi = LatSum_.phi2();
            phi1 = LatSum_.phi1();
         }
         else if (dt==PairPotentials::DT)
         {
            phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
               NTemp_,LatSum_.r2(),PairPotentials::D2Y,dt);
            phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
               NTemp_,LatSum_.r2(),PairPotentials::DY,dt);
         }
         else
         {
            cerr << "Error in MultiLatticeTPP::stiffness" << "\n";
            exit(-1);
         }
         
         //Upper Diag Block (CBK_->Fsize(),CBK_->Fsize())
         for (i=0;i<DIM3;i++)
         {
            for (j=0;j<DIM3;j++)
            {
               for (k=0;k<DIM3;k++)
               {
                  for (l=0;l<DIM3;l++)
                  {
                     Phi2_static[CBK_->INDF(i,j)][CBK_->INDF(k,l)]+=
                        phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                             *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l))
                        +phi1*CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l);
                  }
               }
            }
         }
         //Lower Diag Block (CBK_->Ssize(),CBK_->Ssize())
         for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         {
            for (j=0;j<DIM3;j++)
            {
               for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               {
                  for (l=0;l<DIM3;l++)
                  {
                     Phi2_static[CBK_->INDS(i,j)][CBK_->INDS(k,l)]+=
                        phi*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
                             *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l))
                        +phi1*CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
                  }
               }
            }
         }
         //Off Diag Blocks
         for (i=0;i<DIM3;i++)
         {
            for (j=0;j<DIM3;j++)
            {
               for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               {
                  for (l=0;l<DIM3;l++)
                  {
                     Phi2_static[CBK_->INDF(i,j)][CBK_->INDS(k,l)] =
                        Phi2_static[CBK_->INDS(k,l)][CBK_->INDF(i,j)] +=
                        phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                             *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l))
                        +phi1*CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                           LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
                  }
               }
            }
         }
      }
      
      // Phi2_static = Phi2_static/(2*Vr*NormModulus)
      Phi2_static *= 1.0/((2.0*(Density_ ? CBK_->RefVolume() : 1.0)*NormModulus_));
   }
   else if (dl==DL)
   {
      // Nothing to do: Phi2_static is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stiffness()" << "\n";
      exit(-1);
   }
   return Phi2_static;
}

Matrix const& MultiLatticeTPP::E3() const
{
   double phi,phi1,phi2;
   int i,j,k,l,m,n;
   
   Phi3_static.Resize(CBK_->DOFS()*CBK_->DOFS(),CBK_->DOFS(),0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
         NTemp_,LatSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi1=LatSum_.phi2();
      phi2=LatSum_.phi1();
      
      // DF^3 block
      for (i=0;i<DIM3;i++)
         for (j=0;j<DIM3;j++)
            for (k=0;k<DIM3;k++)
               for (l=0;l<DIM3;l++)
                  for (m=0;m<DIM3;m++)
                     for (n=0;n<DIM3;n++)
                     {
                        Phi3_static[CBK_->INDFF(i,j,k,l)][CBK_->INDF(m,n)] +=
                           phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n))
                           +phi1*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                  *CBK_->D2yDFF(LatSum_.pDX(),k,l,m,n) +
                                  CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                  *CBK_->D2yDFF(LatSum_.pDX(),i,j,m,n) +
                                  CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                  *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l));
                     }
      // DS^3 block
      for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         for (j=0;j<DIM3;j++)
            for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               for (l=0;l<DIM3;l++)
                  for (m=CBK_->NoTrans();m<InternalAtoms_;m++)
                     for (n=0;n<DIM3;n++)
                     {
                        Phi3_static[CBK_->INDSS(i,j,k,l)][CBK_->INDS(m,n)] +=
                           phi*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
                                *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l)
                                *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),m,n))
                           +phi1*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                             LatSum_.Atom(1),i,j)
                                  *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                  +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                              LatSum_.Atom(1),k,l)
                                  *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                  +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                              LatSum_.Atom(1),m,n)
                                  *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l));
                     }
      // DF^2DS blocks
      for (i=0;i<DIM3;i++)
         for (j=0;j<DIM3;j++)
            for (k=0;k<DIM3;k++)
               for (l=0;l<DIM3;l++)
                  for (m=CBK_->NoTrans();m<InternalAtoms_;m++)
                     for (n=0;n<DIM3;n++)
                     {
                        Phi3_static[CBK_->INDFF(i,j,k,l)][CBK_->INDS(m,n)] =
                           Phi3_static[CBK_->INDFS(i,j,m,n)][CBK_->INDF(k,l)] =
                           Phi3_static[CBK_->INDSF(m,n,i,j)][CBK_->INDF(k,l)] += (
                              phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                   *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                   *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                               LatSum_.Atom(1),m,n))
                              +phi1*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                     *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                   LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                     +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                     *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                   LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                     +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                 LatSum_.Atom(1),m,n)
                                     *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l))
                              +phi2*CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),i,j,k,l,m,n));
                     }
      // DS^2DF blocks
      for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         for (j=0;j<DIM3;j++)
            for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               for (l=0;l<DIM3;l++)
                  for (m=0;m<DIM3;m++)
                     for (n=0;n<DIM3;n++)
                     {
                        Phi3_static[CBK_->INDSS(i,j,k,l)][CBK_->INDF(m,n)] =
                           Phi3_static[CBK_->INDSF(i,j,m,n)][CBK_->INDS(k,l)] =
                           Phi3_static[CBK_->INDFS(m,n,i,j)][CBK_->INDS(k,l)] += (
                              phi*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                              LatSum_.Atom(1),i,j)
                                   *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                               LatSum_.Atom(1),k,l)
                                   *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n))
                              +phi1*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),i,j)
                                     *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                   LatSum_.Atom(0),LatSum_.Atom(1),m,n,k,l)
                                     +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                 LatSum_.Atom(1),k,l)
                                     *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                   LatSum_.Atom(0),LatSum_.Atom(1),m,n,i,j)
                                     + CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                     *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
                              +phi2*CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l,m,n));
                     }
   }
   
   
   // Phi3_static = Phi3_static/(2*Vr*NormModulus)
   Phi3_static *= 1.0/(2.0*((Density_ ? CBK_->RefVolume() : 1.0)*NormModulus_));
   
   return Phi3_static;
}

Matrix const& MultiLatticeTPP::E4() const
{
   double phi,phi1,phi2,phi3;
   int i,j,k,l,m,n,s,t;
   
   Phi4_static.Resize(CBK_->DOFS()*CBK_->DOFS(),CBK_->DOFS()*CBK_->DOFS(),0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
         NTemp_,LatSum_.r2(),PairPotentials::D4Y,PairPotentials::T0);
      phi1=Potential_[LatSum_.Atom(0)][LatSum_.Atom(1)]->PairPotential(
         NTemp_,LatSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi2=LatSum_.phi2();
      phi3=LatSum_.phi1();
      
      // DF^4 block
      for (i=0;i<DIM3;i++)
         for (j=0;j<DIM3;j++)
            for (k=0;k<DIM3;k++)
               for (l=0;l<DIM3;l++)
                  for (m=0;m<DIM3;m++)
                     for (n=0;n<DIM3;n++)
                        for (s=0;s<DIM3;s++)
                           for (t=0;t<DIM3;t++)
                           {
                              Phi4_static[CBK_->INDFF(i,j,k,l)][CBK_->INDFF(m,n,s,t)]+=
                                 phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                      *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                      *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                      *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)) +
                                 phi1*(
                                    CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                    *CBK_->D2yDFF(LatSum_.pDX(),m,n,s,t)
                                    +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                    *CBK_->D2yDFF(LatSum_.pDX(),k,l,s,t)
                                    +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                    *CBK_->D2yDFF(LatSum_.pDX(),i,j,s,t)
                                    +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),k,l,m,n)
                                    +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),i,j,m,n)
                                    +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                    *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l))
                                 +phi2*(
                                    CBK_->D2yDFF(LatSum_.pDX(),i,j,s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),k,l,m,n)
                                    +CBK_->D2yDFF(LatSum_.pDX(),k,l,s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),i,j,m,n)
                                    +CBK_->D2yDFF(LatSum_.pDX(),m,n,s,t)
                                    *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l));
                           }
      // DS^4 block
      for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         for (j=0;j<DIM3;j++)
            for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               for (l=0;l<DIM3;l++)
                  for (m=CBK_->NoTrans();m<InternalAtoms_;m++)
                     for (n=0;n<DIM3;n++)
                        for (s=CBK_->NoTrans();s<InternalAtoms_;s++)
                           for (t=0;t<DIM3;t++)
                           {
                              Phi4_static[CBK_->INDSS(i,j,k,l)][CBK_->INDSS(m,n,s,t)] +=
                                 phi*(CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                 LatSum_.Atom(1),i,j)
                                      *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),k,l)
                                      *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),m,n)
                                      *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),s,t))
                                 +phi1*(
                                    CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                               LatSum_.Atom(1),i,j)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),k,l)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
                                    +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),i,j)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),m,n)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                    +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),k,l)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),m,n)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
                                    +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),i,j)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),s,t)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                    +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),k,l)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),s,t)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                    +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),m,n)
                                    *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                LatSum_.Atom(1),s,t)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
                                 +phi2*(
                                    CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                    +CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                    +CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l)
                                    *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t));
                           }
      // DF^3DS blocks
      for (i=0;i<DIM3;i++)
         for (j=0;j<DIM3;j++)
            for (k=0;k<DIM3;k++)
               for (l=0;l<DIM3;l++)
                  for (m=0;m<DIM3;m++)
                     for (n=0;n<DIM3;n++)
                        for (s=CBK_->NoTrans();s<InternalAtoms_;s++)
                           for (t=0;t<DIM3;t++)
                           {
                              Phi4_static[CBK_->INDFF(i,j,k,l)][CBK_->INDFS(m,n,s,t)] =
                                 Phi4_static[CBK_->INDFF(i,j,k,l)][CBK_->INDSF(s,t,m,n)] =
                                 Phi4_static[CBK_->INDFS(i,j,s,t)][CBK_->INDFF(k,l,m,n)] =
                                 Phi4_static[CBK_->INDSF(s,t,i,j)][CBK_->INDFF(k,l,m,n)] += (
                                    phi*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t))
                                    +phi1*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
                                       +(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                         *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                     LatSum_.Atom(1),s,t)
                                         *CBK_->D2yDFF(LatSum_.pDX(),k,l,m,n)
                                         +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                         *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                     LatSum_.Atom(1),s,t)
                                         *CBK_->D2yDFF(LatSum_.pDX(),i,j,m,n)
                                         +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                         *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                     LatSum_.Atom(1),s,t)
                                         *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l)))
                                    +phi2*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                      LatSum_.Atom(1),k,l,m,n,s,t)
                                       +CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),
                                                     i,j,s,t)
                                       *CBK_->D2yDFF(LatSum_.pDX(),k,l,m,n)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                      LatSum_.Atom(1),i,j,m,n,s,t)
                                       +CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                       *CBK_->D2yDFF(LatSum_.pDX(),i,j,m,n)
                                       + CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),m,n)
                                       *CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                      LatSum_.Atom(1),i,j,k,l,s,t)
                                       + CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                      LatSum_.Atom(0),LatSum_.Atom(1),
                                                      m,n,s,t)
                                       *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l)));
                           }
      // DS^3DF blocks
      for (i=CBK_->NoTrans();i<InternalAtoms_;i++)
         for (j=0;j<DIM3;j++)
            for (k=CBK_->NoTrans();k<InternalAtoms_;k++)
               for (l=0;l<DIM3;l++)
                  for (m=CBK_->NoTrans();m<InternalAtoms_;m++)
                     for (n=0;n<DIM3;n++)
                        for (s=0;s<DIM3;s++)
                           for (t=0;t<DIM3;t++)
                           {
                              Phi4_static[CBK_->INDSS(i,j,k,l)][CBK_->INDSF(m,n,s,t)] =
                                 Phi4_static[CBK_->INDSS(i,j,k,l)][CBK_->INDFS(s,t,m,n)] =
                                 Phi4_static[CBK_->INDSF(i,j,s,t)][CBK_->INDSS(k,l,m,n)] =
                                 Phi4_static[CBK_->INDFS(s,t,i,j)][CBK_->INDSS(k,l,m,n)] += (
                                    phi*(
                                       CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),i,j)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),k,l)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t))
                                    +phi1*(
                                       CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),i,j)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),k,l)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,m,n)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),i,j)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,k,l)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),k,l)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,i,j)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),k,l)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),s,t)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l))
                                    +phi2*(
                                       CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                  LatSum_.Atom(1),i,j)
                                       *CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),
                                                      k,l,m,n,s,t)
                                       +CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,i,j)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                       + CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,k,l)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),k,l)
                                       *CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),
                                                      i,j,m,n,s,t)
                                       +CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),s,t,m,n)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),
                                                      i,j,k,l,s,t)));
                           }
      // DF^2DS^2 blocks
      for (i=0;i<DIM3;i++)
         for (j=0;j<DIM3;j++)
            for (k=0;k<DIM3;k++)
               for (l=0;l<DIM3;l++)
                  for (m=CBK_->NoTrans();m<InternalAtoms_;m++)
                     for (n=0;n<DIM3;n++)
                        for (s=CBK_->NoTrans();s<InternalAtoms_;s++)
                           for (t=0;t<DIM3;t++)
                           {
                              Phi4_static[CBK_->INDFF(i,j,k,l)][CBK_->INDSS(m,n,s,t)] =
                                 Phi4_static[CBK_->INDFS(i,j,m,n)][CBK_->INDFS(k,l,s,t)] =
                                 Phi4_static[CBK_->INDFS(i,j,m,n)][CBK_->INDSF(s,t,k,l)] =
                                 Phi4_static[CBK_->INDSF(m,n,i,j)][CBK_->INDFS(k,l,s,t)] =
                                 Phi4_static[CBK_->INDSF(m,n,i,j)][CBK_->INDSF(s,t,k,l)] =
                                 Phi4_static[CBK_->INDSS(m,n,s,t)][CBK_->INDFF(i,j,k,l)] += (
                                    phi*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t))
                                    +phi1*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                       + CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                    LatSum_.Atom(1),m,n)
                                       *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t)
                                       *CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l))
                                    +phi2*(
                                       CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
                                       *CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),
                                                      m,n,s,t,k,l)
                                       +CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),i,j,s,t)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,m,n)
                                       +CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),i,j,m,n)
                                       *CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
                                                     LatSum_.Atom(0),LatSum_.Atom(1),k,l,s,t)
                                       +CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l)
                                       *CBK_->D3yDSSF(LatSum_.Atom(0),LatSum_.Atom(1),
                                                      m,n,s,t,i,j)
                                       +CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l)
                                       *CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),m,n,s,t)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),m,n)
                                       *CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                      LatSum_.Atom(1),i,j,k,l,s,t)
                                       +CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
                                                   LatSum_.Atom(1),s,t)
                                       *CBK_->D3yDFFS(LatSum_.pDX(),LatSum_.Atom(0),
                                                      LatSum_.Atom(1),i,j,k,l,m,n))
                                    +phi3*CBK_->D4yDFFSS(LatSum_.Atom(0),LatSum_.Atom(1),
                                                         i,j,k,l,m,n,s,t));
                           }
   }
   
   
   // Phi4_static = Phi4_static/(2*Vr*NormModulus)
   Phi4_static *= 1.0/(2.0*((Density_ ? CBK_->RefVolume() : 1.0)*NormModulus_));
   
   return Phi4_static;
}

Matrix const& MultiLatticeTPP::CondensedModuli() const
{
   Matrix const& stiff = stiffness();
   int intrn = CBK_->Ssize();
   double factor = 1.0/(intrn/DIM3);
   int fsz = CBK_->Fsize();
   Matrix IM(intrn,intrn);
   CM_static.Resize(fsz,fsz);
   
   for (int i=0;i<fsz;i++)
      for (int j=0;j<fsz;j++)
      {
         CM_static[i][j] = stiff[i][j];
      }
   
   // Make sure there are internal DOF's
   if (intrn)
   {
      for (int i=0;i<intrn;i++)
         for (int j=0;j<intrn;j++)
         {
            IM[i][j] = stiff[fsz+i][fsz+j];
            
            // add translational stiffness to regularize IM, if needed
            if ((!CBK_->NoTrans()) && (i%DIM3 == j%DIM3))
            {
               IM[i][j] += factor;
            }
         }
      IM = IM.Inverse();
      
      // Set up Condensed Moduli
      for (int i=0;i<fsz;i++)
         for (int j=0;j<fsz;j++)
         {
            for (int m=0;m<intrn;m++)
               for (int n=0;n<intrn;n++)
               {
                  CM_static[i][j] -= stiff[i][fsz+m]*IM[m][n]*stiff[fsz+n][j];
               }
         }
   }
   
   // If using symmetrized F, assume standard Voigt notation
   if (fsz==6)
   {
      // Remove 2's and 4's
      for (int i=3;i<6;i++)
      {
         for (int j=0;j<3;j++)
         {
            CM_static[i][j] /= 2.0;
            CM_static[j][i] /= 2.0;
         }
         
         for (int j=3;j<6;j++)
         {
            CM_static[i][j] /= 4.0;
         }
      }
   }
   
   return CM_static;
}

Vector const& MultiLatticeTPP::ThermalExpansion() const
{
   ThermalExp_static.Resize(1,CBK_->DOFS());
#ifdef SOLVE_SVD
   return ThermalExp_static = SolveSVD(E2(),-StressDT());
#else
   return ThermalExp_static = SolvePLU(E2(),-StressDT());
#endif
}

int MultiLatticeTPP::comp(void const* const a,void const* const b)
{
   double t;
   if( *((double*) a) == *((double*) b)) return 0;
   else
   {
      t= *((double*) a) - *((double*) b);
      t/=fabs(t);
      return int(t);
   }
}

int MultiLatticeTPP::abscomp(void const* const a,void const* const b)
{
   double t;
   if( fabs(*((double*) a)) == fabs(*((double*) b))) return 0;
   else
   {
      t= fabs(*((double*) a)) - fabs(*((double*) b));
      t/=fabs(t);
      return int(t);
   }
}

void MultiLatticeTPP::interpolate(Matrix* const EigVals,int const& zero,int const& one,
                                  int const& two)
{
   // Calculate expected value for eigvals and store in zero position
   EigVals[zero] = 2.0*EigVals[one] - EigVals[zero];
   
   double delta,dtmp;
   int i,j,pos;
   
   for (i=0;i<EigVals[0].Cols();++i)
   {
      pos = i;
      delta = fabs(EigVals[zero][0][i] - EigVals[two][0][i]);
      for (j=i+1;j<EigVals[0].Cols();++j)
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

CMatrix const& MultiLatticeTPP::ReferenceDynamicalStiffness(Vector const& K) const
{
   double pi = 4.0*atan(1.0);
   MyComplexDouble Ic(0,1);
   MyComplexDouble A = 2.0*pi*Ic;
   int i,j;
   
   Dk_static.Resize(InternalAtoms_*DIM3,InternalAtoms_*DIM3,0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Dk_static
      if (LatSum_.Atom(0) != LatSum_.Atom(1))
      {
         for (i=0;i<DIM3;++i)
            for (j=0;j<DIM3;++j)
            {
               // y != y' terms (i.e., off block (3x3) diagonal terms)
               Dk_static[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
                  (-2.0*Del(i,j)*LatSum_.phi1()
                   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2())
                  *exp(A*
                       (K[0]*LatSum_.DX(0) + K[1]*LatSum_.DX(1) + K[2]*LatSum_.DX(2)));
               
               // y==y' components (i.e., Phi(0,y,y) term)
               Dk_static[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(0)+j] +=
                  (2.0*Del(i,j)*LatSum_.phi1()
                   +4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2());
            }
      }
      else
      {
         for (i=0;i<DIM3;++i)
            for (j=0;j<DIM3;++j)
            {
               Dk_static[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
                  (-2.0*Del(i,j)*LatSum_.phi1()
                   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2())
                  *(exp(A*
                        (K[0]*LatSum_.DX(0) + K[1]*LatSum_.DX(1)
                         + K[2]*LatSum_.DX(2)))
                    - 1.0);
            }
      }
   }
   // Normalize through the Mass Matrix
   for (int p=0;p<InternalAtoms_;++p)
      for (int q=0;q<InternalAtoms_;++q)
         for (i=0;i<DIM3;++i)
            for (j=0;j<DIM3;++j)
            {
               Dk_static[DIM3*p+i][DIM3*q+j] /= sqrt(AtomicMass_[p]*AtomicMass_[q]);
            }
   
   return Dk_static;
}

void MultiLatticeTPP::ReferenceDispersionCurves(Vector const& K,int const& NoPTS,
                                                char const* const prefix,ostream& out) const
{
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);
   
   Matrix InverseLat(DIM3,DIM3);
   InverseLat = (CBK_->RefLattice()).Inverse();
   
   Matrix EigVal[DIM3];
   for (int i=0;i<DIM3;++i) EigVal[i].Resize(1,InternalAtoms_*DIM3);
   
   Vector Z1(DIM3),Z2(DIM3);
   for (int k=0;k<DIM3;++k)
   {
      Z1[k] = K[k];
      Z2[k] = K[DIM3 + k];
   }
   Z1 = InverseLat*Z1;
   Z2 = InverseLat*Z2;
   
   Vector Z(DIM3),
      DZ=Z2-Z1;
   double dz = 1.0/(NoPTS-1);
   for (int k=0;k<2;++k)
   {
      Z = Z1 + (k*dz)*DZ;
      EigVal[k] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[k][0],InternalAtoms_*DIM3,sizeof(double),&comp);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<InternalAtoms_*DIM3;++i)
      {
         out << setw(w) << EigVal[k][0][i];
         if (Echo_) cout << setw(w) << EigVal[k][0][i];
      }
      out << "\n";
      if (Echo_) cout << "\n";
   }
   int zero=0,one=1,two=2;
   for (int k=2;k<NoPTS;++k)
   {
      Z = Z1 + (k*dz)*DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[two][0],InternalAtoms_*DIM3,sizeof(double),&comp);
      interpolate(EigVal,zero,one,two);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<InternalAtoms_*DIM3;++i)
      {
         out << setw(w) << EigVal[two][0][i];;
         if (Echo_) cout << setw(w) << EigVal[two][0][i];;
      }
      out << "\n";
      if (Echo_) cout << "\n";
      
      zero = (zero+1)%3; one = (zero+1)%3; two = (one+1)%3;
   }
}

int MultiLatticeTPP::ReferenceBlochWave(Vector& K) const
{
   InverseLat_static = (CBK_->RefLattice()).Inverse();
   
   // Iterate over points in cubic unit cell
   for (UCIter_.Reset();!UCIter_.Done();++UCIter_)
   {
      for (int i=0;i<DIM3;++i)
      {
         K[i] = UCIter_[i];
      }
      
      Z_static = InverseLat_static*K;
      A_static = ReferenceDynamicalStiffness(Z_static);
      
      EigVals_static = HermiteEigVal(A_static);
      
      for (int i=0;i<InternalAtoms_*DIM3;++i)
      {
         // if w^2 <= 0.0 --> Re(i*w*x) > 0 --> growing solutions --> unstable
         if ( EigVals_static[0][i] <= 0.0 )
         {
            return 0;
         }
      }
   }
   return 1;
}

//---- needs to be updated----//
void MultiLatticeTPP::LongWavelengthModuli(double const& dk,int const& gridsize,
                                           char const* const prefix,ostream& out) const
{
   double pi = 4*atan(1.0);
   double twopi = 2*pi;
   double GS = double(gridsize);
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);
   
   Matrix
      Lp=CondensedModuli(),
      Ap(DIM3,DIM3),
      A(DIM3,DIM3);
   
   //----------------  setup L condensed moduli wrt F NOT U ----------------------
   Matrix Phi(9,9,0.0),
      Dpp((InternalAtoms_-1)*3,(InternalAtoms_-1)*3,0.0),
      Dfp((InternalAtoms_-1)*3,9,0.0);
   double phi,phi1,tmp[3][3][3];
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      phi = LatSum_.phi2();
      phi1 = LatSum_.phi1();
      
      // upper 9x9 block
      for (int i=0;i<DIM3;i++)
      {
         for (int j=0;j<DIM3;j++)
         {
            for (int k=0;k<DIM3;k++)
            {
               for (int l=0;l<DIM3;l++)
               {
                  Phi[3*i+j][3*k+l]+=
                     4.0*phi*(LatSum_.Dx(i)*LatSum_.DX(j))
                     *(LatSum_.Dx(k)*LatSum_.DX(l))
                     +2*phi1*(Del(i,k)*LatSum_.DX(j)*LatSum_.DX(l));
               }
            }
         }
      }
      
      // lower block
      for (int i=1;i<InternalAtoms_;++i)
      {
         for (int j=0;j<DIM3;++j)
         {
            for (int k=1;k<InternalAtoms_;++k)
            {
               for (int l=0;l<DIM3;++l)
               {
                  Dpp[3*(i-1)+j][3*(k-1)+l]+=
                     phi*CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),i,j)
                     *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l)
                     +phi1*CBK_->D2yDSS(LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
               }
            }
         }
      }
      
      // off-diagonal block
      for (int i=0;i<DIM3;++i)
         for (int j=0;j<DIM3;++j)
            for (int l=0;l<DIM3;++l)
            {
               tmp[i][j][l] = 0.0;
               for (int k=0;k<DIM3;++k)
               {
                  tmp[i][j][l] += (CBK_->RefLattice())[l][k]*(CBK_->DOF())[CBK_->INDF(k,i)]*LatSum_.DX(j)
                     +(CBK_->RefLattice())[l][k]*(CBK_->DOF())[CBK_->INDF(i,k)]*LatSum_.DX(j);
               }
            }
      for (int k=1;k<InternalAtoms_;++k)
         for (int l=0;l<DIM3;++l)
            for (int i=0;i<DIM3;++i)
               for (int j=0;j<DIM3;++j)
               {
                  Dfp[3*(k-1)+l][3*i+j] +=
                     phi*(2.0*LatSum_.Dx(i)*LatSum_.DX(j))
                     *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l) +
                     phi1*0.5*(Del(k,LatSum_.Atom(1)) - Del(k,LatSum_.Atom(0)))
                     *( (CBK_->RefLattice())[l][i]*LatSum_.Dx(j) +
                        (CBK_->RefLattice())[l][j]*LatSum_.Dx(i) +
                        tmp[i][j][l]);
               }
   }
   
   // Condense moduli.
   Phi -= (Dfp.Transpose())*(Dpp.Inverse())*Dfp;
   // Phi = Phi/(2*Vr*NormModulus)
   Phi *= 1.0/(2.0*((Density_ ? CBK_->RefVolume() : 1.0)*NormModulus_));
   
   //-----------------------------------------------------------------------------
   
   Matrix SymPhi(6,6,0.0);
   
   for (int i=0;i<DIM3;++i)
      for (int j=i;j<DIM3;++j)
         for (int k=0;k<DIM3;++k)
            for (int l=k;l<DIM3;++l)
            {
               SymPhi[CBK_->INDF(i,j)][CBK_->INDF(k,l)] = 0.25*(
                  Phi[3*i+j][3*k+l] +
                  Phi[3*j+i][3*k+l] +
                  Phi[3*i+j][3*l+k] +
                  Phi[3*j+i][3*l+k]);
            }
   
   cout << "Condensed Moduli Check!:" << "\n";
   cout << setw(w) << Lp-SymPhi << "\n";
   
   
   
   Vector K(DIM3),Z(DIM3,0.0);
   Matrix BlkEigVal(1,InternalAtoms_*DIM3);
   Matrix ModEigVal(1,DIM3);
   
   double Mc=0.0;
   double Vc= Density_ ? CBK_->RefVolume() : 1.0;
   for (int i=0;i<InternalAtoms_;++i)
   {
      Mc += AtomicMass_[i];
   }
   
   for (int phi=0;phi<gridsize;++phi)
   {
      for (int theta=0;theta<gridsize;++theta)
      {
         K[0] = sin(pi*(phi/GS))*cos(twopi*(theta/GS));
         K[1] = sin(pi*(phi/GS))*sin(twopi*(theta/GS));
         K[2] = cos(pi*(phi/GS));
         
         Z=dk*K;
         BlkEigVal = HermiteEigVal(ReferenceDynamicalStiffness(Z));
         
         // sort by absolute value
         qsort(BlkEigVal[0],InternalAtoms_*DIM3,sizeof(double),&abscomp);
         
         for (int i=0;i<DIM3;++i)
         {
            // wave speed squared
            BlkEigVal[0][i] /= (twopi*dk*twopi*dk);
         }
         
         for (int i=0;i<DIM3;++i)
            for (int j=0;j<DIM3;++j)
            {
               A[i][j] = 0.0;
               for (int k=0;k<DIM3;++k)
                  for (int l=0;l<DIM3;++l)
                  {
                     A[i][j] += Phi[3*i+k][3*j+l]*K[k]*K[l];
                  }
            }
         
         ModEigVal = SymEigVal(A);
         qsort(ModEigVal[0],DIM3,sizeof(double),&abscomp);
         for (int i=0;i<3;++i)
         {
            // normalize by G/(Mc/Vc)
            ModEigVal[0][i] *= NormModulus_/(Mc/Vc);
         }
         
         out << prefix << setw(w/2) << phi << setw(w/2) << theta;
         if (Echo_) cout << prefix << setw(w/2) << phi << setw(w/2) << theta;
         for (int i=0;i<DIM3;++i)
         {
            out << setw(w) << ModEigVal[0][i];
            if (Echo_) cout << setw(w) << ModEigVal[0][i];
         }
         for (int i=0;i<DIM3;++i)
         {
            out << setw(w) << BlkEigVal[0][i];
            if (Echo_) cout << setw(w) << BlkEigVal[0][i];
         }
         for (int i=0;i<DIM3;++i)
         {
            out << setw(w) << (ModEigVal[0][i]-BlkEigVal[0][i])/ModEigVal[0][i];
            if (Echo_) cout << setw(w)
                            << (ModEigVal[0][i]-BlkEigVal[0][i])/ModEigVal[0][i];
         }
         out << "\n";
         if (Echo_) cout << "\n";
      }
      out << "\n";
      if (Echo_) cout << "\n";
   }
}

void MultiLatticeTPP::NeighborDistances(int const& cutoff,ostream& out) const
{
   Matrix NeighborDist =
      LatSum_.NeighborDistances(cutoff,pow(double(10),double(-(out.precision()-1))));
   
   int W=out.width();
   int types = (InternalAtoms_*(InternalAtoms_+1))/2;
   for (int i=0;i<cutoff;++i)
   {
      out << setw(W) << NTemp_ << setw(W) << NeighborDist[i][0];
      for (int j=0;j<types;++j)
      {
         out << setw(W/4) << int(NeighborDist[i][1+j]);
      }
      out << "\n";
   }
   out << "\n";
}

void MultiLatticeTPP::Print(ostream& out,PrintDetail const& flag)
{
   int W;
   int NoNegTestFunctions;
   double engy,entropy,heatcapacity;
   int RankOneConvex;
   int BlochWaveStable;
   

   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   engy = energy();
   TE_static = ThermalExpansion();
   entropy = Entropy();
   heatcapacity = HeatCapacity();
   str_static = stress();
   stiff_static = stiffness();

   NoNegTestFunctions=TestFunctions(TestFunctVals_static,LHS);
   CondModuli_static = CondensedModuli();
   
   CondEV_static=SymEigVal(CondModuli_static);
   RankOneConvex = FullScanRank1Convex3D(CBK_,CondModuli_static,ConvexityDX_);
   
   K_static.Resize(DIM3,0.0);
   if (RankOneConvex)
   {
      BlochWaveStable = BlochWave(K_static);
   }
   else
   {
      BlochWaveStable = -1;
   }
   
   switch (flag)
   {
      case PrintLong:
         out << "MultiLatticeTPP:" << "\n" << "\n";
         out << "Density_ = " << Density_ << "\n";
         out << "Using: " << (*CBK_) << " Kinematics" << "\n";
         out << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
         for (int i=0;i<InternalAtoms_;++i)
         {
            out << "Atom_" << i << "          "
                << "Species : " << setw(5) << AtomSpecies_[i]
                << "          Position : " << setw(W) << CBK_->AtomPositions(i) << "\n";
         }
         out << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
         for (int i=0;i<NumberofSpecies_;++i)
         {
            out << "Atomic Mass " << i << "  : "
                << setw(W) << SpeciesMass_[i] << "\n";
         }
         out << "Tref = " << setw(W) << Tref_ << "\n";
         //<< "PhiRef = " << setw(W) << PhiRef_ << "; "
         //<< "EntropyRef = " << setw(W) << EntropyRef_ << "; "
         //<< "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << "\n";
         out << "Potential Parameters : " << "\n";
         for (int i=0;i<NumberofSpecies_;++i)
         {
            for (int j=i;j<NumberofSpecies_;j++)
            {
               out << "[" << i << "][" << j << "] -- "
                   << setw(W) << *SpeciesPotential_[i][j] << "\n";
            }
         }
         out << "Normalization Modulus : " << setw(W) << NormModulus_ << "\n";
         out << "EulerAngles : " << setw(W) << EulerAng_[0]
             << setw(W) << EulerAng_[1] << setw(W) << EulerAng_[2] << "\n";
         out << "Loading Proportions : " << setw(W) << LoadingProportions_ << "\n";
         // also send to cout
         if (Echo_)
         {
            cout << "MultiLatticeTPP:" << "\n" << "\n";
            cout << "Density_ = " << Density_ << "\n";
            cout << "Using: " << (*CBK_) << " Kinematics" << "\n";
            cout << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
            for (int i=0;i<InternalAtoms_;++i)
            {
               cout << "Atom_" << i << "          "
                    << "Species : " <<setw(5) << AtomSpecies_[i]
                    << "          Position : " << setw(W) << CBK_->AtomPositions(i) << "\n";
            }
            cout << "Influence Distance   : " << setw(W) << InfluenceDist_ << "\n";
            for (int i=0;i<NumberofSpecies_;++i)
            {
               cout << "Atomic Mass " << i << "  : "
                    << setw(W) << SpeciesMass_[i] << "\n";
            }
            cout << "Tref = " << setw(W) << Tref_ << "\n";
            //<< "PhiRef = " << setw(W) << PhiRef_ << "; "
            //<< "EntropyRef = " << setw(W) << EntropyRef_ << "; "
            //<< "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << "\n";
            cout << "Potential Parameters : " << "\n";
            for (int i=0;i<NumberofSpecies_;++i)
            {
               for (int j=i;j<NumberofSpecies_;j++)
               {
                  cout << "[" << i << "][" << j << "] -- "
                       << setw(W) << *SpeciesPotential_[i][j] << "\n";
               }
            }
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
             << "Potential Value (Normalized):" << setw(W) << engy << "\n"
             << "Thermal Expansion:" << setw(W) << TE_static << "\n"
             << "Entropy:" << setw(W) << entropy << "\n"
             << "HeatCapacity:" << setw(W) << heatcapacity << "\n";
         for (int i=0;i<InternalAtoms_;++i)
         {
            out << "BodyForce Value " << i << " (Inf Normalized):"
                << setw(W) << BodyForce_[i] << "\n";
         }
         out << "Stress (Normalized):" << setw(W) << str_static << "\n"
             << "Stiffness (Normalized):" << setw(W) << stiff_static
            //<< "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):"  << setw(W) << EigenValues
            //<< "Bifurcation Info:" << setw(W) << MinEigVal
            //<< setw(W) << NoNegEigVal << "\n"
             << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n"<<setw(W)
             << TestFunctVals_static<< "\n"
             << "Bifurcation Info:" << setw(W) << NoNegTestFunctions << "\n"
             << "Condensed Moduli (Normalized):" << setw(W) << CondModuli_static
             << "CondEV Info:" << setw(W) << CondEV_static
             << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
             << "BlochWave Stability (GridSize=" << GridSize_ << "):"
             << setw(W) << BlochWaveStable << ", "
             << setw(W) << K_static << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << "\n"
                 << "Lambda (Normalized): " << setw(W) << Lambda_ << "\n"
                 << "DOF's :" << "\n" << setw(W) << CBK_->DOF() << "\n"
                 << "Potential Value (Normalized):" << setw(W) << engy << "\n"
                 << "Thermal Expansion:" << setw(W) << TE_static << "\n"
                 << "Entropy:" << setw(W) << entropy << "\n"
                 << "HeatCapacity:" << setw(W) << heatcapacity << "\n";
            for (int i=0;i<InternalAtoms_;++i)
            {
               cout << "BodyForce Value " << i << " (Inf Normalized):"
                    << setw(W) << BodyForce_[i] << "\n";
            }
            cout << "Stress (Normalized):" << setw(W) << str_static << "\n"
                 << "Stiffness (Normalized):" << setw(W) << stiff_static
               //<< "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):"  << setw(W) << EigenValues
               //<< "Bifurcation Info:" << setw(W) << MinEigVal
               //<< setw(W) << NoNegEigVal << "\n"
                 << "Eigenvalue Info (Rots->1,2,3; Trans->4,5,6):" << "\n"<<setw(W)
                 << TestFunctVals_static << "\n"
                 << "Bifurcation Info:" << setw(W) << NoNegTestFunctions << "\n"
                 << "Condensed Moduli (Normalized):" << setw(W) << CondModuli_static
                 << "CondEV Info:" << setw(W) << CondEV_static
                 << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << "\n"
                 << "BlochWave Stability (GridSize=" << GridSize_ << "):"
                 << setw(W) << BlochWaveStable << ", "
                 << setw(W) << K_static << "\n";
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

ostream& operator<<(ostream& out,MultiLatticeTPP& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}


//---------------------- Debug Mode Handler --------------------------


void MultiLatticeTPP::DebugMode()
{
   char *Commands[] = {
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
      "Potential_",
      "ConvexityDX_",
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
      "SetParameters"
   };
   int NOcommands=52;
   
   char response[LINELENGTH];
   char prompt[] = "Debug > ";
   int W=cout.width();
   
   cout << setw(0) << prompt;
   
   cin.getline(response,LINELENGTH);
   
   int indx;
   while (strcasecmp(response,"q") &&
          strcasecmp(response,"quit") &&
          strcasecmp(response,"exit"))
   {
      indx=0;
      if (!strcmp(response,Commands[indx++]))
         cout << "InternalAtoms_ = " << InternalAtoms_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "CBK_->DOFS() = " << CBK_->DOFS() << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "InfluenceDist_ = " << InfluenceDist_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "NTemp_ = " << NTemp_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
      {
         for (int i=0;i<CBK_->DOFS();++i)
            cout << "DOF_[" << i << "] = " << (CBK_->DOF())[i] << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
         cout << "RefLattice_= " << setw(W) << CBK_->RefLattice();
      else if (!strcmp(response,Commands[indx++]))
         cout << "Density_= " << Density_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "NormModulus_= " << NormModulus_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "Lambda_= " << Lambda_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
      {
         for (int i=0;i<InternalAtoms_;++i)
         {
            cout << "BodyForce_[" << i << "]= " << setw(W)
                 << BodyForce_[i] << "\n";
         }
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         for (int i=0;i<InternalAtoms_;++i)
         {
            cout << "AtomicMass_[" << i << "]= " << setw(W)
                 << AtomicMass_[i] << "\n";
         }
      }
      else if (!strcmp(response,Commands[indx++]))
         cout << "GridSize_= " << GridSize_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
      {
         for (int i=0;i<InternalAtoms_;++i)
            for (int j=i;j<InternalAtoms_;++j)
            {
               cout << "Potential_[" << i << "][" << j << "]= "
                    << setw(W) << Potential_[i][j] << "\n";
            }
      }
      else if (!strcmp(response,Commands[indx++]))
         cout << "ConvexityDX_= " << ConvexityDX_ << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "stress= " << setw(W) << stress();
      else if (!strcmp(response,Commands[indx++]))
         cout << "stiffness= " << setw(W) << stiffness();
      else if (!strcmp(response,Commands[indx++]))
         cout << "CondensedModuli= " << setw(W) << CondensedModuli();
      else if (!strcmp(response,Commands[indx++]))
      {
         Vector K(DIM3,0.0);
         int NoPTS;
         char prefix[LINELENGTH];
         int oldEcho_=Echo_;
         cout << "\tK > ";
         cin >> K;
         cin.sync(); // clear input
         cout << "\tNoPTS > ";
         cin >> NoPTS;
         cin.sync(); // clear input
         cout << "\tprefix > ";
         cin >> prefix;
         cin.sync(); // clear input
         Echo_=0;
         cout << "ReferenceDispersionCurves= ";
         ReferenceDispersionCurves(K,NoPTS,prefix,cout);
         Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         Vector K(DIM3,0.0);
         cout << "ReferenceBlochWave= " << ReferenceBlochWave(K) << "\t" << K << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "\tK > ";
         Vector K(DIM3,0.0);
         cin >> K;
         cin.sync(); // clear input
         cout << "ReferenceDynamicalStiffness= "
              << setw(W) << ReferenceDynamicalStiffness(K) << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         Vector DOF(CBK_->DOFS(),0.0);
         cout << "\tDOF > ";
         cin >> DOF;
         cin.sync(); // clear input
         SetDOF(DOF);
      }
      else if (!strcmp(response,Commands[indx++]))
         cout << "StressDT= " << setw(W) << StressDT();
      else if (!strcmp(response,Commands[indx++]))
         cout << "StiffnessDT= " << setw(W) << StiffnessDT();
      else if (!strcmp(response,Commands[indx++]))
      {
         double Temp;
         cout << "\tTemp > ";
         cin >> Temp;
         cin.sync(); // clear input
         SetTemp(Temp);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         double dist;
         cout << "\tInfluenceDist > ";
         cin >> dist;
         cin.sync(); // clear input
         SetInfluenceDist(dist);
      }
      else if (!strcmp(response,Commands[indx++]))
         cout << "energy= " << energy() << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "E0= " << E0() << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "E1= " << setw(W) << E1() << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "E2= " << setw(W) << E2() << "\n";
      else if (!strcmp(response,Commands[indx++]))
         cout << "E3= " << setw(W) << E3();
      else if (!strcmp(response,Commands[indx++]))
         cout << "E4= " << setw(W) << E4();
      else if (!strcmp(response,Commands[indx++]))
      {
         int GridSize;
         cout << "\tGridSize > ";
         cin >> GridSize;
         cin.sync(); // clear input
         SetGridSize(GridSize);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int oldEcho_=Echo_;
         int cutoff;
         cout << "\tcutoff > ";
         cin >> cutoff;
         cin.sync(); // clear input
         Echo_ = 0;
         cout << setw(W);
         NeighborDistances(cutoff,cout);
         Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int oldEcho_=Echo_;
         Echo_=0;
         cout << setw(W) << *this;
         Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int oldEcho_=Echo_;
         Echo_=0;
         cout << setw(W);
         Print(cout,PrintLong);
         Echo_=oldEcho_;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         double lambda;
         cout << "\tLambda > ";
         cin >> lambda;
         cin.sync(); // clear input
         SetLambda(lambda);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "StressDL= " << setw(W) << StressDL();
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "StiffnessDL= " << setw(W) << StiffnessDL();
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int iter;
         cout << "\titer > ";
         cin >> iter;
         cin.sync(); // clear input
         FindLatticeSpacing(iter);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int width;
         int oldEcho=Echo_;
         double epsilon;
         cout << "\tConsistencyEpsilon > ";
         cin >> epsilon;
         cout << "\tWidth > ";
         cin >> width;
         Echo_=0;
         ConsistencyCheck(epsilon,width,cout);
         Echo_=oldEcho;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "dbg_ = " << dbg_ << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         double Tol;
         int MaxItr;
         cout << "\tTolerence > ";
         cin >> Tol;
         cout << "\tMaxItr > ";
         cin >> MaxItr;
         RefineEqbm(Tol,MaxItr,&cout);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "EulerAng_ = "
              << setw(W) << EulerAng_[0]
              << setw(W) << EulerAng_[1]
              << setw(W) << EulerAng_[2] << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "Rotation_ = "
              << setw(W) << Rotation_ << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "Loading_ = "
              << setw(W) << Loading_ << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << setw(W);
         PrintCurrentCrystalParamaters(cout);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "ThermalExpansion = " << setw(W) << ThermalExpansion() << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "Entropy = " << setw(W) << Entropy() << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         cout << "HeatCapacity = " << setw(W) << HeatCapacity() << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int n;
         cout << "\tNoAtoms > ";
         cin >> n;
         cout << "P.Transpose()\n"
              << setw(20) << TranslationProjection1D(n).Transpose() << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int n,f;
         cout << "\tFsize > ";
         cin >> f;
         cout << "\tNoAtoms > ";
         cin >> n;
         cout << "P.Transpose()\n"
              << setw(20) << TranslationProjection3D(f,n).Transpose() << "\n";
      }
      else if (!strcmp(response,Commands[indx++]))
      {
         int no = SpeciesPotential_[0][0]->GetNoParameters();
         double *vals;
         int sze = no*((NumberofSpecies_+1)*(NumberofSpecies_)/2);
         vals = new double[sze];
         cout << "\tEnter new Parameter values > ";
         for (int i=0;i<sze;++i)
            cin >> vals[i];
         SetParameters(vals);
      }
      else if (!strcmp(response,"?") ||
               !strcasecmp(response,"help"))
      {
         cout << setiosflags(ios::left);
         for (int i=0;i<NOcommands/2 + NOcommands%2;++i)
         {
            cout << "  " << setw(30) << Commands[i];
            if (i==NOcommands/2 && !NOcommands%2)
               cout << "\n";
            else
               cout << setw(30) << Commands[NOcommands/2+i] << "\n";
            
            if (!((i+1)%30))
            {
               cout << "more...." << "\n";
               char ans;
               cin.sync(); // clear input
               ans=kbhitWait();
               if (ans=='q') break;
            }
         }
         cout << resetiosflags(ios::left) << "\n";
      }
      else if (!strcmp(response,"\n") ||
               !strcmp(response,""))
      {
      }
      else
      {
         cout << "!--- Error - Unknown command ---!" << "\n" << "\n";
      }
      
      cout << "\n" << prompt;
      cin.getline(response,LINELENGTH);
   }
}


void MultiLatticeTPP::RefineEqbm(double const& Tol,int const& MaxItr,ostream* const out)
{
   Vector dx(CBK_->DOFS(),0.0);
   Vector Stress=E1();
   int itr=0;
   while ((itr < MaxItr) && Stress.Norm() > Tol)
   {
      ++itr;

#ifdef SOLVE_SVD
      dx = SolveSVD(E2(),Stress,MAXCONDITION,Echo_);
#else
      dx = SolvePLU(E2(),Stress);
#endif

      SetDOF(CBK_->DOF()-dx);
      Stress=E1();
      if (out != 0)
      {
         *out << setw(20) << Stress;
         
         *out << itr << "\tdx " << dx.Norm() << "\tstress " << Stress.Norm() << "\n";
      }
   }
}

void MultiLatticeTPP::PrintCurrentCrystalParamaters(ostream& out) const
{
   Matrix F(DIM3,DIM3,0.0);
   Vector CurrentLattice[DIM3];
   int W=out.width();
   out.width(0);
   
   for (int i=0;i<DIM3;++i)
   {
      CurrentLattice[i].Resize(DIM3);
      CurrentLattice[i] = CBK_->CurrentLatticeVec(i);
   }
   
   
   out << "TITLE LatticeStatics crystal structure scaled by 10.0" << "\n";
   out << "DIMENSION 3" << "\n";
   out << "CELL" << setw(W) << 10.0*CurrentLattice[0].Norm()
       << setw(W) << 10.0*CurrentLattice[1].Norm()
       << setw(W) << 10.0*CurrentLattice[2].Norm();
   
   double alpha,beta,gamma;
   double pi = 4.0*atan(1.0);
   
   alpha = acos(CurrentLattice[1]*CurrentLattice[2]
                /(CurrentLattice[1].Norm()*CurrentLattice[2].Norm()));
   beta  = acos(CurrentLattice[0]*CurrentLattice[2]
                /(CurrentLattice[0].Norm()*CurrentLattice[2].Norm()));
   gamma = acos(CurrentLattice[0]*CurrentLattice[1]
                /(CurrentLattice[0].Norm()*CurrentLattice[1].Norm()));
   
   out << setw(W) << alpha*180.0/pi
       << setw(W) << beta*180.0/pi
       << setw(W) << gamma*180.0/pi << "\n";
   
   out << "SYMMETRY  NUMBER 1  LABEL P1  " << "\n";
   out << "SYM MAT  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0 0.0000 0.0000 0.0000" << "\n";
   out << "\n" << "ATOMS" << "\n"
       << "NAME" << setw(W) << "X" << setw(W) << "Y" << setw(W) << "Z" << "\n";
   char const *species[] = {"Ni","Ti","C"};
   out << setw(4) << species[(AtomSpecies_[0] > 3)?3:AtomSpecies_[0]]
       << setw(W) << CBK_->FractionalPosVec(0) << "\n";
   for (int i=1;i<InternalAtoms_;++i)
   {
      out << setw(4) << species[(AtomSpecies_[i] > 3)?3:AtomSpecies_[i]];
      out << setw(W) << CBK_->FractionalPosVec(i) << "\n";
   }
   
   out << "EOF" << "\n";
   
   out << "\n"
       << "Temperature : " << setw(W) << NTemp_ << "\n"
       << "Lambda : " << setw(W) << Lambda_ << "\n"
       << "DOFs : " << setw(W) << CBK_->DOF() << "\n";
}
