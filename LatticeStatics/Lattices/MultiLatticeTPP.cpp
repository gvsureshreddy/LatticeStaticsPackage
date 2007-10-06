#include "MultiLatticeTPP.h"
#include <cmath>

#include "UtilityFunctions.h"

using namespace std;

MultiLatticeTPP::~MultiLatticeTPP()
{
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

MultiLatticeTPP::MultiLatticeTPP(char *datafile,const char *prefix,int Echo,int Width,int Debug)
{
   Echo_ = Echo;
   dbg_ = Debug;
   // Get Lattice definition
   char tmp[LINELENGTH];
   if(!GetParameter(prefix,"InternalAtoms",datafile,"%u",&INTERNAL_ATOMS)) exit(-1);
   
   // Initiate the CBK object (default to SymLagrangeCB)
   const char *CBKin[] = {"SymLagrangeCB","SymMixedCB","SymEulerCB",
			  "LagrangeCB","MixedCB","EulerCB"};
   switch (GetStringParameter(prefix,"CBKinematics",datafile,CBKin,6,0))
   {
      case 5:
	 CBK_ = new EulerCB(INTERNAL_ATOMS,prefix,datafile);
	 break;
      case 4:
	 CBK_ = new MixedCB(INTERNAL_ATOMS,prefix,datafile);
	 break;
      case 3:
	 CBK_ = new LagrangeCB(INTERNAL_ATOMS,prefix,datafile);
	 break;
      case 2:
	 CBK_ = new SymEulerCB(INTERNAL_ATOMS,prefix,datafile);
	 break;
      case 1:
	 CBK_ = new SymMixedCB(INTERNAL_ATOMS,prefix,datafile);
	 break;
      case 0:
      default:
	 CBK_ = new SymLagrangeCB(INTERNAL_ATOMS,prefix,datafile);
   }
   
   if (DOFMAX < CBK_->DOFS())
   {
      cerr << "Error (MultiLatticeTPP()): DOFMAX < " << CBK_->DOFS() << " in Lattice.h" << endl;
      exit(-5);
   }
   
   // Setup Bodyforce_
   BodyForce_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
      BodyForce_[i].Resize(DIM3,0.0);

   // Get Thermo parameters
   if (!GetParameter(prefix,"Tref",datafile,"%lf",&Tref_)) exit(-1);
   //if (!GetParameter(prefix,"PhiRef",datafile,"%lf",&PhiRef_)) exit(-1);
   //if (!GetParameter(prefix,"EntropyRef",datafile,"%lf",&EntropyRef_)) exit(-1);
   //if (!GetParameter(prefix,"HeatCapacityRef",datafile,"%lf",&HeatCapacityRef_)) exit(-1);

   if (!GetIntVectorParameter(prefix,"AtomSpecies",datafile,INTERNAL_ATOMS,AtomSpecies_))
      exit(-1);
   NumberofSpecies_ = AtomSpecies_[0];
   for (int i=1;i<INTERNAL_ATOMS;++i)
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
   Potential_ = new PairPotentials**[INTERNAL_ATOMS];
   Potential_[0] = new PairPotentials*[INTERNAL_ATOMS*INTERNAL_ATOMS];
   for (int i=1;i<INTERNAL_ATOMS;++i)
   {
      Potential_[i] = Potential_[i-1] + INTERNAL_ATOMS;
   }

   SpeciesMass_ = new double[NumberofSpecies_];
   AtomicMass_ = new double[INTERNAL_ATOMS];

   for (int i=0;i<NumberofSpecies_;++i)
   {
      for (int j=i;j<NumberofSpecies_;++j)
      {
	 SpeciesPotential_[i][j] = SpeciesPotential_[j][i]
	    = InitializePairPotential(datafile,prefix,i,j);
      }
      sprintf(tmp,"AtomicMass_%u",i);
      if (!GetParameter(prefix,tmp,datafile,"%lf",&(SpeciesMass_[i]))) exit(-1);
   }
   
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      for (int j=i;j<INTERNAL_ATOMS;++j)
      {
	 Potential_[i][j] = Potential_[j][i]
	    = SpeciesPotential_[AtomSpecies_[i]][AtomSpecies_[j]];
      }

      AtomicMass_[i] = SpeciesMass_[AtomSpecies_[i]];
   }
	 
   // Get Lattice parameters
   NTemp_ = 1.0;
   if(!GetParameter(prefix,"InfluanceDist",datafile,"%lf",&InfluenceDist_)) exit(-1);
   if(!GetParameter(prefix,"NormModulus",datafile,"%lf",&NormModulus_)) exit(-1);
   if(!GetParameter(prefix,"ConvexityDX",datafile,"%lf",&ConvexityDX_)) exit(-1);

   // Get Loading parameters
   const char *loadparams[] = {"Temperature","Load"};
   int NoParams=2;
   switch (GetStringParameter(prefix,"LoadingParameter",datafile,loadparams,NoParams))
   {
      case 0: LoadParameter_ = Temperature; break;
      case 1: LoadParameter_ = Load; break;
      case -1: cerr << "Unknown Loading Parameter" << endl; exit(-1); break;
   }
   Lambda_ = 0.0;
   if(!GetParameter(prefix,"EulerAngle_X",datafile,"%lf",&(EulerAng_[0]))) exit(-1);
   if(!GetParameter(prefix,"EulerAngle_Y",datafile,"%lf",&(EulerAng_[1]))) exit(-1);
   if(!GetParameter(prefix,"EulerAngle_Z",datafile,"%lf",&(EulerAng_[2]))) exit(-1);
   LoadingProportions_.Resize(DIM3);
   if(!GetVectorParameter(prefix,"LoadProportions",datafile,&LoadingProportions_)) exit(-1);
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
   if(!GetParameter(prefix,"MaxIterations",datafile,"%u",&iter)) exit(-1);
   if(!GetParameter(prefix,"BlochWaveGridSize",datafile,"%u",&GridSize_)) exit(-1);

   // Initiate the Lattice Sum object
   LatSum_(CBK_,INTERNAL_ATOMS,Potential_,&InfluenceDist_,&NTemp_);

   int err=0;
   err=FindLatticeSpacing(datafile,prefix,iter);
   if (err)
   {
      cerr << "unable to find initial lattice spacing!" << endl;
      exit(-1);
   }

   // Setup initial status for parameters
   if(!GetParameter(prefix,"NTemp",datafile,"%lf",&NTemp_)) exit(-1);
   if(!GetParameter(prefix,"Lambda",datafile,"%lf",&Lambda_)) exit(-1);
   // Make any changes to atomic potentials that might be required
   strcpy(tmp,prefix); strcat(tmp,"Update-");
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      for (int j=i;j<INTERNAL_ATOMS;++j)
      {
	 if (AtomSpecies_[i] < AtomSpecies_[j])
	    UpdatePairPotential(datafile,tmp,AtomSpecies_[i],AtomSpecies_[j],Potential_[i][j]);
	 else
	    UpdatePairPotential(datafile,tmp,AtomSpecies_[j],AtomSpecies_[i],Potential_[j][i]);
      }
   }
   LatSum_.Recalc();

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   UCIter_(GridSize_);
}

int MultiLatticeTPP::FindLatticeSpacing(char *datafile,const char *prefix,int iter)
{
   Lambda_=0.0;
   NTemp_=1.0;

   CBK_->SetReferenceDOFs();
   LatSum_.Recalc();

   if (Echo_)
      RefineEqbm(1.0e-13,iter,&cout);
   else
      RefineEqbm(1.0e-13,iter,NULL);

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

// Lattice Routines
double MultiLatticeTPP::energy(PairPotentials::TDeriv dt)
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
   Vr = CBK_->RefVolume();
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
      //	      (NTemp_*Tref_)*EntropyRef_ -
      //	      HeatCapacityRef_*(NTemp_*Tref_)*(log(NTemp_*Tref_) - 1.0)
      //	 )/NormModulus_;
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
      cerr << "Error in MultiLatticeTPP::energy" << endl;
      exit(-1);
   }
   
   return Phi;
}

Matrix MultiLatticeTPP::stress(PairPotentials::TDeriv dt,LDeriv dl)
{
   static Matrix S;
   double ForceNorm = 0.0;
   double phi,Vr;
   int i,j;

   S.Resize(1,CBK_->DOFS(),0.0);

   Vr = CBK_->RefVolume();
   
   if (dl==L0)
   {
      for (i=0;i<INTERNAL_ATOMS;++i)
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
	    cerr << "Error in MultiLatticeTPP::stress" << endl;
	    exit(-1);
	 }
	 
	 for (i=0;i<DIM3;i++)
	 {
	    for (j=0;j<DIM3;j++)
	    {
	       S[0][CBK_->INDF(i,j)] += phi*CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j);
	    }
	 }
	 for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 {
	    for (j=0;j<DIM3;j++)
	    {
	       S[0][CBK_->INDS(i,j)] += phi*CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),
						LatSum_.Atom(1),i,j);
	    }
	 }
      }

      // BodyForce[i] = BodyForce[i] / ForceNorm
      for (i=0;i<INTERNAL_ATOMS;i++)
      {
	 for (j=0;j<DIM3;j++)
	 {
	    BodyForce_[i][j] /= ForceNorm;
	 }
      }

      // S = S/(2*Vr*NormModulus)
      S *= 1.0/(2.0*(Vr*NormModulus_));

      // Load terms
      if (dt == PairPotentials::T0)
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       S[0][CBK_->INDF(i,j)] -= Lambda_*Loading_[j][i];
	    }
      }

   }
   else if (dl==DL)
   {
      // dl=DL
      for (i=0;i<DIM3;++i)
	 for (j=0;j<DIM3;++j)
	    S[0][CBK_->INDF(i,j)] -= Loading_[j][i];
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stress()" << endl;
      exit(-1);
   }
   
   return S;
}
      
Matrix MultiLatticeTPP::stiffness(PairPotentials::TDeriv dt,LDeriv dl)
{
   static Matrix Phi;
   Matrix F(DIM3,DIM3);
   double phi,phi1;
   int i,j,k,l;

   Phi.Resize(CBK_->DOFS(),CBK_->DOFS(),0.0);

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
	    cerr << "Error in MultiLatticeTPP::stiffness" << endl;
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
		     Phi[CBK_->INDF(i,j)][CBK_->INDF(k,l)]+=
			phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
			     *CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),k,l))
			+phi1*CBK_->D2yDFF(LatSum_.pDX(),i,j,k,l);
		  }
	       }
	    }
	 }
	 //Lower Diag Block (CBK_->Ssize(),CBK_->Ssize())
	 for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 {
	    for (j=0;j<DIM3;j++)
	    {
	       for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       {
		  for (l=0;l<DIM3;l++)
		  {
		     Phi[CBK_->INDS(i,j)][CBK_->INDS(k,l)]+=
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
	       for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       {
		  for (l=0;l<DIM3;l++)
		  {
		     Phi[CBK_->INDF(i,j)][CBK_->INDS(k,l)] =
			Phi[CBK_->INDS(k,l)][CBK_->INDF(i,j)] +=
			phi*(CBK_->DyDF(LatSum_.pDx(),LatSum_.pDX(),i,j)
			     *CBK_->DyDS(LatSum_.pDx(),LatSum_.Atom(0),LatSum_.Atom(1),k,l))
			+phi1*CBK_->D2yDFS(LatSum_.pDx(),LatSum_.pDX(),
					  LatSum_.Atom(0),LatSum_.Atom(1),i,j,k,l);
		  }
	       }
	    }
	 }
      }
      
      // Phi = Phi/(2*Vr*NormModulus)
      Phi *= 1.0/(2.0*(CBK_->RefVolume()*NormModulus_));
   }
   else if (dl==DL)
   {
      // Nothing to do: Phi is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiLatticeTpp::stiffness()" << endl;
      exit(-1);
   }   
   return Phi;
}

Matrix MultiLatticeTPP::E3()
{
   static Matrix Phi;
   double phi,phi1,phi2;
   int i,j,k,l,m,n;

   Phi.Resize(CBK_->DOFS()*CBK_->DOFS(),CBK_->DOFS(),0.0);

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
			Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDF(m,n)] +=
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
      for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 for (j=0;j<DIM3;j++)
	    for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       for (l=0;l<DIM3;l++)
		  for (m=CBK_->NoTrans();m<INTERNAL_ATOMS;m++)
		     for (n=0;n<DIM3;n++)
		     {
			Phi[CBK_->INDSS(i,j,k,l)][CBK_->INDS(m,n)] +=
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
		  for (m=CBK_->NoTrans();m<INTERNAL_ATOMS;m++)
		     for (n=0;n<DIM3;n++)
		     {
			Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDS(m,n)] =
			   Phi[CBK_->INDFS(i,j,m,n)][CBK_->INDF(k,l)] =
			   Phi[CBK_->INDSF(m,n,i,j)][CBK_->INDF(k,l)] += (
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
      for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 for (j=0;j<DIM3;j++)
	    for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       for (l=0;l<DIM3;l++)
		  for (m=0;m<DIM3;m++)
		     for (n=0;n<DIM3;n++)
		     {
			Phi[CBK_->INDSS(i,j,k,l)][CBK_->INDF(m,n)] =
			   Phi[CBK_->INDSF(i,j,m,n)][CBK_->INDS(k,l)] =
			   Phi[CBK_->INDFS(m,n,i,j)][CBK_->INDS(k,l)] += (
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

   
   // Phi = Phi/(2*Vr*NormModulus)
   Phi *= 1.0/(2.0*(CBK_->RefVolume()*NormModulus_));

   return Phi;
}

Matrix MultiLatticeTPP::E4()
{
   static Matrix Phi;
   double phi,phi1,phi2,phi3;
   int i,j,k,l,m,n,s,t;

   Phi.Resize(CBK_->DOFS()*CBK_->DOFS(),CBK_->DOFS()*CBK_->DOFS(),0.0);

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
			      Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDFF(m,n,s,t)]+=
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
      for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 for (j=0;j<DIM3;j++)
	    for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       for (l=0;l<DIM3;l++)
		  for (m=CBK_->NoTrans();m<INTERNAL_ATOMS;m++)
		     for (n=0;n<DIM3;n++)
			for (s=CBK_->NoTrans();s<INTERNAL_ATOMS;s++)
			   for (t=0;t<DIM3;t++)
			   {
			      Phi[CBK_->INDSS(i,j,k,l)][CBK_->INDSS(m,n,s,t)] +=
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
			for (s=CBK_->NoTrans();s<INTERNAL_ATOMS;s++)
			   for (t=0;t<DIM3;t++)
			   {
			      Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDFS(m,n,s,t)] =
				 Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDSF(s,t,m,n)] =
				 Phi[CBK_->INDFS(i,j,s,t)][CBK_->INDFF(k,l,m,n)] =
				 Phi[CBK_->INDSF(s,t,i,j)][CBK_->INDFF(k,l,m,n)] += (
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
      for (i=CBK_->NoTrans();i<INTERNAL_ATOMS;i++)
	 for (j=0;j<DIM3;j++)
	    for (k=CBK_->NoTrans();k<INTERNAL_ATOMS;k++)
	       for (l=0;l<DIM3;l++)
		  for (m=CBK_->NoTrans();m<INTERNAL_ATOMS;m++)
		     for (n=0;n<DIM3;n++)
			for (s=0;s<DIM3;s++)
			   for (t=0;t<DIM3;t++)
			   {
			      Phi[CBK_->INDSS(i,j,k,l)][CBK_->INDSF(m,n,s,t)] =
				 Phi[CBK_->INDSS(i,j,k,l)][CBK_->INDFS(s,t,m,n)] =
				 Phi[CBK_->INDSF(i,j,s,t)][CBK_->INDSS(k,l,m,n)] =
				 Phi[CBK_->INDFS(s,t,i,j)][CBK_->INDSS(k,l,m,n)] += (
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
		  for (m=CBK_->NoTrans();m<INTERNAL_ATOMS;m++)
		     for (n=0;n<DIM3;n++)
			for (s=CBK_->NoTrans();s<INTERNAL_ATOMS;s++)
			   for (t=0;t<DIM3;t++)
			   {
			      Phi[CBK_->INDFF(i,j,k,l)][CBK_->INDSS(m,n,s,t)] =
				 Phi[CBK_->INDFS(i,j,m,n)][CBK_->INDFS(k,l,s,t)] =
				 Phi[CBK_->INDFS(i,j,m,n)][CBK_->INDSF(s,t,k,l)] =
				 Phi[CBK_->INDSF(m,n,i,j)][CBK_->INDFS(k,l,s,t)] =
				 Phi[CBK_->INDSF(m,n,i,j)][CBK_->INDSF(s,t,k,l)] =
				 Phi[CBK_->INDSS(m,n,s,t)][CBK_->INDFF(i,j,k,l)] += (
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

   
   // Phi = Phi/(2*Vr*NormModulus)
   Phi *= 1.0/(2.0*(CBK_->RefVolume()*NormModulus_));

   return Phi;
}

Matrix MultiLatticeTPP::CondensedModuli()
{
   Matrix stiff = stiffness();
   int intrn = CBK_->Ssize();
   double factor = 1.0/(intrn/DIM3);
   int fsz = CBK_->Fsize();
   Matrix CM(fsz,fsz), IM(intrn,intrn);
   
   for (int i=0;i<fsz;i++)
      for (int j=0;j<fsz;j++)
      {
	 CM[i][j] = stiff[i][j];
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
		  CM[i][j] -= stiff[i][fsz+m]*IM[m][n]*stiff[fsz+n][j];
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
	    CM[i][j] /= 2.0;
	    CM[j][i] /= 2.0;
	 }
	 
	 for (int j=3;j<6;j++)
	 {
	    CM[i][j] /= 4.0;
	 }
      }
   }

   return CM;
}

int MultiLatticeTPP::comp(const void *a,const void *b)
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

int MultiLatticeTPP::abscomp(const void *a,const void *b)
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

void MultiLatticeTPP::interpolate(Matrix *EigVals,int zero,int one,int two)
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

CMatrix MultiLatticeTPP::ReferenceDynamicalStiffness(Vector &K)
{
   static CMatrix Dk;
   static double pi = 4.0*atan(1.0);
   static MyComplexDouble Ic(0,1);
   static MyComplexDouble A = 2.0*pi*Ic;
   int i,j;

   Dk.Resize(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3,0.0);
   
   for (LatSum_.Reset();!LatSum_.Done();++LatSum_)
   {
      // Calculate Dk
      if (LatSum_.Atom(0) != LatSum_.Atom(1))
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       // y != y' terms (i.e., off block (3x3) diagonal terms)
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
		  (-2.0*Del(i,j)*LatSum_.phi1()
		   -4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2())
		  *exp(A*
		       (K[0]*LatSum_.DX(0) + K[1]*LatSum_.DX(1) + K[2]*LatSum_.DX(2)));

	       // y==y' components (i.e., Phi(0,y,y) term)
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(0)+j] +=
		  (2.0*Del(i,j)*LatSum_.phi1()
		   +4.0*LatSum_.Dx(i)*LatSum_.Dx(j)*LatSum_.phi2());
	    }
      }
      else
      {
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       Dk[DIM3*LatSum_.Atom(0)+i][DIM3*LatSum_.Atom(1)+j] +=
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
   for (int p=0;p<INTERNAL_ATOMS;++p)
      for (int q=0;q<INTERNAL_ATOMS;++q)
	 for (i=0;i<DIM3;++i)
	    for (j=0;j<DIM3;++j)
	    {
	       Dk[DIM3*p+i][DIM3*q+j] /= sqrt(AtomicMass_[p]*AtomicMass_[q]);
	    }
   
   return Dk;
}

void MultiLatticeTPP::ReferenceDispersionCurves(Vector K,int NoPTS,const char *prefix,
						ostream &out)
{
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);

   Matrix InverseLat(DIM3,DIM3);
   InverseLat = (CBK_->RefLattice()).Inverse();

   Matrix EigVal[3];
   for (int i=0;i<3;++i) EigVal[i].Resize(1,INTERNAL_ATOMS*DIM3);

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
      qsort(EigVal[k][0],INTERNAL_ATOMS*DIM3,sizeof(double),&comp);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 out << setw(w) << EigVal[k][0][i];
	 if (Echo_) cout << setw(w) << EigVal[k][0][i];
      }
      out << endl;
      if (Echo_) cout << endl;
   }
   int zero=0,one=1,two=2;
   for (int k=2;k<NoPTS;++k)
   {
      Z = Z1 + (k*dz)*DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[two][0],INTERNAL_ATOMS*DIM3,sizeof(double),&comp);
      interpolate(EigVal,zero,one,two);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 out << setw(w) << EigVal[two][0][i];;
	 if (Echo_) cout << setw(w) << EigVal[two][0][i];;
      }
      out << endl;
      if (Echo_) cout << endl;

      zero = (++zero)%3; one = (zero+1)%3; two = (one+1)%3;
   }
}

int MultiLatticeTPP::ReferenceBlochWave(Vector &K)
{
   static CMatrix A(INTERNAL_ATOMS*DIM3,INTERNAL_ATOMS*DIM3);
   static Matrix EigVals(1,INTERNAL_ATOMS*DIM3);
   static Matrix InverseLat(DIM3,DIM3);
   static Vector Z(DIM3);

   InverseLat = (CBK_->RefLattice()).Inverse();

   // Iterate over points in cubic unit cell
   for (UCIter_.Reset();!UCIter_.Done();++UCIter_)
   {
      for (int i=0;i<DIM3;++i)
      {
	 K[i] = UCIter_[i];
      }

      Z = InverseLat*K;
      A = ReferenceDynamicalStiffness(Z);

      EigVals = HermiteEigVal(A);
      
      for (int i=0;i<INTERNAL_ATOMS*DIM3;++i)
      {
	 // if w^2 <= 0.0 --> Re(i*w*x) > 0 --> growing solutions --> unstable
	 if ( EigVals[0][i] <= 0.0 )
	 {
	    return 0;
	 }
      }
   }
   return 1;
}

//---- needs to be updated----//
void MultiLatticeTPP::LongWavelengthModuli(double dk, int gridsize,const char *prefix,
					   ostream &out)
{
   static double pi = 4*atan(1.0);
   static double twopi = 2*pi;
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
      Dpp((INTERNAL_ATOMS-1)*3,(INTERNAL_ATOMS-1)*3,0.0),
      Dfp((INTERNAL_ATOMS-1)*3,9,0.0);
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
      for (int i=1;i<INTERNAL_ATOMS;++i)
      {
	 for (int j=0;j<DIM3;++j)
	 {
	    for (int k=1;k<INTERNAL_ATOMS;++k)
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
      for (int k=1;k<INTERNAL_ATOMS;++k)
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
   Phi *= 1.0/(2.0*(CBK_->RefVolume()*NormModulus_));
   
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

   cout << "Condensed Moduli Check!:" << endl;
   cout << setw(w) << Lp-SymPhi << endl;

   

   Vector K(DIM3),Z(DIM3,0.0);
   Matrix BlkEigVal(1,INTERNAL_ATOMS*DIM3);
   Matrix ModEigVal(1,DIM3);

   double Mc=0.0;
   double Vc=CBK_->RefVolume();
   for (int i=0;i<INTERNAL_ATOMS;++i)
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
	 qsort(BlkEigVal[0],INTERNAL_ATOMS*DIM3,sizeof(double),&abscomp);

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
	 out << endl;
	 if (Echo_) cout << endl;
      }
      out << endl;
      if (Echo_) cout << endl;
   }
}

void MultiLatticeTPP::NeighborDistances(int cutoff,ostream &out)
{
   Matrix NeighborDist =
      LatSum_.NeighborDistances(cutoff,pow(double(10),double(-(out.precision()-1))));
   
   int W=out.width();
   int types = (INTERNAL_ATOMS*(INTERNAL_ATOMS+1))/2;
   for (int i=0;i<cutoff;++i)
   {
      out << setw(W) << NTemp_ << setw(W) << NeighborDist[i][0];
      for (int j=0;j<types;++j)
      {
	 out << setw(W/4) << int(NeighborDist[i][1+j]);
      }
      out << endl;
   }
   out << endl;
}

void MultiLatticeTPP::Print(ostream &out,PrintDetail flag)
{
   static int W;
   static int NoNegEigVal;
   static double MinEigVal;
   static double engy,entropy,heatcapacity;
   static Matrix
      str(1,CBK_->DOFS()),
      stiff(CBK_->DOFS(),CBK_->DOFS()),
      EigenValues(1,CBK_->DOFS()),
      CondEV(1,CBK_->Fsize());
   static Matrix
      CondModuli(CBK_->Fsize(),CBK_->Fsize());
   static int RankOneConvex;
   static Vector K(DIM3);
   static int BlochWaveStable;
   
   W=out.width();

   out.width(0);
   if (Echo_) cout.width(0);

   NoNegEigVal=0;

   engy = energy();
   entropy = Entropy();
   heatcapacity = HeatCapacity();
   str = stress();
   stiff = stiffness();
   
   EigenValues=SymEigVal(stiff);
   MinEigVal = EigenValues[0][0];
   for (int i=0;i<CBK_->DOFS();i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   CondModuli = CondensedModuli();

   CondEV=SymEigVal(CondModuli);
   RankOneConvex = FullScanRank1Convex3D(CondModuli,ConvexityDX_);

   K.Resize(DIM3,0.0);
   if (RankOneConvex)
   {
      BlochWaveStable = BlochWave(K);
   }
   else
   {
      BlochWaveStable = -1;
   }


   switch (flag)
   {
      case PrintLong:
	 out << "MultiLatticeTPP:" << endl << endl;
	 out << "Using: " << (*CBK_) << " Kinematics" << endl;
	 out << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "Atom_" << i << "          "
		<< "Species : " << setw(5) << AtomSpecies_[i]
		<< "          Position : " << setw(W) << CBK_->AtomPositions(i) << endl;
	 }
	 out << "Influence Distance   : " << setw(W) << InfluenceDist_ << endl;
	 for (int i=0;i<NumberofSpecies_;++i)
	 {
	    out << "Atomic Mass " << i << "  : "
		<< setw(W) << SpeciesMass_[i] << endl;
	 }
	 out << "Tref = " << setw(W) << Tref_ << endl;
	 //<< "PhiRef = " << setw(W) << PhiRef_ << "; "
	 //<< "EntropyRef = " << setw(W) << EntropyRef_ << "; "
	 //<< "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << endl;
	 out << "Potential Parameters : " << endl;
	 for (int i=0;i<NumberofSpecies_;++i)
	 {
	    for (int j=i;j<NumberofSpecies_;j++)
	    {
	       out << "[" << i << "][" << j << "] -- "
		   << setw(W) << SpeciesPotential_[i][j] << endl;
	    }
	 }
	 out << "Normalization Modulus : " << setw(W) << NormModulus_ << endl;
	 out << "EulerAngles : " << setw(W) << EulerAng_[0]
	     << setw(W) << EulerAng_[1] << setw(W) << EulerAng_[2] << endl;
	 out << "Loading Proportions : " << setw(W) << LoadingProportions_ << endl;
	 // also send to cout
	 if (Echo_)
	 {
	    cout << "MultiLatticeTPP:" << endl << endl;
	    cout << "Using: " << (*CBK_) << " Kinematics" << endl;
	    cout << "RefLattice_ : " << setw(W) << CBK_->RefLattice();
	    for (int i=0;i<INTERNAL_ATOMS;++i)
	    {
	       cout << "Atom_" << i << "          "
		    << "Species : " <<setw(5) << AtomSpecies_[i]
		    << "          Position : " << setw(W) << CBK_->AtomPositions(i) << endl;
	    }
	    cout << "Influence Distance   : " << setw(W) << InfluenceDist_ << endl;
	    for (int i=0;i<NumberofSpecies_;++i)
	    {
	       cout << "Atomic Mass " << i << "  : "
		    << setw(W) << SpeciesMass_[i] << endl;
	    }
	    cout << "Tref = " << setw(W) << Tref_ << endl;
	    //<< "PhiRef = " << setw(W) << PhiRef_ << "; "
	    //<< "EntropyRef = " << setw(W) << EntropyRef_ << "; "
	    //<< "HeatCapacityRef = " << setw(W) << HeatCapacityRef_ << endl;
	    cout << "Potential Parameters : " << endl;
	    for (int i=0;i<NumberofSpecies_;++i)
	    {
	       for (int j=i;j<NumberofSpecies_;j++)
	       {
		  cout << "[" << i << "][" << j << "] -- "
		       << setw(W) << SpeciesPotential_[i][j] << endl;
	       }
	    }
	    cout << "Normalization Modulus : " << setw(W) << NormModulus_ << endl;
	    cout << "EulerAngles : " << setw(W) << EulerAng_[0]
		 << setw(W) << EulerAng_[1] << setw(W) << EulerAng_[2] << endl;
	    cout << "Loading Proportions : " << setw(W) << LoadingProportions_ << endl;
	 }
	 // passthrough to short
      case PrintShort:
	 out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	     << "Lambda (Normalized): " << setw(W) << Lambda_ << endl
	     << "DOF's :" << endl << setw(W) << CBK_->DOF() << endl
	     << "Potential Value (Normalized):" << setw(W) << engy << endl
	     << "Entropy:" << setw(W) << entropy << endl
	     << "HeatCapacity:" << setw(W) << heatcapacity << endl;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "BodyForce Value " << i << " (Inf Normalized):"
		<< setw(W) << BodyForce_[i] << endl;
	 }
	 out << "Stress (Normalized):" << setw(W) << str << endl
	     << "Stiffness (Normalized):" << setw(W) << stiff
	     << "Eigenvalue Info:"  << setw(W) << EigenValues
	     << "Bifurcation Info:" << setw(W) << MinEigVal
	     << setw(W) << NoNegEigVal << endl
	     << "Condensed Moduli (Normalized):" << setw(W) << CondModuli
	     << "CondEV Info:" << setw(W) << CondEV
	     << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl
	     << "BlochWave Stability (GridSize=" << GridSize_ << "):"
	     << setw(W) << BlochWaveStable << ", "
	     << setw(W) << K << endl;
	 // send to cout also
	 if (Echo_)
	 {
	    cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
		 << "Lambda (Normalized): " << setw(W) << Lambda_ << endl
		 << "DOF's :" << endl << setw(W) << CBK_->DOF() << endl
		 << "Potential Value (Normalized):" << setw(W) << engy << endl
		 << "Entropy:" << setw(W) << entropy << endl
		 << "HeatCapacity:" << setw(W) << heatcapacity << endl;
	    for (int i=0;i<INTERNAL_ATOMS;++i)
	    {
	       cout << "BodyForce Value " << i << " (Inf Normalized):"
		    << setw(W) << BodyForce_[i] << endl;
	    }
	    cout << "Stress (Normalized):" << setw(W) << str << endl
		 << "Stiffness (Normalized):" << setw(W) << stiff
		 << "Eigenvalue Info:"  << setw(W) << EigenValues
		 << "Bifurcation Info:" << setw(W) << MinEigVal
		 << setw(W) << NoNegEigVal << endl
		 << "Condensed Moduli (Normalized):" << setw(W) << CondModuli
		 << "CondEV Info:" << setw(W) << CondEV
		 << "Condensed Moduli Rank1Convex:" << setw(W) << RankOneConvex << endl
		 << "BlochWave Stability (GridSize=" << GridSize_ << "):"
		 << setw(W) << BlochWaveStable << ", "
		 << setw(W) << K << endl;
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

ostream &operator<<(ostream &out,MultiLatticeTPP &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}


//---------------------- Debug Mode Handler --------------------------


void MultiLatticeTPP::DebugMode()
{
   char *Commands[] = {
      "INTERNAL_ATOMS",
      "DOFS",
      "InfluenceDist_",
      "NTemp_",
      "DOF_",
      "RefLattice_",
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
      "Entropy",
      "TranslationProjection1D",
      "TranslationProjection3D",
   };
   int NOcommands=48;
   
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
	 cout << "INTERNAL_ATOMS = " << INTERNAL_ATOMS << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "CBK_->DOFS() = " << CBK_->DOFS() << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "InfluenceDist_ = " << InfluenceDist_ << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "NTemp_ = " << NTemp_ << endl;
      else if (!strcmp(response,Commands[indx++]))
      {
	 for (int i=0;i<CBK_->DOFS();++i)
	    cout << "DOF_[" << i << "] = " << (CBK_->DOF())[i] << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
	 cout << "RefLattice_= " << setw(W) << CBK_->RefLattice();
      else if (!strcmp(response,Commands[indx++]))
	 cout << "NormModulus_= " << NormModulus_ << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "Lambda_= " << Lambda_ << endl;
      else if (!strcmp(response,Commands[indx++]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "BodyForce_[" << i << "]= " << setw(W)
		 << BodyForce_[i] << endl;
	 }
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    cout << "AtomicMass_[" << i << "]= " << setw(W)
		 << AtomicMass_[i] << endl;
	 }
      }
      else if (!strcmp(response,Commands[indx++]))
	 cout << "GridSize_= " << GridSize_ << endl;
      else if (!strcmp(response,Commands[indx++]))
      {
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	    for (int j=i;j<INTERNAL_ATOMS;++j)
	    {
	       cout << "Potential_[" << i << "][" << j << "]= "
		    << setw(W) << Potential_[i][j] << endl;
	    }
      }
      else if (!strcmp(response,Commands[indx++]))
	 cout << "ConvexityDX_= " << ConvexityDX_ << endl;
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
	 cout << "ReferenceBlochWave= " << ReferenceBlochWave(K) << "\t" << K << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << "\tK > ";
	 Vector K(DIM3,0.0);
	 cin >> K;
	 cin.sync(); // clear input
	 cout << "ReferenceDynamicalStiffness= "
	      << setw(W) << ReferenceDynamicalStiffness(K) << endl;
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
	 cout << "energy= " << energy() << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "E0= " << E0() << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "E1= " << E1() << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "E2= " << E2() << endl;
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
	 char datafl[265],prefix[265];
	 cout << "\tdatafile > ";
	 cin >> datafl;
	 cin.sync(); // clear input
	 cout << "\tprefix > ";
	 cin >> prefix;
	 cin.sync(); // clear input
	 cout << "\titer > ";
	 cin >> iter;
	 cin.sync(); // clear input
	 FindLatticeSpacing(datafl,prefix,iter);
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
	 cout << "dbg_ = " << dbg_ << endl;
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
	      << setw(W) << EulerAng_[2] << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << "Rotation_ = "
	      << setw(W) << Rotation_ << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << "Loading_ = "
	      << setw(W) << Loading_ << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << setw(W);
	 PrintCurrentCrystalParamaters(cout);
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << "Entropy = " << setw(W) << Entropy() << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 int n;
	 cout << "\tNoAtoms > ";
	 cin >> n;
	 cout << setw(15) << TranslationProjection1D(n).Transpose() << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 int n,f;
	 cout << "\tFsize > ";
	 cin >> f;
	 cout << "\tNoAtoms > ";
	 cin >> n;
	 cout << setw(15) << TranslationProjection3D(f,n).Transpose() << endl;
      }
      else if (!strcmp(response,"?") ||
	       !strcasecmp(response,"help"))
      {
	 cout << setiosflags(ios::left);
	 for (int i=0;i<NOcommands/2 + NOcommands%2;++i)
	 {
	    cout << "  " << setw(30) << Commands[i];
	    if (i==NOcommands/2 && !NOcommands%2)
	       cout << endl;
	    else
	       cout << setw(30) << Commands[NOcommands/2+i] << endl;
	    
	    if (!((i+1)%30))
	    {
	       cout << "more...." << endl;
	       char ans;
	       cin.sync(); // clear input
	       ans=kbhitWait();
	       if (ans=='q') break;
	    }
	 }
	 cout << resetiosflags(ios::left) << endl;
      }
      else if (!strcmp(response,"\n") ||
	       !strcmp(response,""))
      {
      }
      else
      {
	 cout << "!--- Error - Unknown command ---!" << endl << endl;
      }
      
      cout << endl << prompt;
      cin.getline(response,LINELENGTH);
   }  
}


void MultiLatticeTPP::RefineEqbm(double Tol,int MaxItr,ostream *out)
{
   Vector dx(CBK_->DOFS(),0.0);
   Vector Stress=stress();
   int itr=0;
   
   while ((itr < MaxItr) && Stress.Norm() > Tol)
   {
      ++itr;

#ifdef SOLVE_SVD
      dx = SolveSVD(stiffness(),Stress,MAXCONDITION,Echo_);
#else
      dx = SolvePLU(stiffness(),Stress);
#endif
      SetDOF(CBK_->DOF()-dx);
      
      Stress=stress();

      if (out != NULL)
      {
	 *out << setw(20) << Stress;
	 
	 *out << itr << "\tdx " << dx.Norm() << "\tstress " << Stress.Norm() << endl;
      }
   }
}

void MultiLatticeTPP::PrintCurrentCrystalParamaters(ostream &out)
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

   
   out << "TITLE LatticeStatics crystal structure scaled by 10.0" << endl;
   out << "DIMENSION 3" << endl;
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
       << setw(W) << gamma*180.0/pi << endl;

   out << "SYMMETRY  NUMBER 1  LABEL P1  " << endl;
   out << "SYM MAT  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0 0.0000 0.0000 0.0000" << endl;
   out << endl << "ATOMS" << endl
       << "NAME" << setw(W) << "X" << setw(W) << "Y" << setw(W) << "Z" << endl;
   char const *species[] = {"Ni","Ti","C"};
   out << setw(4) << species[(AtomSpecies_[0] > 3)?3:AtomSpecies_[0]]
       << setw(W) << CBK_->FractionalPosVec(0) << endl;
   for (int i=1;i<INTERNAL_ATOMS;++i)
   {
      out << setw(4) << species[(AtomSpecies_[i] > 3)?3:AtomSpecies_[i]];
      out << setw(W) << CBK_->FractionalPosVec(i) << endl;
   }

   out << "EOF" << endl;

   out << endl
       << "Temperature : " << setw(W) << NTemp_ << endl
       << "Lambda : " << setw(W) << Lambda_ << endl
       << "DOFs : " << setw(W) << CBK_->DOF() << endl;
}
