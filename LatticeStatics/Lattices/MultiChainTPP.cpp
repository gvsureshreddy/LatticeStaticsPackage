#include "MultiChainTPP.h"
#include <cmath>

#include "UtilityFunctions.h"

using namespace std;

MultiChainTPP::~MultiChainTPP()
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
   delete [] AtomPositions_;
}

MultiChainTPP::MultiChainTPP(char *datafile,const char *prefix,int Echo,int Width,int Debug)
{
   Echo_ = Echo;
   dbg_ = Debug;
   // Get Lattice definition
   char tmp[LINELENGTH];
   if(!GetParameter(prefix,"InternalAtoms",datafile,"%u",&INTERNAL_ATOMS)) exit(-1);
   DOFS = INTERNAL_ATOMS;
   if (DOFMAX < DOFS)
   {
      cerr << "Error (MultiChainTPP()): DOFMAX < " << DOFS << " in Lattice.h" << endl;
      exit(-5);
   }
   
   // Set RefLattice_
   RefLattice_.Resize(DIM1,DIM1);
   if(!GetParameter(prefix,"LatticeBasis",datafile,"%lf",&(RefLattice_[0][0]))) exit(-1);
   
   // First Size DOF
   DOF_.Resize(DOFS,0.0);
   DOF_[0] = 1.0;
   // Set AtomPositions_
   AtomPositions_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
   {
      AtomPositions_[i].Resize(DIM1);
      sprintf(tmp,"AtomPosition_%u",i);
      if(!GetParameter(prefix,tmp,datafile,"%lf",&(AtomPositions_[i][0]))) exit(-1);
   }
   
   // Setup Bodyforce_
   BodyForce_ = new Vector[INTERNAL_ATOMS];
   for (int i=0;i<INTERNAL_ATOMS;++i)
      BodyForce_[i].Resize(DIM1,0.0);

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
   if(!GetParameter(prefix,"InfluenceDist",datafile,"%lf",&InfluenceDist_)) exit(-1);
   if(!GetParameter(prefix,"NormModulus",datafile,"%lf",&NormModulus_)) exit(-1);

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
   
   // needed to initialize reference length
   int iter;
   if(!GetParameter(prefix,"MaxIterations",datafile,"%u",&iter)) exit(-1);
   if(!GetParameter(prefix,"BlochWaveGridSize",datafile,"%u",&GridSize_)) exit(-1);

   //set LagrangeCB_
   const char *CBKin[] = {"LagrangeCB","MixedCB"};
   switch (GetStringParameter(prefix,"CBKinematics",datafile,CBKin,2,0))
   {
      case 1:
	 LagrangeCB_ = 0;
	 break;
      case 0:
      default:
	 LagrangeCB_ = 1;
	 break;
   }
   
   // Initiate the Lattice Sum object
   ChainSum_(&DOF_,LagrangeCB_,0,&RefLattice_,INTERNAL_ATOMS,AtomPositions_,Potential_,
	     &InfluenceDist_,&NTemp_);

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
   ChainSum_.Recalc();

   // Initiate the Unit Cell Iterator for Bloch wave calculations.
   ChainIter_(GridSize_);
   if (dbg_)
   {
      if (EnterDebugMode())
      {
	 cout << setw(Width);
	 DebugMode();
      }
   }
   
}

int MultiChainTPP::FindLatticeSpacing(char *datafile,const char *prefix,int iter)
{
   Lambda_=0.0;
   NTemp_=1.0;
   DOF_[0] = 1.0;
   for (int i=1;i<DOFS;i++)
   {
      DOF_[i] = 0.0;
   }

   ChainSum_.Recalc();

   if (Echo_)
      RefineEqbm(1.0e-13,iter,&cout);
   else
      RefineEqbm(1.0e-13,iter,NULL);

   // Clean up numerical round off (at least for zero values)
   for (int i=0;i<DOFS;++i)
   {
      if (fabs(DOF_[i]) < 1.0e-13) DOF_[i] = 0.0;
   }

   // Update RefLattice_
   RefLattice_ = RefLattice_*DOF_[0];

   // update atom pos
   for (int i=1;i<INTERNAL_ATOMS;++i)
   {
      AtomPositions_[i][0] += DOF_[i];
   }

   // reset DOF
   DOF_[0] = 1.0;
   for (int i=1;i<DOFS;i++)
   {
      DOF_[i] = 0.0;
   }

   ChainSum_.Recalc();
   return 0;
}

// Lattice Routines

double MultiChainTPP::PI(double *Dx,double *DX)
{
   return 2.0*Dx[0]*DX[0];
}

double MultiChainTPP::PSI(double *DX)
{
   return 2.0*DX[0]*DX[0];
}

double MultiChainTPP::OMEGA(double *Dx,int p,int q,int i)
{
   return LagrangeCB_ ?
      2.0*DOF_[0]*RefLattice_[0][0]*DELTA(i,p,q)*Dx[0] :
      2.0*DELTA(i,p,q)*Dx[0];
}

double MultiChainTPP::SIGMA(int p,int q,int i,int j)
{
   return LagrangeCB_ ?
      2.0*DOF_[0]*DOF_[0]*RefLattice_[0][0]*DELTA(i,p,q)*RefLattice_[0][0]*DELTA(j,p,q) :
      2.0*DELTA(i,p,q)*DELTA(j,p,q);
}

double MultiChainTPP::GAMMA(double *Dx,double *DX,int p,int q,int i)
{
   return LagrangeCB_ ?
      4.0*RefLattice_[0][0]*DELTA(i,p,q)*Dx[0] :
      2.0*DELTA(i,p,q)*DX[0];
}

double MultiChainTPP::THETA(double *DX,int p,int q,int i)
{
   return LagrangeCB_ ?
      4.0*RefLattice_[0][0]*DELTA(i,p,q)*DX[0] :
      0.0;
}

double MultiChainTPP::XI(int p,int q,int i,int j)
{
   return LagrangeCB_ ?
      4.0*DOF_[0]*RefLattice_[0][0]*DELTA(i,p,q)*RefLattice_[0][0]*DELTA(j,p,q) :
      0.0;
}

double MultiChainTPP::LAMDA(int p,int q,int i,int j)
{
   return LagrangeCB_ ?
      4.0*RefLattice_[0][0]*DELTA(i,p,q)*RefLattice_[0][0]*DELTA(j,p,q) :
      0.0;
}


double MultiChainTPP::energy(PairPotentials::TDeriv dt)
{
   double Phi = 0.0;
   double Vr;

   for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
   {
      // Calculate Phi
      Phi += Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	 NTemp_,ChainSum_.r2(),PairPotentials::Y0,dt);
   }

   // Phi = Phi/(2*Vr*NormModulus)
   Vr = RefLattice_.Det();
   Phi *= 1.0/(2.0*(Vr*NormModulus_));

   // Apply loading potential and Thermal term
   if (dt == PairPotentials::T0)
   {
      // Loading
      Phi -= Lambda_*(DOF_[0] - 1.0);

      // Thermal term
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
      cerr << "Error in MultiChainTPP::energy" << endl;
      exit(-1);
   }
   
   return Phi;
}

Matrix MultiChainTPP::stress(PairPotentials::TDeriv dt,LDeriv dl)
{
   static Matrix S;
   double ForceNorm = 0.0;
   double phi,Vr;
   int i;

   S.Resize(1,DOFS,0.0);

   Vr = RefLattice_.Det();
   
   if (dl==L0)
   {
      for (i=0;i<INTERNAL_ATOMS;++i)
      {
	 BodyForce_[i][0] = 0.0;
      }

      for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
      {
	 // Calculate bodyforce
	 // NOTE: phi1 = d(phi)/d(r2)
	 // We need d(phi)/dr = 2*r*d(phi)/d(r2)
	 phi = 2.0*sqrt(ChainSum_.r2())*ChainSum_.phi1();
	 if (ForceNorm < fabs(-phi/2.0))
	 {
	    ForceNorm = fabs(-phi/2.0);
	 }
	 BodyForce_[ChainSum_.Atom(0)][0] += -phi*ChainSum_.Dx(0)/(2.0*sqrt(ChainSum_.r2()));


	 // Claculate the Stress
	 if (dt == PairPotentials::T0)
	    phi=ChainSum_.phi1();
	 else if (dt == PairPotentials::DT)
	    phi=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	       NTemp_,ChainSum_.r2(),PairPotentials::DY,dt);
	 else
	 {
	    cerr << "Error in MultiChainTPP::stress" << endl;
	    exit(-1);
	 }
	 
	 S[0][0] += phi*PI(ChainSum_.pDx(),ChainSum_.pDX());
	 for (i=1;i<INTERNAL_ATOMS;i++)
	 {
	    S[0][i] += phi*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i);
	 }
      }

      // BodyForce[i] = BodyForce[i] / ForceNorm
      for (i=0;i<INTERNAL_ATOMS;i++)
      {
	 BodyForce_[i][0] /= ForceNorm;
      }

      // S = S/(2*Vr*NormModulus)
      S *= 1.0/(2.0*(Vr*NormModulus_));

      // Load terms
      if (dt == PairPotentials::T0)
      {
	 S[0][0] -= Lambda_;
      }

   }
   else if (dl==DL)
   {
      // dl=DL
      S[0][0] -= 1.0;
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiChainTpp::stress()" << endl;
      exit(-1);
   }
   
   return S;
}
      
Matrix MultiChainTPP::stiffness(PairPotentials::TDeriv dt,LDeriv dl)
{
   static Matrix Phi;
   double phi,phi1;
   int i,j;

   Phi.Resize(DOFS,DOFS,0.0);

   if (dl==L0)
   {
      for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
      {
	 if (dt==PairPotentials::T0)
	 {
	    phi = ChainSum_.phi2();
	    phi1 = ChainSum_.phi1();
	 }
	 else if (dt==PairPotentials::DT)
	 {
	    phi=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	       NTemp_,ChainSum_.r2(),PairPotentials::D2Y,dt);
	    phi1=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	       NTemp_,ChainSum_.r2(),PairPotentials::DY,dt);
	 }
	 else
	 {
	    cerr << "Error in MultiChainTPP::stiffness" << endl;
	    exit(-1);
	 }
	 
	 //Upper Diag Block (1,1)
	 Phi[0][0] += phi*(PI(ChainSum_.pDx(),ChainSum_.pDX())
			   *PI(ChainSum_.pDx(),ChainSum_.pDX()))
	    +phi1*PSI(ChainSum_.pDX());
	 
	 //Lower Diag Block (INTERNAL_ATOMS-1,INTERNAL_ATOMS-1)
	 for (i=1;i<INTERNAL_ATOMS;i++)
	 {
	    for (j=1;j<INTERNAL_ATOMS;j++)
	    {
	       Phi[i][j] +=
		  phi*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j))
		  +phi1*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j);
	    }
	 }
	 
	 //Off Diag Blocks
	 for (i=1;i<INTERNAL_ATOMS;i++)
	 {
	    Phi[0][i] = Phi[i][0] +=
	       phi*(PI(ChainSum_.pDx(),ChainSum_.pDX())
		    *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i))
	       +phi1*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			   ChainSum_.Atom(0),ChainSum_.Atom(1),i);
	 }
      }
      Phi *= 1.0/(2.0*(RefLattice_.Det()*NormModulus_));
   }
   else if (dl==DL)
   {
      // Nothing to do: Phi is zero
   }
   else
   {
      cerr << "Unknown LDeriv dl in MultiChainTpp::stiffness()" << endl;
      exit(-1);
   }   
   return Phi;
}

Matrix MultiChainTPP::E3()
{
   static Matrix Phi;
   double phi,phi1,phi2;
   int i,j,k;

   Phi.Resize(DOFS*DOFS,DOFS,0.0);

   for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
   {
      phi=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	 NTemp_,ChainSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi1=ChainSum_.phi2();
      phi2=ChainSum_.phi1();
	 
      // DF^3 block
      Phi[0][0] +=
	 phi*(PI(ChainSum_.pDx(),ChainSum_.pDX())
	      *PI(ChainSum_.pDx(),ChainSum_.pDX())
	      *PI(ChainSum_.pDx(),ChainSum_.pDX()))
	 +phi1*(3.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
		*PSI(ChainSum_.pDX()));
      // DV^3 block
      for (i=1;i<INTERNAL_ATOMS;i++)
	 for (j=1;j<INTERNAL_ATOMS;j++)
	    for (k=1;k<INTERNAL_ATOMS;k++)
	    {
	       Phi[i*DOFS + j][k] +=
		  phi*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k))
		  +phi1*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			 *SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			 + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			 *SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,k)
			 + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			 *SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j));
	    }
      // DU^2DV blocks
      for (i=1;i<INTERNAL_ATOMS;i++)
      {
	 Phi[0][i] = Phi[i*DOFS][0] = Phi[i][0] += (
	    phi*(PI(ChainSum_.pDx(),ChainSum_.pDX())
		 *PI(ChainSum_.pDx(),ChainSum_.pDX())
		 *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i))
	    +phi1*(PI(ChainSum_.pDx(),ChainSum_.pDX())
		   *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		   + PI(ChainSum_.pDx(),ChainSum_.pDX())
		   *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		   + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		   *PSI(ChainSum_.pDX()))
	    +phi2*THETA(ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i));
      }
      // DV^2DU blocks
      for (i=1;i<INTERNAL_ATOMS;i++)
	 for (j=1;j<INTERNAL_ATOMS;j++)
	 {
	    Phi[i*DOFS + j][0] = Phi[i*DOFS][j] =
	       Phi[i][j] += (
		  phi*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		       *PI(ChainSum_.pDx(),ChainSum_.pDX()))
		  +phi1*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			 *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
				ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			 + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			 *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
				ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			 + PI(ChainSum_.pDx(),ChainSum_.pDX())
			 *SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j))
		  +phi2*XI(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j));
	 }
   }

   
   // Phi = Phi/(2*Vr*NormModulus)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*NormModulus_));

   return Phi;
}

Matrix MultiChainTPP::E4()
{
   static Matrix Phi;
   double phi,phi1,phi2,phi3;
   int i,j,k,m;

   Phi.Resize(DOFS*DOFS,DOFS*DOFS,0.0);

   for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
   {  
      phi=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	 NTemp_,ChainSum_.r2(),PairPotentials::D4Y,PairPotentials::T0);
      phi1=Potential_[ChainSum_.Atom(0)][ChainSum_.Atom(1)]->PairPotential(
	 NTemp_,ChainSum_.r2(),PairPotentials::D3Y,PairPotentials::T0);
      phi2=ChainSum_.phi2();
      phi3=ChainSum_.phi1();
      
      // DU^4 block 
      Phi[0][0]+=
	 phi*(PI(ChainSum_.pDx(),ChainSum_.pDX())
	      *PI(ChainSum_.pDx(),ChainSum_.pDX())
	      *PI(ChainSum_.pDx(),ChainSum_.pDX())
	      *PI(ChainSum_.pDx(),ChainSum_.pDX()))
	 +phi1*(
	    6.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
	    *PI(ChainSum_.pDx(),ChainSum_.pDX())
	    *PSI(ChainSum_.pDX()))
	 +phi2*(
	    3.0*PSI(ChainSum_.pDX())
	    *PSI(ChainSum_.pDX()));      
      // DV^4 block
      for (i=1;i<INTERNAL_ATOMS;i++)
	 for (j=1;j<INTERNAL_ATOMS;j++)
	    for (k=1;k<INTERNAL_ATOMS;k++)
	       for (m=1;m<INTERNAL_ATOMS;m++)
	       {
		  Phi[i*DOFS+j][k*DOFS+m] +=
		     phi*(OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			  *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			  *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			  *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),m))
		     +phi1*(
			OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),k,m)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,m)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),m)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,m)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),m)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),m)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j))
		     +phi2*(
			SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),k,m)
			+ SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,m)
			+ SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,m)
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k));
	       }

      // DU^3DV blocks
      for (i=1;i<INTERNAL_ATOMS;i++)
      {
	 Phi[0][i] = Phi[0][i*DOFS] = Phi[i][0] = Phi[i*DOFS][0] += (
	    phi*(
	       PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i))
	    +phi1*(
	       3.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
	       +3.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
	       *PSI(ChainSum_.pDX()))
	    +phi2*(
	       3.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
	       *THETA(ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
	       +3.0*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
	       *PSI(ChainSum_.pDX())));
      }
      // DV^3DU blocks
      for (i=1;i<INTERNAL_ATOMS;i++)
	 for (j=1;j<INTERNAL_ATOMS;j++)
	    for (k=1;k<INTERNAL_ATOMS;k++)
	    {
	       Phi[i*DOFS+j][k*DOFS] = Phi[i*DOFS+j][k] = Phi[i*DOFS][j*DOFS+k] =
		  Phi[i][j*DOFS+k] += (
		     phi*(
			OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*PI(ChainSum_.pDx(),ChainSum_.pDX()))
		     +phi1*(
			OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*PI(ChainSum_.pDx(),ChainSum_.pDX())
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*PI(ChainSum_.pDx(),ChainSum_.pDX())
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*PI(ChainSum_.pDx(),ChainSum_.pDX())
			*SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j))
		     +phi2*(
			OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
			*XI(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			*XI(ChainSum_.Atom(0),ChainSum_.Atom(1),i,k)
			+ OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			*XI(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
			+ SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),k)
			+ SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,k)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),j)
			+ SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),j,k)
			*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			       ChainSum_.Atom(0),ChainSum_.Atom(1),i)));
	    }
      // DU^2DV^2 blocks
      for (i=1;i<INTERNAL_ATOMS;i++)
	 for (j=1;j<INTERNAL_ATOMS;j++)
	 {
	    Phi[0][i*DOFS+j] =
	       Phi[j*DOFS][i] =
	       Phi[i*DOFS+j][0] =
	       Phi[i][j*DOFS] += (
		  phi*(
		     PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j))
		  +phi1*(
		     PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
		     + 2.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			    ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     + 2.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			    ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     *OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     *PSI(ChainSum_.pDX()))
		  +phi2*(
		     2.0*PI(ChainSum_.pDx(),ChainSum_.pDX())
		     *XI(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
		     + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     *THETA(ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     + OMEGA(ChainSum_.pDx(),ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     *THETA(ChainSum_.pDX(),ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     + 2.0*GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
				 ChainSum_.Atom(0),ChainSum_.Atom(1),i)
		     *GAMMA(ChainSum_.pDx(),ChainSum_.pDX(),
			    ChainSum_.Atom(0),ChainSum_.Atom(1),j)
		     + SIGMA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j)
		     *PSI(ChainSum_.pDX()))
		  +phi3*LAMDA(ChainSum_.Atom(0),ChainSum_.Atom(1),i,j));
	 }
   }
   
   
   // Phi = Phi/(2*Vr*NormModulus)
   Phi *= 1.0/(2.0*(RefLattice_.Det()*NormModulus_));
   
   return Phi;
}

Matrix MultiChainTPP::CondensedModuli()
{
   Matrix stiff = stiffness();
   int intrn = DOFS-1;
   Matrix CM(1,1), IM(intrn,intrn);
   
   CM[0][0] = stiff[0][0];

   // Make sure there are internal DOF's
   if (intrn)
   {
      for (int i=0;i<intrn;i++)
	 for (int j=0;j<intrn;j++)
	 {
	    IM[i][j] = stiff[1+i][1+j];
	 }
      IM = IM.Inverse();
      
      // Set up Condensed Moduli
      for (int m=0;m<intrn;m++)
	 for (int n=0;n<intrn;n++)
	 {
	    CM[0][0] -= stiff[0][1+m]*IM[m][n]*stiff[1+n][0];
	 }
   }
   
   return CM;
}

int MultiChainTPP::comp(const void *a,const void *b)
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

int MultiChainTPP::abscomp(const void *a,const void *b)
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

void MultiChainTPP::interpolate(Matrix *EigVals,int zero,int one,int two)
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

CMatrix MultiChainTPP::ReferenceDynamicalStiffness(Vector &K)
{
   static CMatrix Dk;
   static double pi = 4.0*atan(1.0);
   static MyComplexDouble Ic(0,1);
   static MyComplexDouble A = 2.0*pi*Ic;

   Dk.Resize(INTERNAL_ATOMS,INTERNAL_ATOMS,0.0);
   
   for (ChainSum_.Reset();!ChainSum_.Done();++ChainSum_)
   {
      // Calculate Dk
      if (ChainSum_.Atom(0) != ChainSum_.Atom(1))
      {
	 // y != y' terms (i.e., off diagonal terms)
	 Dk[ChainSum_.Atom(0)][ChainSum_.Atom(1)] +=
	    (-2.0*ChainSum_.phi1()
	     -4.0*ChainSum_.Dx(0)*ChainSum_.Dx(0)*ChainSum_.phi2())
	    *exp(A*K[0]*ChainSum_.DX(0));

	 // y==y' components (i.e., Phi(0,y,y) term)
	 Dk[ChainSum_.Atom(0)][ChainSum_.Atom(0)] +=
	    (2.0*ChainSum_.phi1()+4.0*ChainSum_.Dx(0)*ChainSum_.Dx(0)*ChainSum_.phi2());
      }
      else
      {
	 Dk[ChainSum_.Atom(0)][ChainSum_.Atom(1)] +=
	    (-2.0*ChainSum_.phi1()
	     -4.0*ChainSum_.Dx(0)*ChainSum_.Dx(0)*ChainSum_.phi2())
	    *(exp(A*K[0]*ChainSum_.DX(0)) - 1.0);
      }
   }
   // Normalize through the Mass Matrix
   for (int p=0;p<INTERNAL_ATOMS;++p)
      for (int q=0;q<INTERNAL_ATOMS;++q)
      {
	 Dk[p][q] /= sqrt(AtomicMass_[p]*AtomicMass_[q]);
      }
   
   return Dk;
}

void MultiChainTPP::ReferenceDispersionCurves(Vector K,int NoPTS,const char *prefix,
						ostream &out)
{
   int w=out.width();
   out.width(0);
   if (Echo_) cout.width(0);

   double InverseLat;
   InverseLat = 1.0/RefLattice_[0][0];

   Matrix EigVal[3];
   for (int i=0;i<3;++i) EigVal[i].Resize(1,INTERNAL_ATOMS);

   double Z1,Z2;
   Z1 = K[0];
   Z2 = K[3];

   Z1 = InverseLat*Z1;
   Z2 = InverseLat*Z2;
   
   Vector Z(1);
   double
      DZ=Z2-Z1;
   double dz = 1.0/(NoPTS-1);
   for (int k=0;k<2;++k)
   {
      Z[0] = Z1 + (k*dz)*DZ;
      EigVal[k] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[k][0],INTERNAL_ATOMS,sizeof(double),&comp);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS;++i)
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
      Z[0] = Z1 + (k*dz)*DZ;
      EigVal[two] = HermiteEigVal(ReferenceDynamicalStiffness(Z));
      qsort(EigVal[two][0],INTERNAL_ATOMS,sizeof(double),&comp);
      interpolate(EigVal,zero,one,two);
      
      out << prefix << setw(w) << k*dz;
      if (Echo_) cout << prefix << setw(w) << k*dz;
      for (int i=0;i<INTERNAL_ATOMS;++i)
      {
	 out << setw(w) << EigVal[two][0][i];;
	 if (Echo_) cout << setw(w) << EigVal[two][0][i];;
      }
      out << endl;
      if (Echo_) cout << endl;

      zero = (++zero)%3; one = (zero+1)%3; two = (one+1)%3;
   }
}

int MultiChainTPP::ReferenceBlochWave(Vector &K)
{
   static CMatrix A(INTERNAL_ATOMS,INTERNAL_ATOMS);
   static Matrix EigVals(1,INTERNAL_ATOMS);
   static double InverseLat;
   static Vector Z(1);

   for (int i=0;i<K.Dim();++i) K[i]=0.0;

   InverseLat = 1.0/RefLattice_[0][0];

   // Iterate over points in linear unit cell
   for (ChainIter_.Reset();!ChainIter_.Done();++ChainIter_)
   {
      K[0] = ChainIter_[0];
      
      Z = InverseLat*K;
      A = ReferenceDynamicalStiffness(Z);
      
      EigVals = HermiteEigVal(A);
      
      for (int i=0;i<INTERNAL_ATOMS;++i)
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

void MultiChainTPP::LongWavelengthModuli(double dk, int gridsize,const char *prefix,
					   ostream &out)
{
}

void MultiChainTPP::NeighborDistances(int cutoff,ostream &out)
{
   Matrix NeighborDist =
      ChainSum_.NeighborDistances(cutoff,pow(double(10),double(-(out.precision()-1))));
   
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

void MultiChainTPP::Print(ostream &out,PrintDetail flag)
{
   static int W;
   static int NoNegEigVal;
   static double MinEigVal;
   static double engy,entropy,heatcapacity;
   static Matrix
      str(1,DOFS),
      stiff(DOFS,DOFS),
      EigenValues(1,DOFS),
      CondEV(1,1);
   static Matrix
      CondModuli(1,1);
   static int RankOneConvex;
   static Vector K(1);
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
   for (int i=0;i<DOFS;i++)
   {
      if (EigenValues[0][i] < 0)
	 NoNegEigVal++;

      if (MinEigVal > EigenValues[0][i])
	 MinEigVal = EigenValues[0][i];
   }

   CondModuli = CondensedModuli();

   CondEV=CondModuli;
   RankOneConvex = (CondEV[0][0] > 0) ? 1 : 0;

   K.Resize(1,0.0);
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
	 out << "MultiChainTPP:" << endl << endl;
	 out << "LagrangeCB: " << LagrangeCB_ << endl;
	 out << "RefLattice_ : " << setw(W) << RefLattice_;
	 for (int i=0;i<INTERNAL_ATOMS;++i)
	 {
	    out << "Atom_" << i << "          "
		<< "Species : " << setw(5) << AtomSpecies_[i]
		<< "          Position : " << setw(W) << AtomPositions_[i] << endl;
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
	 // also send to cout
	 if (Echo_)
	 {
	    cout << "MultiChainTPP:" << endl << endl;
	    cout << "LagrangeCB: " << LagrangeCB_ << endl;
	    cout << "RefLattice_ : " << setw(W) << RefLattice_;
	    for (int i=0;i<INTERNAL_ATOMS;++i)
	    {
	       cout << "Atom_" << i << "          "
		<< "Species : " << setw(5) << AtomSpecies_[i]
		<< "          Position : " << setw(W) << AtomPositions_[i] << endl;
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
	 }
	 // passthrough to short
      case PrintShort:
	 out << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
	     << "Lambda (Normalized): " << setw(W) << Lambda_ << endl
	     << "DOF's :" << endl << setw(W) << DOF_ << endl
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
	     << "BlochWave Stability:" << setw(W) << BlochWaveStable << ", "
	     << setw(W) << K << endl;
	 // send to cout also
	 if (Echo_)
	 {
	    cout << "Temperature (Ref Normalized): " << setw(W) << NTemp_ << endl
		 << "Lambda (Normalized): " << setw(W) << Lambda_ << endl
		 << "DOF's :" << endl << setw(W) << DOF_ << endl
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

ostream &operator<<(ostream &out,MultiChainTPP &A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}


//---------------------- Debug Mode Handler --------------------------


void MultiChainTPP::DebugMode()
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
      "Entropy"
   };
   int NOcommands=41;
   
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
	 cout << "DOFS = " << DOFS << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "InfluenceDist_ = " << InfluenceDist_ << endl;
      else if (!strcmp(response,Commands[indx++]))
	 cout << "NTemp_ = " << NTemp_ << endl;
      else if (!strcmp(response,Commands[indx++]))
      {
	 for (int i=0;i<DOFS;++i)
	    cout << "DOF_[" << i << "] = " << DOF_[i] << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
	 cout << "RefLattice_= " << setw(W) << RefLattice_;
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
	 cout << "stress= " << setw(W) << stress();
      else if (!strcmp(response,Commands[indx++]))
	 cout << "stiffness= " << setw(W) << stiffness();
      else if (!strcmp(response,Commands[indx++]))
	 cout << "CondensedModuli= " << setw(W) << CondensedModuli();
      else if (!strcmp(response,Commands[indx++]))
      {
	 Vector K(6,0.0);
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
	 Vector K(1,0.0);
	 cout << "ReferenceBlochWave= " << ReferenceBlochWave(K) << "\t" << K << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 cout << "\tK > ";
	 Vector K(1,0.0);
	 cin >> K;
	 cin.sync(); // clear input
	 cout << "ReferenceDynamicalStiffness= "
	      << setw(W) << ReferenceDynamicalStiffness(K) << endl;
      }
      else if (!strcmp(response,Commands[indx++]))
      {
	 Vector DOF(DOFS,0.0);
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
	 cout << "E0= " << setw(W) << E0();
      else if (!strcmp(response,Commands[indx++]))
	 cout << "E1= " << setw(W) << E1();
      else if (!strcmp(response,Commands[indx++]))
	 cout << "E2= " << setw(W) << E2();
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
	 cout << "Entropy = " << setw(W) << Entropy() << endl;
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


void MultiChainTPP::RefineEqbm(double Tol,int MaxItr,ostream *out)
{
   Vector dx(DOFS,0.0);
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

      SetDOF(DOF_-dx);

      Stress=stress();

      if (out != NULL)
      {
	 *out << setw(20) << Stress;
	 
	 *out << itr << "\tdx " << dx.Norm() << "\tstress " << Stress.Norm() << endl;
      }
   }
}
