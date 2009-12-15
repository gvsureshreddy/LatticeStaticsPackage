#include "DFTExternal.h"

using namespace std;

fstream DFTExternal::dbug;

const double DFTExternal::Alt[3][3][3] = {{{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 1.0},
                                           {0.0, -1.0, 0.0}},
                                          {{0.0, 0.0, -1.0},
                                           {0.0, 0.0, 0.0},
                                           {1.0, 0.0, 0.0}},
                                          {{0.0, 1.0, 0.0},
                                           {-1.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}}};

const double DFTExternal::Del[3][3] = {{1.0, 0.0, 0.0},
                                       {0.0, 1.0, 0.0},
                                       {0.0, 0.0, 1.0}};

DFTExternal::~DFTExternal()
{
   cout << "DFTExternal Function Calls:\n"
        << "\tE0 calls - " << CallCount_[1] << "\n"
        << "\tE1 calls - " << CallCount_[2] << "\n"
        << "\tE1DLoad calls - " << CallCount_[3] << "\n"
        << "\tE2 calls - " << CallCount_[4] << "\n";
   dbug.close();
}

DFTExternal::DFTExternal(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input,Echo),
   Lambda_(0.0),
   Width_(Width)
{
   dbug.open("DFTExternal-dof-force.data", ios::out);
   dbug << fixed << setprecision(15);

   
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash,"DFTExternal");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   DOF_.Resize(DOFS_,0.0);
   E1CachedValue_.Resize(DOFS_);
   E1DLoadCachedValue_.Resize(DOFS_);
   E2CachedValue_.Resize(DOFS_,DOFS_);
   EmptyV_.Resize(2,0.0);
   EmptyM_.Resize(2,2,0.0);
   
   LoadParameter_ = Load;
   for (int i=0;i<cachesize;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
}

void DFTExternal::UpdateValues(UpdateFlag flag) const
{
   int retid;
   fstream in;
   fstream out;
   out.open("DFTModelInput",ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << "DFTModelInput " << " for write"
           << "\n";
      exit(-1);
   }
   out << setiosflags(ios::scientific) << setprecision(20);
   for (int i=0;i<6;++i)
   {
      out << setw(30) << DOF_[i];
   }
   out << "\n";
   int q=0;
   for (int i=6;i<DOFS_;++i)
   {
      out << setw(30) << DOF_[i];
      q = (q+1)%3;
      if (q==0) out << "\n";
   }

   out.close();

   // update with correct command.
   if (flag == NoStiffness)
   {
      retid=system("./script_main 0 >& /dev/null");
   }
   else
   {
      retid=system("./script_main 2 >& /dev/null");
   }
   cerr << "DFTExternal (flag=" << flag << ") system() call returned with id: " 
        << retid << endl;

   // calculate pressure terms.
   Matrix B(3,3);
   B[0][0] = DOF_[0];
   B[1][1] = DOF_[1];
   B[2][2] = DOF_[2];
   B[1][2]=B[2][1] = DOF_[3];
   B[2][0]=B[0][2] = DOF_[4];
   B[0][1]=B[1][0] = DOF_[5];
   Matrix U(3,3);
   U[0][0] = 1.0+DOF_[0];
   U[1][1] = 1.0+DOF_[1];
   U[2][2] = 1.0+DOF_[2];
   U[1][2]=U[2][1] = DOF_[3];
   U[2][0]=U[0][2] = DOF_[4];
   U[0][1]=U[1][0] = DOF_[5];
   double UDet = U.Det();
   Matrix Uinv = U.Inverse();
   double PressureEnergy = Lambda_*UDet;
   Matrix PressureTerm = PressureEnergy*(U.Inverse()).Transpose();
   Vector PressureStress(DOFS_,0.0);
   PressureStress[0] = PressureTerm[0][0];
   PressureStress[1] = PressureTerm[1][1];
   PressureStress[2] = PressureTerm[2][2];
   PressureStress[3] = PressureTerm[1][2]+PressureTerm[2][1];
   PressureStress[4] = PressureTerm[2][0]+PressureTerm[0][2];
   PressureStress[5] = PressureTerm[0][1]+PressureTerm[1][0];
   Matrix PressureStiffness(DOFS_,DOFS_,0.0);
   if (flag==NeedStiffness)
   {
      for (int q=0;q<3;++q)
      {
         for (int s=0;s<3;++s)
         {
            PressureStiffness[0][0] += Lambda_*0.25*
               (Alt[0][0][q]*Alt[0][0][s] +
                Alt[0][0][q]*Alt[0][0][s] +
                Alt[0][0][q]*Alt[0][0][s] +
                Alt[0][0][q]*Alt[0][0][s])*U[q][s];
            PressureStiffness[0][1]=PressureStiffness[1][0] += Lambda_*0.25*
               (Alt[0][1][q]*Alt[0][1][s] +
                Alt[0][1][q]*Alt[0][1][s] +
                Alt[0][1][q]*Alt[0][1][s] +
                Alt[0][1][q]*Alt[0][1][s])*U[q][s];
            PressureStiffness[0][2]=PressureStiffness[2][0] += Lambda_*0.25*
               (Alt[0][2][q]*Alt[0][2][s] +
                Alt[0][2][q]*Alt[0][2][s] +
                Alt[0][2][q]*Alt[0][2][s] +
                Alt[0][2][q]*Alt[0][2][s])*U[q][s];
            PressureStiffness[0][3]=PressureStiffness[3][0] += Lambda_*0.25*
               (Alt[0][1][q]*Alt[0][2][s] +
                Alt[0][1][q]*Alt[0][2][s] +
                Alt[0][2][q]*Alt[0][1][s] +
                Alt[0][2][q]*Alt[0][1][s])*U[q][s];
            PressureStiffness[0][4]=PressureStiffness[4][0] += Lambda_*0.25*
               (Alt[0][2][q]*Alt[0][0][s] +
                Alt[0][2][q]*Alt[0][0][s] +
                Alt[0][0][q]*Alt[0][2][s] +
                Alt[0][0][q]*Alt[0][2][s])*U[q][s];
            PressureStiffness[0][5]=PressureStiffness[5][0] += Lambda_*0.25*
               (Alt[0][0][q]*Alt[0][1][s] +
                Alt[0][0][q]*Alt[0][1][s] +
                Alt[0][1][q]*Alt[0][0][s] +
                Alt[0][1][q]*Alt[0][0][s])*U[q][s];

            PressureStiffness[1][1] += Lambda_*0.25*
               (Alt[1][1][q]*Alt[1][1][s] +
                Alt[1][1][q]*Alt[1][1][s] +
                Alt[1][1][q]*Alt[1][1][s] +
                Alt[1][1][q]*Alt[1][1][s])*U[q][s];
            PressureStiffness[1][2]=PressureStiffness[2][1] += Lambda_*0.25*
               (Alt[1][2][q]*Alt[1][2][s] +
                Alt[1][2][q]*Alt[1][2][s] +
                Alt[1][2][q]*Alt[1][2][s] +
                Alt[1][2][q]*Alt[1][2][s])*U[q][s];
            PressureStiffness[1][3]=PressureStiffness[3][1] += Lambda_*0.25*
               (Alt[1][1][q]*Alt[1][2][s] +
                Alt[1][1][q]*Alt[1][2][s] +
                Alt[1][2][q]*Alt[1][1][s] +
                Alt[1][2][q]*Alt[1][1][s])*U[q][s];
            PressureStiffness[1][4]=PressureStiffness[4][1] += Lambda_*0.25*
               (Alt[1][2][q]*Alt[1][0][s] +
                Alt[1][2][q]*Alt[1][0][s] +
                Alt[1][0][q]*Alt[1][2][s] +
                Alt[1][0][q]*Alt[1][2][s])*U[q][s];
            PressureStiffness[1][5]=PressureStiffness[5][1] += Lambda_*0.25*
               (Alt[1][0][q]*Alt[1][1][s] +
                Alt[1][0][q]*Alt[1][1][s] +
                Alt[1][1][q]*Alt[1][0][s] +
                Alt[1][1][q]*Alt[1][0][s])*U[q][s];
            
            PressureStiffness[2][2] += Lambda_*0.25*
               (Alt[2][2][q]*Alt[2][2][s] +
                Alt[2][2][q]*Alt[2][2][s] +
                Alt[2][2][q]*Alt[2][2][s] +
                Alt[2][2][q]*Alt[2][2][s])*U[q][s];
            PressureStiffness[2][3]=PressureStiffness[3][2] += Lambda_*0.25*
               (Alt[2][1][q]*Alt[2][2][s] +
                Alt[2][1][q]*Alt[2][2][s] +
                Alt[2][2][q]*Alt[2][1][s] +
                Alt[2][2][q]*Alt[2][1][s])*U[q][s];
            PressureStiffness[2][4]=PressureStiffness[4][2] += Lambda_*0.25*
               (Alt[2][2][q]*Alt[2][0][s] +
                Alt[2][2][q]*Alt[2][0][s] +
                Alt[2][0][q]*Alt[2][2][s] +
                Alt[2][0][q]*Alt[2][2][s])*U[q][s];
            PressureStiffness[2][5]=PressureStiffness[5][2] += Lambda_*0.25*
               (Alt[2][0][q]*Alt[2][1][s] +
                Alt[2][0][q]*Alt[2][1][s] +
                Alt[2][1][q]*Alt[2][0][s] +
                Alt[2][1][q]*Alt[2][0][s])*U[q][s];

            PressureStiffness[3][3] += Lambda_*0.25*
               (Alt[1][1][q]*Alt[2][2][s] +
                Alt[2][1][q]*Alt[1][2][s] +
                Alt[1][2][q]*Alt[2][1][s] +
                Alt[2][2][q]*Alt[1][1][s])*U[q][s];
            PressureStiffness[3][4]=PressureStiffness[4][3] += Lambda_*0.25*
               (Alt[1][2][q]*Alt[2][0][s] +
                Alt[2][2][q]*Alt[1][0][s] +
                Alt[1][0][q]*Alt[2][2][s] +
                Alt[2][0][q]*Alt[1][2][s])*U[q][s];
            PressureStiffness[3][5]=PressureStiffness[5][3] += Lambda_*0.25*
               (Alt[1][0][q]*Alt[2][1][s] +
                Alt[2][0][q]*Alt[1][1][s] +
                Alt[1][1][q]*Alt[2][0][s] +
                Alt[2][1][q]*Alt[1][0][s])*U[q][s];
            
            PressureStiffness[4][4] += Lambda_*0.25*
               (Alt[2][2][q]*Alt[0][0][s] +
                Alt[0][2][q]*Alt[2][0][s] +
                Alt[2][0][q]*Alt[0][2][s] +
                Alt[0][0][q]*Alt[2][2][s])*U[q][s];
            PressureStiffness[4][5]=PressureStiffness[5][4] += Lambda_*0.25*
               (Alt[2][0][q]*Alt[0][1][s] +
                Alt[0][0][q]*Alt[2][1][s] +
                Alt[2][1][q]*Alt[0][0][s] +
                Alt[0][1][q]*Alt[2][0][s])*U[q][s];

            PressureStiffness[5][5] += Lambda_*0.25*
               (Alt[0][0][q]*Alt[1][1][s] +
                Alt[1][0][q]*Alt[0][1][s] +
                Alt[0][1][q]*Alt[1][0][s] +
                Alt[1][1][q]*Alt[0][0][s])*U[q][s];
         }
      }
   }
   
   in.open("in.latpak",ios::in);
   if (in.fail())
   {
      cerr << "Error: Unable to open file : " << "in.latpak " << "for read" << "\n";
      exit(-2);
   }


   // Read in DFT information
   Matrix Stresses(3,3);
   Matrix Forces(DOFS_-6,1);
   
   in >> DFTEnergyCachedValue_;
   in >> Stresses;
   in >> Forces;

   // set dwdc
   Matrix dwdc(3,3);

   for (int s=0;s<3;++s)
      for (int t=0;t<3;++t)
      {
         dwdc[s][t] = 0.0;
         for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
               dwdc[s][t] += 0.5*UDet*Uinv[s][i]*Stresses[i][j]*Uinv[t][j];
      }

   // set BFB energy value
   E0CachedValue_ = UDet*DFTEnergyCachedValue_;
   // add phantom energy for translations
   double TrEig[3] = {1.0, 2.0, 3.0};
   Vector Tsq(3);
   int InternalAtoms = (DOFS_-6)/3;
   for (int j=0;j<3;++j)
   {
      Tsq[j]=0.0;
      for (int i=0;i<InternalAtoms;++i)
      {
         Tsq[j]+=DOF_[6+3*i+j];
      }
      Tsq[j] = (Tsq[j]*Tsq[j])/InternalAtoms;
   }
   E0CachedValue_ += -PressureEnergy + 0.5*(TrEig[0]*Tsq[0] + TrEig[1]*Tsq[1] + TrEig[2]*Tsq[2]);
   Cached_[1] = 1;

   // set BFB stress values
   E1CachedValue_.Resize(DOFS_,0.0);
   double bfbstresses[3][3];
   for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
      {
         bfbstresses[i][j] =0.0;
         for (int r=0;r<3;++r)
            bfbstresses[i][j] += U[i][r]*dwdc[r][j] + U[j][r]*dwdc[r][i];
      }
   for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
         E1CachedValue_[((i==j)? i : 6-(i+j))] += bfbstresses[i][j];

   // set BFB force values
   for (int i=0;i<DOFS_-6;++i)
   {
      E1CachedValue_[6+i] = UDet*Forces[i][0];
   }

   // Phantom energy terms for the gradient.
   Vector T(3);
   Vector ME1(DOFS_,0.0);
   for (int j=0;j<3;++j)
   {
      T[j]=0.0;
      for (int i=0;i<InternalAtoms;++i)
         T[j]+=DOF_[6+3*i+j];
      T[j]/=InternalAtoms;
   }
   for (int i=0;i<InternalAtoms;++i)
   {
      for (int j=0;j<3;++j)
      {
         ME1[6+3*i+j] += TrEig[j]*T[j];
      }
   }
   E1CachedValue_ += -PressureStress + ME1;
   Cached_[2] = 1;
   

   // Stiffness terms if needed
   if (flag==NeedStiffness)
   {
      // read in DFT information
      Matrix StrainStiffnesses(6,6);
      Matrix CrossStiffnesses(6,DOFS_-6);
      Matrix ShiftStiffnesses(DOFS_-6,DOFS_-6);
      
      in >> StrainStiffnesses;
      in >> CrossStiffnesses;
      in >> ShiftStiffnesses;
      
      // set d2wdcds and d2wd2c
      Matrix d2wdcds[3]; for (int i=0;i<3;++i) d2wdcds[i].Resize(3,DOFS_-6);
      double strainstiffnesses[3][3][3][3];
      double d2wd2c[3][3][3][3];
      
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<DOFS_-6;++k)
            {
               d2wdcds[i][j][k] = 0.0;
               for (int a=0;a<3;++a)
                  for (int b=0;b<3;++b)
                  {
                     d2wdcds[i][j][k] += 0.5*UDet*Uinv[i][a]*
                        CrossStiffnesses[((a==b)? a : 6-(a+b))][k]*Uinv[j][b];
                  }
            }

      Matrix bfbcrossstiffness[3]; for (int i=0;i<3;++i) bfbcrossstiffness[i].Resize(3,DOFS_-6);
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<DOFS_-6;++k)
            {
               bfbcrossstiffness[i][j][k] = 0.0;
               for (int r=0;r<3;++r)
                  bfbcrossstiffness[i][j][k]
                     += U[i][r]*d2wdcds[r][j][k] + U[j][r]*d2wdcds[r][i][k];
            }
      
      for (int a=0;a<3;++a)
         for (int b=0;b<3;++b)
            for (int d=0;d<3;++d)
               for (int g=0;g<3;++g)
               {
                  strainstiffnesses[a][b][d][g]
                     = StrainStiffnesses[((a==b)? a : 6-(a+b))][((d==g)? d : 6-(d+g))];
               }
      
      for (int a=0;a<3;++a)
         for (int b=0;b<3;++b)
            for (int d=0;d<3;++d)
               for (int g=0;g<3;++g)
               {
                  d2wd2c[a][b][d][g] = 0.0;
                  for (int i=0;i<3;++i)
                     for (int j=0;j<3;++j)
                        for (int k=0;k<3;++k)
                           for (int l=0;l<3;++l)
                              d2wd2c[a][b][d][g] += 0.25*UDet*(
                                 Uinv[a][i]*strainstiffnesses[i][j][k][l]
                                 *Uinv[b][j]*Uinv[d][k]*Uinv[g][l]);
                  
                  for (int j=0;j<3;++j)
                     for (int k=0;k<3;++k)
                        for (int l=0;l<3;++l)
                           d2wd2c[a][b][d][g] += 0.25*UDet*(
                              - 0.5*(Uinv[a][k]*Uinv[d][k]*Stresses[l][j]
                                     *Uinv[g][l]*Uinv[b][j]
                                     + Uinv[a][l]*Uinv[g][l]*Stresses[k][j]
                                     *Uinv[d][k]*Uinv[b][j]));
               }


      double bfbstrainstiffness[3][3][3][3];
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
               for (int l=0;l<3;++l)
               {
                  bfbstrainstiffness[i][j][k][l] = 0.5*(Del[i][k]*dwdc[l][j]
                                                        + Del[i][l]*dwdc[k][j]
                                                        + Del[j][k]*dwdc[l][i]
                                                        + Del[j][l]*dwdc[k][i]);
                  for (int r=0;r<3;++r)
                     for (int a=0;a<3;++a)
                     {
                        bfbstrainstiffness[i][j][k][l] +=
                           U[i][r]*d2wd2c[r][j][a][l]*U[k][a]
                           + U[i][r]*d2wd2c[r][j][a][k]*U[l][a]
                           + U[j][r]*d2wd2c[r][i][a][l]*U[k][a]
                           + U[j][r]*d2wd2c[r][i][a][k]*U[l][a];
                     }
               }

      // Set bfb stiffness value
      E2CachedValue_.Resize(DOFS_,DOFS_,0.0);
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
               for (int l=0;l<3;++l)
               {
                  E2CachedValue_[((i==j)? i : 6-(i+j))][((k==l)? k : 6-(k+l))] +=
                     bfbstrainstiffness[i][j][k][l];
               }
      
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<DOFS_-6;++k)
            {
               E2CachedValue_[((i==j)? i : 6-(i+j))][6+k] += bfbcrossstiffness[i][j][k];
               E2CachedValue_[6+k][((i==j)? i : 6-(i+j))] += bfbcrossstiffness[i][j][k];
            }

      for (int i=0;i<DOFS_-6;++i)
         for (int j=0;j<DOFS_-6;++j)
         {
            E2CachedValue_[6+i][6+j] = UDet*ShiftStiffnesses[i][j];
         }

      // Phantom Energy terms
      Matrix ME2(DOFS_,DOFS_,0.0);
      for (int i=0;i<InternalAtoms;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<InternalAtoms;++k)
            {
               ME2[6+3*i+j][6+3*k+j] += TrEig[j]/InternalAtoms;
            }

      E2CachedValue_ += -PressureStiffness + ME2;
      Cached_[4] = 1;
   }
   in.close();

   //output force/dof data
   dbug << setw(Width_) << DOF_ << setw(Width_) << Lambda_  << setw(Width_) << flag << endl;
   dbug << setw(Width_) << E1CachedValue_ << endl;
}

double DFTExternal::E0() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
      CallCount_[1]++;
   }
   
   return E0CachedValue_;
}

Vector const& DFTExternal::E1() const
{
   if (!Cached_[2])
   {
      UpdateValues(NoStiffness);
      CallCount_[2]++;
   }

   return E1CachedValue_;
}

Vector const& DFTExternal::E1DLoad() const
{
   if (!Cached_[3])
   {
      // calculate pressure terms.
      Matrix U(3,3);
      U[0][0] = 1.0+DOF_[0];
      U[1][1] = 1.0+DOF_[1];
      U[2][2] = 1.0+DOF_[2];
      U[1][2]=U[2][1] = DOF_[3];
      U[2][0]=U[0][2] = DOF_[4];
      U[0][1]=U[1][0] = DOF_[5];
      double UDet = U.Det();
      Matrix PressureTerm = UDet*U.Inverse();
      E1DLoadCachedValue_.Resize(DOFS_,0.0);
      E1DLoadCachedValue_[0] = -PressureTerm[0][0];
      E1DLoadCachedValue_[1] = -PressureTerm[1][1];
      E1DLoadCachedValue_[2] = -PressureTerm[2][2];
      E1DLoadCachedValue_[3] = -(PressureTerm[1][2]+PressureTerm[2][1]);
      E1DLoadCachedValue_[4] = -(PressureTerm[2][0]+PressureTerm[0][2]);
      E1DLoadCachedValue_[5] = -(PressureTerm[0][1]+PressureTerm[1][0]);
      
      Cached_[3] = 1;
      CallCount_[3]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& DFTExternal::E2() const
{
   if (!Cached_[4])
   {
      UpdateValues(NeedStiffness);
      CallCount_[4]++;
   }
   
   return E2CachedValue_;
}

void DFTExternal::Print(ostream& out,PrintDetail const& flag,
                        PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   double mintestfunct;
   double J;
   Matrix
      stiff(DOFS_,DOFS_);
   Vector str(DOFS_);
   Vector TestFunctVals(NumTestFunctions());
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   stiff = E2();
   str = E1();
   engy = E0();
   J = (1.0+DOF_[0])*( (1.0+DOF_[1])*(1.0+DOF_[2]) - DOF_[3]*DOF_[3] )
      -(DOF_[5])*( (1.0+DOF_[0])*(1.0+DOF_[2]) - DOF_[4]*DOF_[4] )
      +(DOF_[4])*( (1.0+DOF_[0])*(1.0+DOF_[1]) - DOF_[5]*DOF_[5] );
   
   NoNegTestFunctions=TestFunctions(TestFunctVals,LHS);
   mintestfunct = TestFunctVals[0];
   for (int i=0;i<DOFS_;++i)
   {
      if (mintestfunct > TestFunctVals[i])
         mintestfunct = TestFunctVals[i];
   }
   
   switch (flag)
   {
      case PrintLong:
         out << "DFTExternal:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "DFTExternal:" << "\n" << "\n";
         }
         // passthrough to short
      case PrintShort:
         out << "Lambda: " << setw(W) << Lambda_ << "\n"
             << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
             << "J (V/Vo) :" << setw(W) << J << "\n"
             << "DFT Energy Value:" << setw(W) << DFTEnergyCachedValue_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Stress:" << "\n" << setw(W) << str << "\n\n"
             << "Stiffness:" << setw(W) << stiff
             << "Eigenvalue Info:"  << "\n"<< setw(W) << TestFunctVals << "\n"
             << "Bifurcation Info:" << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda: " << setw(W) << Lambda_ << "\n"
                 << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                 << "J (V/Vo) :" << setw(W) << J << "\n"
                 << "DFT Energy Value:" << setw(W) << DFTEnergyCachedValue_ << "\n"
                 << "Potential Value:" << setw(W) << engy << "\n";

            cout << "Stress:" << "\n" << setw(W) << str << "\n\n"
                 << "Stiffness:" << setw(W) << stiff
                 << "Eigenvalue Info:"  << "\n" << setw(W) << TestFunctVals <<"\n"
                 << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }

   // insert blank line in dbug file to identify each record
   dbug << endl;
}

ostream& operator<<(ostream& out,DFTExternal& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
