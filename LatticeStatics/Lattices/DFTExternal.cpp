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
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
   dbug.close();
}

DFTExternal::DFTExternal(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input),
   Lambda_(0.0),
   Echo_(Echo),
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
      q = (q+1)%4;
      out << setw(30) << DOF_[i];
      if (q==3) out << "\n";
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
   cerr << "DFTExternal system() call returned with id: " << retid << endl;

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
   double PressureEnergy = Lambda_*UDet;
   Matrix PressureTerm = PressureEnergy*U.Inverse();
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

   in >> E0CachedValue_;
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
   Cached_[0] = 1;

   Matrix H(9,9);
   for (int k=0;k<3;++k)
      for (int l=0;l<3;++l)
         for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
            {
               H[3*i+j][3*k+l] = Del[k][i]*Del[l][j] + Del[k][i]*B[j][l];
            }
   // get stresses
   Matrix tmp(3,3);
   in >> tmp;
   Vector Tmp(9);
   for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
      {
         Tmp[3*i+j] = tmp[i][j];
      }
   Tmp = SolvePLU(H,Tmp);
   // initialize
   E1CachedValue_.Resize(DOFS_,0.0);
   for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
      {
         int a = ((i==j)? i : (6-i+j));
         E1CachedValue_[a] += UDet*Tmp[3*i+j];
      }
   // get forces
   for (int i=6;i<DOFS_;++i)
   {
      in >> E1CachedValue_[i];
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
   Cached_[1] = 1;
   
   if (flag==NeedStiffness)
   {
      Matrix CSym(6,6);
      in >> CSym;
      Matrix C(9,9);
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
               for (int l=0;l<3;++l)
               {
                  int a = ((i==j)? i : (6-(i+j)));
                  int b = ((k==l)? k : (6-(k+l)));
                  C[3*i+j][3*k+l] = CSym[a][b];
               }
      Matrix HInv = H.Inverse();
      C = HInv*C*(HInv.Transpose());
      Matrix L(6,6,0.0);
      for (int i=0;i<3;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
               for (int l=0;l<3;++l)
               {
                  int a = ((i==j)? i : (6-(i+j)));
                  int b = ((k==l)? k : (6-(k+l)));
                  L[a][b] += UDet*C[3*i+j][3*k+l];
               }
      for (int i=0;i<6;++i)
         for (int j=0;j<6;++j)
         {
            E2CachedValue_[i][j] = L[i][j];
         }
      Matrix D(DOFS_-6,6);
      in >> D;
      Matrix DTmp(DOFS_-6,9);
      for (int i=0;i<DOFS_-6;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
            {
               int a = ((j==k)? j : (6-(j+k)));
               DTmp[i][3*j+k] = D[i][a];
            }
      DTmp = UDet*DTmp*(HInv.Transpose());
      Matrix DFinal(DOFS_-6,6,0.0);
      for (int i=0;i<DOFS_-6;++i)
         for (int j=0;j<3;++j)
            for (int k=0;k<3;++k)
            {
               int a = ((j==k)? j : (6-(j+k)));
               DFinal[i][a] += DTmp[i][3*j+k];
            }
      for (int i=6;i<DOFS_;++i)
         for (int j=0;j<6;++j)
         {
            E2CachedValue_[i][j] = E2CachedValue_[j][i] = DFinal[i-6][j];
         }
      Matrix P(DOFS_-6,DOFS_-6);
      in >> P;
      for (int i=6;i<DOFS_;++i)
         for (int j=6;j<DOFS_;++j)
         {
            E2CachedValue_[i][j] = P[i-6][j-6];
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
      Cached_[3] = 1;
   }
   in.close();

   //output force/dof data
   dbug << setw(Width_) << DOF_ << setw(Width_) << Lambda_  << setw(Width_) << NoStiffness << endl;
   dbug << setw(Width_) << E1CachedValue_ << endl;
}

double DFTExternal::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues(NoStiffness);
      CallCount_[0]++;
   }
   
   return E0CachedValue_;
}

Vector const& DFTExternal::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& DFTExternal::E1DLoad() const
{
   if (!Cached_[2])
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
      
      Cached_[2] = 1;
      CallCount_[2]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& DFTExternal::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
      CallCount_[3]++;
   }
   
   return E2CachedValue_;
}

void DFTExternal::Print(ostream& out,PrintDetail const& flag)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   double mintestfunct;
   Matrix
      stiff(DOFS_,DOFS_);
   Vector str(DOFS_);
   Vector TestFunctVals(DOFS_);
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   stiff = E2();
   str = E1();
   engy = E0();
   
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
