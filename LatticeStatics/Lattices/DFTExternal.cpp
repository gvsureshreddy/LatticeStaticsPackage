#include "DFTExternal.h"
#include <fstream>

using namespace std;

const double DFTExternal::Alt[3][3][3] = {{{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 1.0},
                                           {0.0, -1.0, 0.0}},
                                          {{0.0, 0.0, -1.0},
                                           {0.0, 0.0, 0.0},
                                           {1.0, 0.0, 0.0}},
                                          {{0.0, 1.0, 0.0},
                                           {-1.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}}};

DFTExternal::~DFTExternal()
{
   cout << "DFTExternal Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

DFTExternal::DFTExternal(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input),
   Lambda_(0.0),
   Echo_(Echo),
   Width_(Width)
{
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
   if (flag==NeedStiffness)
   {
      out << 1 << "\n";
   }
   else
   {
      out << 0 << "\n";
   }
   for (int i=0;i<6;++i)
   {
      out << setw(30) << DOF_[i];
   }
   out << "\n";
   for (int i=6;i<DOFS_;++i)
   {
      out << setw(30) << DOF_[i];
   }
   out << "\n";

   out.close();

   // update with correct command.
   retid=system("DFTModel DFTModelInput DFTModelOutput >& /dev/null");
   cerr << "DFTExternal system() call returned with id: " << retid << endl;

   // calculate pressure terms.
   Matrix U(3,3);
   U[0][0] = 1.0+DOF_[0];
   U[1][1] = 1.0+DOF_[1];
   U[2][2] = 1.0+DOF_[2];
   U[1][2]=U[2][1] = DOF_[3];
   U[2][0]=U[0][2] = DOF_[4];
   U[0][1]=U[1][0] = DOF_[5];
   double PressureEnergy = Lambda_*U.Det();
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
   
   in.open("DFTModelOutput",ios::in);
   if (in.fail())
   {
      cerr << "Error: Unable to open file : " << "DFTModelOutput " << "for read" << "\n";
      exit(-2);
   }

   in >> E0CachedValue_;
   E0CachedValue_ += PressureEnergy;
   Cached_[0] = 1;

   for (int i=0;i<6;++i)
   {
      in >> E1CachedValue_[i];
   }
   for (int i=6;i<DOFS_;++i)
   {
      in >> E1CachedValue_[i];
   }
   E1CachedValue_ += PressureStress;
   Cached_[1] = 1;
   
   if (flag==NeedStiffness)
   {
      in >> E2CachedValue_;
      E2CachedValue_ += PressureStiffness;
      Cached_[3] = 1;
   }

   in.close();
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
      Matrix PressureTerm = U.Det()*U.Inverse();
      E1DLoadCachedValue_.Resize(DOFS_,0.0);
      E1DLoadCachedValue_[0] = PressureTerm[0][0];
      E1DLoadCachedValue_[1] = PressureTerm[1][1];
      E1DLoadCachedValue_[2] = PressureTerm[2][2];
      E1DLoadCachedValue_[3] = PressureTerm[1][2]+PressureTerm[2][1];
      E1DLoadCachedValue_[4] = PressureTerm[2][0]+PressureTerm[0][2];
      E1DLoadCachedValue_[5] = PressureTerm[0][1]+PressureTerm[1][0];
      
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
}

ostream& operator<<(ostream& out,DFTExternal& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
