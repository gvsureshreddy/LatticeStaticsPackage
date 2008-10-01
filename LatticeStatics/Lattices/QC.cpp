#include <fstream>
#include <string>
#include <sstream>
#include "QC.h"

using namespace std;

extern "C" void qcbfb_energy_(int& mode,int& nfree,double* u,double& t,double& E,double* Eu,
                             double* Euu,double* Eut);
extern "C" void qcbfb_restart_(char* filename);

QC::~QC()
{
   cout << "QC Function Calls:\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

QC::QC(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input),
   Lambda_(0.0),
   Echo_(Echo),
   Width_(Width)
{
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash,"QC");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   DOF_.Resize(DOFS_,0.0);
   E1CachedValue_.Resize(DOFS_);
   E1DLoadCachedValue_.Resize(DOFS_);
   E2CachedValue_.Resize(DOFS_,DOFS_);
   EmptyV_.Resize(DOFS_,0.0);
   EmptyM_.Resize(DOFS_,DOFS_,0.0);

   LoadParameter_ = Load;
   for (int i=0;i<cachesize;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
}

void QC::UpdateValues(UpdateFlag flag) const
{
   if (NoStiffness==flag)
   {
      int mode=0;
      qcbfb_energy_(mode,DOFS_,&(DOF_[0]),Lambda_,E0CachedValue_,&(E1CachedValue_[0]),0,0);
      Cached_[0]=1;
      Cached_[1]=1;
   }
   else if (NeedStiffness==flag)
   {
      int mode=1;
      qcbfb_energy_(mode,DOFS_,&(DOF_[0]),Lambda_,E0CachedValue_,&(E1CachedValue_[0]),
                   &(E2CachedValue_[0][0]),&(E1DLoadCachedValue_[0]));
      Cached_[0]=1;
      Cached_[1]=1;
      Cached_[2]=1;
      Cached_[3]=1;
   }
   else
   {
      cerr << "Error in QC::UpdateValues(), unknown UpdateFlag.\n";
      exit(-45);
   }
}

double QC::E0() const
{
   if (!Cached_[0])
   {
      UpdateValues(NoStiffness);
      CallCount_[0]++;
   }
   
   return E0CachedValue_;
}

Vector const& QC::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
      CallCount_[1]++;
   }

   return E1CachedValue_;
}

Vector const& QC::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues(NeedStiffness);
      CallCount_[2]++;
   }
   
   return E1DLoadCachedValue_;
}

Matrix const& QC::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
      CallCount_[3]++;
   }
   
   return E2CachedValue_;
}

int QC::CriticalPointInfo(int const& CPCrossingNum,char const& CPSubNum,
                          Vector const& DrDt,int const& NumZeroEigenVals,
                          double const& Tolerance,int const& Width,
                          PerlInput const& Input,ostream& out)
{
   Matrix
      D2=E2(),
      EigVec,
      EigVal=SymEigVal(D2,&EigVec);
   Vector D1T(D2.Cols());
   int Bif = 2;
   
   D1T=StressDL();
   
   int dofs=D2.Rows();
   
   Matrix Mode;

   out << "Critical Point Crossing Number: " << CPCrossingNum << CPSubNum << "\n";
   if (Echo_)
   {
      cout << "Critical Point Crossing Number: " << CPCrossingNum << CPSubNum << "\n";
   }
   
   // Find the modes
   int count = 0,
      Ind[DOFMAX];
   
   for (int i=0;i<dofs;i++)
   {
      if (fabs(EigVal[0][i]) < Tolerance)
      {
         Ind[count++]=i;
      }
   }

   // Check for incorrect number of modes
   if (count != NumZeroEigenVals)
   {
      int skp;
      for (int j=count;j<NumZeroEigenVals;++j)
      {
         Ind[count] = 0;
         int a=0;
         while ((Ind[count] == Ind[a]) && (a < j))
         {
            (Ind[count])++;
            a++;
         }
         
         for (int i=1;i<dofs;++i)
         {
            skp=0;
            for (int k=0;k<count;++k)
            {
               if (Ind[k] == i) skp=1;
            }
            
            if (!skp)
            {
               if (fabs(EigVal[0][i]) < fabs(EigVal[0][Ind[count]]))
                  Ind[count] = i;
            }
         }
         count++;
      }
      
      out << "NOTE: Incorrect number of zero eigenvalues found. "
          << "Modes with smallest abs. value used." << "\n";
      if (Echo_)
      {
         cout << "NOTE: Incorrect number of zero eigenvalues found. "
              << "Modes with smallest abs. value used." << "\n";
      }
   }
   
   for (int i=0;i<count;++i)
      out << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
   out << "\n";
   if (Echo_)
   {
      for (int i=0;i<count;++i)
         cout << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
      cout << "\n";
   }
   
   if (BIFMAX < count)
   {
      cerr << "Error: BIFMAX < " << count << " in Lattice.h" << "\n";
      exit(-6);
   }
   
   Mode.Resize(count,dofs);
   
   for (int i=0;i<count;i++)
   {
      for (int j=0;j<dofs;j++)
      {
         Mode[i][j] = EigVec[j][Ind[i]];
      }
   }
   
   // Print out the critical Eigenvalues
   out << "EigenValues: ";
   if (Echo_) cout << "EigenValues: ";
   for (int i=0;i<count;++i)
   {
      out << setw(Width) << EigVal[0][Ind[i]];
      if (Echo_) cout << setw(Width) << EigVal[0][Ind[i]];
   }
   
   if (Echo_) cout << "\n";
   out << "\n";
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "-";
      out << "-";
   }
   if (Echo_) cout << endl;
   out << endl;
   
   // output a QC restart file
   ostringstream cpfilename;   
   cpfilename << Input.LastInputFileName() << ".CP." << setw(2) << setfill('0')
	      << CPCrossingNum << CPSubNum << ".res";
   char fortranstring[80];
   strcpy(fortranstring,cpfilename.str().c_str());
   for (int i=strlen(fortranstring);i<80;++i)
   {
      fortranstring[i] = ' ';
   }
   qcbfb_restart_(fortranstring);

   // output a new input file to help restart at this critical point
   cpfilename.str("");
   fstream cpfile;
   char tmp[2048];
   strcpy(tmp,Input.LastInputFileName());
   tmp[strlen(tmp)-3] = 0;
   cpfilename << tmp;
   if (2 == Bif)
      cpfilename << ".CP.";
   else if (1 == Bif)
      cpfilename << ".BP.";
   else
      cpfilename << ".TP.";
   cpfilename << setw(2) << setfill('0') << CPCrossingNum << CPSubNum << ".bfb";
   cpfile.open(cpfilename.str().c_str(),ios::out);

   cpfile << setprecision(out.precision()) << scientific;
   Vector T(dofs+1);
   Vector M(dofs+1);
   Vector dof=DOF();
   for (int i=0;i<dofs;++i)
   {
      T[i] = dof[i];
   }
   T[dofs] = ((LoadParameter_ == Temperature) ? Temp() : Lambda());

   cpfile << Input.ReconstructedInput();
   cpfile << "\n\n";
   
   Input.writeString(cpfile,"Bifurcation","StartType","Type");
   for (int i=0;i<count;++i)
   {
      for (int j=0;j<dofs;++j)
      {
         M[j] = Mode[i][j];
      }
      M[dofs] = 0.0;
      Input.writeVector(cpfile,M,"StartType","Tangent");
   }
   Input.writeVector(cpfile,T,"StartType","BifurcationPoint");
   cpfile.close();
   
   return Bif;
}


void QC::Print(ostream& out,PrintDetail const& flag)
{
   int W;
   int NoNegTestFunctions;
   double engy;
   double mintestfunct;
   Vector TestFunctVals(DOFS_);
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
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
         out << "QC:" << "\n" << "\n";

         if (Echo_)
         {
            cout << "QC:" << "\n" << "\n";
         }
         // passthrough to short
      case PrintShort:
         out << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "Potential Value:" << setw(W) << engy << "\n";

         out << "Bifurcation Info:" << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "Lambda (t): " << setw(W) << Lambda_ << "\n"
                 << "Potential Value:" << setw(W) << engy << "\n";

            cout << "Bifurcation Info:" << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         break;
   }
}

ostream& operator<<(ostream& out,QC& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
