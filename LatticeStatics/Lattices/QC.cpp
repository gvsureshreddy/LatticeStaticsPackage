 #include <fstream>
#include "QC.h"

using namespace std;

extern "C" void qcbfb_energy_(int& mode,int& nfree,double* u,double& t,double& E,double* Eu,
                             double* Euu,double* Eut);
extern "C" void qcbfb_restart_(char* filename);
extern "C" void qcbfb_output_(int& nfree,double* u,double& prop,int& nint,int* intdata,int& ndouble,double* doubledata);

QC::~QC()
{
   cout << "QC Function Calls:\n"
        << "\tEvaluation - w/o stiffness - " << EvaluationCount_[0] << "\n"
        << "\tEvaluation - w   stiffness - " << EvaluationCount_[1] << "\n"
        << "\tE0 calls - " << CallCount_[0] << "\n"
        << "\tE1 calls - " << CallCount_[1] << "\n"
        << "\tE1DLoad calls - " << CallCount_[2] << "\n"
        << "\tE2 calls - " << CallCount_[3] << "\n";
}

QC::QC(PerlInput const& Input,int const& Echo,int const& Width):
   Lattice(Input,Echo),
   Lambda_(0.0),
   Width_(Width),
   SolutionNumber_(0)
{
   Stable_[0] = 1;
   Stable_[1] = 1;

   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash,"QC");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   RemoveTranslation_ = Input.getPosInt(Hash,"RemoveTranslation");// 0 -off, 1,2,3 - direction
   if (RemoveTranslation_ > 3)
   {
      cerr << "Error: QC - RemoveTranslation parameter too big must be 0, 1, 2, or 3.\n";
      exit(-34);
   }
   
   if (RemoveTranslation_ > 0)
   {
      // determine translation mode.
      char stor[2048];
      char const* const tabfile = Input.getString(Hash,"TranslationTable");
      fstream table(tabfile,ios::in);
      if (!table.is_open())
      {
         cerr << "Error: QC unable to open file " << tabfile << "!" << endl;
         exit(-12);
      }
      table.getline(stor,2048); // header line
      table.getline(stor,2048); // header line
      int bfbdof;
      int qcnode;
      int qcdir;

      TranslationMode_.Resize(DOFS_,0.0);
      for (int i=0;i<DOFS_;++i)
      {
         table >> bfbdof;
         table >> qcnode;
         table >> qcdir;
         table.getline(stor,2048); // rest of line
         if (qcdir == RemoveTranslation_)
         {
            TranslationMode_[bfbdof-1] = 1.0;
         }
      }
      table.close();

      TranslationMode_ /= TranslationMode_.Norm();
   }
   
   if (Input.ParameterOK(Hash,"Tolerance"))
   {
      Tolerance_ = Input.getDouble(Hash,"Tolerance");
   }
   else
   {
      Tolerance_ = Input.useDouble(1.0e-6,Hash,"Tolerance");  // Default Value
   }

   if (NumExtraTFs_ > 0)
   {
      if (Input.ParameterOK(Hash,"ExtraTFs"))
      {
         if (Input.getArrayLength(Hash,"ExtraTFs") == NumExtraTFs_)
         {
            ExtraTestFunctions_.Resize(NumExtraTFs_);
            PreviousExtraTestFunctions_.Resize(NumExtraTFs_);
            ExtraTestFunctionMultipliers_.Resize(NumExtraTFs_,1.0);
            Input.getVector(ExtraTestFunctions_,Hash,"ExtraTFs");
         }
         else
         {
            cerr << "Error: ArrayLength of " << Hash.Name
                 << "{ExtraTFs} is not equal to Lattice{NumExtraTFs}.\n";
            exit(-2);
         }
      }
      else
      {
         cerr << "Error: ExtraTFs not defined but Lattice{NumExtraTFs} = "
              << NumExtraTFs_ << ".\n";
         exit(-3);
      }
   }
   else
   {
      ExtraTestFunctions_.Resize(0);
      PreviousExtraTestFunctions_.Resize(0);
      ExtraTestFunctionMultipliers_.Resize(0);
   }

   // get input file header
   char tmp[2048];
   strcpy(tmp,Input.LastInputFileName());
   int len = strlen(tmp);
   tmp[len-3] = 'i';
   tmp[len-2] = 'n';
   tmp[len-1] = 0;
   fstream infile(tmp,ios::in);
   if (!infile.is_open())
   {
      cerr << "Error: QC unable to open file" << tmp << "!" << endl;
      exit(-12);
   }
   infile.getline(tmp,2048);
   while (strcmp("macros",tmp))
   {
      InFileHeader_ << tmp << "\n";
      infile.getline(tmp,2048);
   }
   infile.close();
   
   DOF_.Resize(DOFS_,0.0);
   E1CachedValue_.Resize(DOFS_);
   E1DLoadCachedValue_.Resize(DOFS_);
   E2CachedValue_.Resize(DOFS_,DOFS_);
   stiffdl_static.Resize(DOFS_,DOFS_);
   E3_static.Resize(DOFS_*DOFS_,DOFS_);
   EmptyV_.Resize(DOFS_,0.0);
   EmptyM_.Resize(DOFS_,DOFS_,0.0);

   LoadParameter_ = Load;
   for (int i=0;i<cachesize;++i)
   {
      Cached_[i] = 0;
      CallCount_[i] = 0;
   }
   EvaluationCount_[0] = 0;
   EvaluationCount_[1] = 0;
}

void QC::UpdateValues(UpdateFlag flag) const
{
   if (NoStiffness==flag)
   {
      int mode=0;
      qcbfb_energy_(mode,DOFS_,&(DOF_[0]),Lambda_,E0CachedValue_,&(E1CachedValue_[0]),0,0);
      Cached_[0]=1;
      Cached_[1]=1;
      EvaluationCount_[0]++;
      if (RemoveTranslation_ > 0)
      {
         double T=0.0;
         for (int i=0;i<DOFS_;++i)
         {
            T += TranslationMode_[i]*DOF_[i];
         }

         E0CachedValue_ += 0.5*T*T;
         for (int i=0;i<DOFS_;++i)
         {
            E1CachedValue_[i] += T*TranslationMode_[i];
         }
      }
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
      EvaluationCount_[1]++;

      if (RemoveTranslation_ > 0)
      {
         double T=0.0;
         for (int i=0;i<DOFS_;++i)
         {
            T += TranslationMode_[i]*DOF_[i];
         }

         E0CachedValue_ += 0.5*T*T;
         for (int i=0;i<DOFS_;++i)
         {
            E1CachedValue_[i] += T*TranslationMode_[i];
            for (int j=0;j<DOFS_;++j)
            {
               E2CachedValue_[i][j] += TranslationMode_[i]*TranslationMode_[j];
            }
         }
      }
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
   }
   CallCount_[0]++;
   
   return E0CachedValue_;
}

Vector const& QC::E1() const
{
   if (!Cached_[1])
   {
      UpdateValues(NoStiffness);
   }
   CallCount_[1]++;

   return E1CachedValue_;
}

Vector const& QC::E1DLoad() const
{
   if (!Cached_[2])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[2]++;
   
   return E1DLoadCachedValue_;
}

Matrix const& QC::E2() const
{
   if (!Cached_[3])
   {
      UpdateValues(NeedStiffness);
   }
   CallCount_[3]++;
   
   return E2CachedValue_;
}

Matrix const& QC::StiffnessDL() const
{
   double load = Lambda_;

   Lambda_ = load+10.0*Tolerance_; for (int i=0;i<cachesize;++i) Cached_[i]=0;
   stiffdl_static = E2();
   Lambda_ = load-10.0*Tolerance_; for (int i=0;i<cachesize;++i) Cached_[i]=0;
   stiffdl_static -= E2();
   stiffdl_static /= 2.0*Tolerance_;

   Lambda_ = load; for (int i=0;i<cachesize;++i) Cached_[i]=0;
   return stiffdl_static;
}

Matrix const& QC::E3() const
{
   Vector OrigDOF = DOF();
   Vector pert;

   for (int i=0;i<DOFS_;++i)
   {
      pert.Resize(DOFS_,0.0);
      pert[i] = Tolerance_;
      DOF_ = OrigDOF + pert; for (int r=0;r<cachesize;++r) Cached_[r]=0;
      {
         Matrix const& stiff = E2();
         for (int j=0;j<DOFS_;++j)
         {
            for (int k=0;k<DOFS_;++k)
            {
               E3_static[DOFS_*j+k][i] = stiff[j][k];
            }
         }
      }
      DOF_ = OrigDOF - pert; for (int r=0;r<cachesize;++r) Cached_[r]=0;
      {
         Matrix const& stiff = E2();
         for (int j=0;j<DOFS_;++j)
         {
            for (int k=0;k<DOFS_;++k)
            {
               E3_static[DOFS_*j+k][i] -= stiff[j][k];
            }
         }
      }
      for (int j=0;j<DOFS_;++j)
      {
         for (int k=0;k<DOFS_;++k)
         {
            E3_static[DOFS_*j+k][i] /= 2.0*Tolerance_;
         }
      }
   }

   DOF_ = OrigDOF; for (int r=0;r<cachesize;++r) Cached_[r]=0;

   return E3_static;
}

int QC::CriticalPointInfo(int* const CPCrossingNum,int const& TFIndex,Vector const& DrDt,
                          int const& CPorBif,int const& NumZeroEigenVals,double const& Tolerance,
                          int const& Width,PerlInput const& Input,ostream& out)
{
   int const IndexZeros = 4;
   int const OccuranceZeros = 3;

   int Bif;

   // do standard CPInfo stuff and output bfb restart file
   Bif = Lattice::CriticalPointInfo(CPCrossingNum,TFIndex,DrDt,CPorBif,NumZeroEigenVals,
                                    Tolerance,Width,Input,out);
   
   ostringstream cpfilename;
   if (1 == Bif)
      cpfilename << ".B";
   else if (0 == Bif)
      cpfilename << ".T";
   else
      cpfilename << ".E";

   // output a QC restart file
   ostringstream qcfilename;
   char tmp[2048];
   strcpy(tmp,Input.LastInputFileName());
   tmp[strlen(tmp)-4] = 0;
   qcfilename.fill('0');
   qcfilename << tmp << cpfilename.str() << setw(IndexZeros)
	      << TFIndex << "-" << setw(OccuranceZeros)
              << CPCrossingNum[TFIndex] << ".res";
   qcfilename.fill(' ');
   char fortranstring[80];
   strcpy(fortranstring,qcfilename.str().c_str());
   for (int i=strlen(fortranstring);i<80;++i)
   {
      fortranstring[i] = ' ';
   }
   qcbfb_restart_(fortranstring);

   // output a qc input file (if bif pt)
   ostringstream bfbfilename;
   bfbfilename.fill('0');
   bfbfilename << tmp << cpfilename.str() << setw(IndexZeros) << TFIndex
               << "-" << setw(OccuranceZeros) << CPCrossingNum[TFIndex];
   bfbfilename.fill(' ');
   if (1 == Bif)
   {
      fstream infile;
      strcpy(tmp,bfbfilename.str().c_str());
      int len = strlen(tmp);
      tmp[len] = '.';
      tmp[len+1] = 'i';
      tmp[len+2] = 'n';
      tmp[len+3] = 0;
      infile.open(tmp,ios::out);
      infile << "% Input file for: " << bfbfilename.str() << "\n";
      infile << InFileHeader_.str();
      infile << "macros\n";
      infile << "restart,read," << bfbfilename.str() << "\n";
      infile << "status\n";
      infile << "tole,,1.0d-6\n";
      infile << "proportional,,2,,-1000.,-1000.,1000.,1000.\n\n";
      infile << "% generate output for initial configuration\n";
      infile << "form\n";
      infile << "report\n\n";
      infile << "% Start bfb solution\n";
      infile << "bfb,rest," << bfbfilename.str() << "\n\n";
      infile << "loop,,100\n";
      infile << "   bfb\n";
      infile << "   touch,,\\#PROCESSING#\n";
      infile << "   conv,bfb\n";
      infile << "next\n\n";
      infile << "% End bfb solution (release memory)\n";
      infile << "bfb,term\n\n";
      infile << "end\n";
      infile << "stop\n";
      infile.close();
   }


   // touch sentinel file to indicate that all bif pt. files have been created.
   if (Bif == 1)
   {
      ostringstream sentinelfilename;
      sentinelfilename << "touch " << bfbfilename.str() << ".sentinel";
      system(sentinelfilename.str().c_str());
   }
   
   return Bif;
}

void QC::ExtraTestFunctions(Vector& TF) const
{
   if ((Stable_[0]) || (Stable_[1]))
   {
      for (int i=0;i<NumExtraTFs_;++i)
      {
         TF[i] = ExtraTestFunctionMultipliers_[i]*(ExtraTestFunctions_[i] - Lambda());
         PreviousExtraTestFunctions_[i] = TF[i];
      }
   }
   else
   {
      for (int i=0;i<NumExtraTFs_;++i)
      {
         TF[i] = ExtraTestFunctionMultipliers_[i]*(ExtraTestFunctions_[i] - Lambda());
         if (TF[i]*PreviousExtraTestFunctions_[i] < 0.0)
         {
            ExtraTestFunctionMultipliers_[i] *= -1.0;
            TF[i] = -TF[i];
         }
         PreviousExtraTestFunctions_[i] = TF[i];
      }
   }
}



void QC::Print(ostream& out,PrintDetail const& flag,
               PrintPathSolutionType const& SolType)
{
   int W;
   int NoNegTestFunctions=0;
   double engy;
   double E1norm;
   double mintestfunct;
   Vector TestFunctVals(NumTestFunctions());
   
   W=out.width();
   
   out.width(0);
   if (Echo_) cout.width(0);
   
   engy = E0();
   E1norm = E1().Norm();
   
   TestFunctions(TestFunctVals,LHS);
   mintestfunct = TestFunctVals[0];
   // check only the EigenValTFs
   for (int i=0;i<DOFS_;++i)
   {
      if ((UseEigenValTFs() == 1) && (TestFunctVals[i] < 0.0)) ++NoNegTestFunctions;
      if (mintestfunct > TestFunctVals[i])
         mintestfunct = TestFunctVals[i];
   }
   Stable_[1] = Stable_[0];
   if (NoNegTestFunctions == 0)
   {
      Stable_[0] = 1;
   }
   else
   {
      Stable_[0] = 0;
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
         out << "QC SolutionNumber = " << SolutionNumber_ << "\n"
             << "Lambda (t): " << setw(W) << Lambda_ << "\n"
             << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
             << "Potential Value: " << setw(W) << engy << "\n"
             << "Force Norm: " << setw(W) << E1norm << "\n";

         out << "ExtraTF Info: ";
         for (int i=DOFS_;i<NumTestFunctions();++i)
         {
            out << setw(W) << TestFunctVals[i];
         }
         out << "\n";
         out << "Bifurcation Info: " << setw(W) << mintestfunct
             << setw(W) << NoNegTestFunctions << "\n";
         // send to cout also
         if (Echo_)
         {
            cout << "QC SolutionNumber = " << SolutionNumber_ << "\n"
                 << "Lambda (t): " << setw(W) << Lambda_ << "\n"
                 << "DOF Norm: " << setw(W) << DOF_.Norm() << "\n"
                 << "Potential Value: " << setw(W) << engy << "\n"
                 << "Force Norm: " << setw(W) << E1norm << "\n";

            cout << "ExtraTF Info: ";
            for (int i=DOFS_;i<NumTestFunctions();++i)
            {
               cout << setw(W) << TestFunctVals[i];
            }
            cout << "\n";
            cout << "Bifurcation Info: " << setw(W) << mintestfunct
                 << setw(W) << NoNegTestFunctions << "\n";
         }
         ++SolutionNumber_;
         
         int nint = 1;
         int ndouble = 1;
         double dummy = 0.0;
         int tpflag = -1;
         int bifflag = -2;
         int extratfflag = -3;
         switch (SolType)
         {
            case NotSolutionPt:
               break;
            case RegularPt:
               qcbfb_output_(DOFS_,&(DOF_[0]),Lambda_,nint,&NoNegTestFunctions,ndouble,&dummy);
               break;
            case TurningPt:
               qcbfb_output_(DOFS_,&(DOF_[0]),Lambda_,nint,&tpflag,ndouble,&dummy);
               break;
            case BifurcationPt:
               qcbfb_output_(DOFS_,&(DOF_[0]),Lambda_,nint,&bifflag,ndouble,&dummy);
               break;
            case ExtraTFPt:
               qcbfb_output_(DOFS_,&(DOF_[0]),Lambda_,nint,&extratfflag,ndouble,&dummy);
               break;
         }
         break;
   }
}

ostream& operator<<(ostream& out,QC& A)
{
   A.Print(out,Lattice::PrintShort);
   return out;
}
