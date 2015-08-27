#ifndef RSE__FEAP
#define RSE__FEAP

#include <string>
#include <sstream>
#include "PerlInput.h"
#include "Lattice.h"
#include "CMatrix.h"

using namespace std;

class FEAP : public Lattice
{
private:
   mutable int DOFS_; // BFB number of degrees of freedom
   mutable int DOFS_F_; // FEAP number of degrees of freedom

   mutable Vector DOF_; // BFB degrees of freedom
   mutable Vector DOF_F_; // FEAP degrees of freedom

   int * BoundNodes_;
   int * PeriodicNodes_;
   int * InnerNodes_;

   mutable Matrix Jacobian_; // d(DOF_F)/d(DOF_)
   mutable Matrix FJacobian_; // d(DOF_F)/d(DOF_+1)
   mutable Matrix DispJacobian_; // d(DOF_F)/d(U)
   int** Map_ ; // Represents the sparse 3D array d2(DOF_F)/d(DOF_)2
   int** FMap_ ; // Represents the sparse 3D array d2(DOF_F)/d(DOF_+1)2

   mutable double Lambda_;
   enum LoadingType {DEAD_LOAD, PRESSURE_LOAD, DISPLACEMENT_CONTROL};
   mutable LoadingType LoadingType_; // 0 - dead-load; 1 - pressure-load; 2 - displacement control
   mutable Matrix Load_; //Bio stress tensor
   mutable Matrix U_; //Right stretch tensor
   mutable Matrix F_; //Deformation gradient
   mutable double StretchRatio_; //Stretch ratio giving U_{22} = StretchRatio_*U_{11}

   int** N_; // Contains info for Dynamical Stiffness matrix. Only half the neighbouring cells is needed
              // Each row i contains : N_[i][0] & N_[i][1] : LV coordinates of neighbour cell
              //                       N_[i][k] & N_[i][l] : pairs of interacting nodes (reference cell node first) (FEAP numbering)
   int N_rows_;
   int N_cols_;
   mutable Vector K_; // Bloch wave vector
   mutable Vector KDirection_; //Direction for

   mutable CMatrix Dk_; //Dynamical matrix (Bloch-wave stiffness matrix)

   int KSpaceResolution_; //Used in Bloch wave analysis, when AnalysisType => Full or KDirection
   Vector KRange_; //Range of KVectors when AnalysisType => KDirection
   int NumKVectors_;
   Matrix KVectorMatrix_;//KVectors when AnalysisType => KVectors
   double TFLoad_;//Load when ExtraTestFunctions => LoadingParameter
   int TFType_;
   void DynamicalMatrix(Vector const& K) const; //Dynamical Matrix, uses the FEAP stiffness matrix
   void DynamicalMatrixBis(Vector const& K) const; //Dynamical Matrix, uses the BFB stiffness matrix


   char const* ffin_; // FEAP input file name
   int ndf_; // number of DOFs per node
   int ndm_; // number of spatial dimensions
   int nuc_; //number of unit cells
   int numnp_; // number of nodes in mesh
   int nel_; // number of elements in mesh
   int nen1_; // number of nodes per element
   int neq_; // number of reduced equations
   int* eqnID_; // equation number ID array
   int* elmConn_; // element connectivity array
   int* bcID_; // displacement boundary condition ID array
   Vector X_; // Reference coordinates of nodes
   Vector X_F_; // Reference coordinates of nodes in FEAP DOF_F_ form (with 0 for theta)
   int nbn_; // number of boundary nodes
   double eps_; //eps of penalty term for translation in energy
   double CellArea_; //Area of unit cell
   double HexSize_; // Hexagon diagonal length

   int Width_;
   double Tolerance_;

   enum UpdateFlag {NoStiffness = 0, NeedStiffness = 1};
   void UpdateValues(UpdateFlag flag) const; //Updates energy, first and second derivatives od energy
   void MapHelper(int const i, int const offset, int const* const Nodes);

   void UpdateDOF_F() const; //Uses the BFB DOFs and maps to FEAP DOFs
   void UpdateJacobian() const;
   void JacobianHelper(int const ii, int const jj, int const shift) const;

   double RankOneConvex(Matrix const& d2WdFdF) const;

   void print_gpl_config(fstream& out) const;
   fstream config_out_; //Prints out in separate file the configuration for plotting
   fstream critical_eig_out_; //Prints out in separate file eigenmode of stiff matrix when successfully bissected
   mutable fstream bloch_wave_out_; //Stream, prints out in separate file bloch wave analysis results. Columns are: k1, k2, minimum eigval, second minimum eigval
   mutable int config_count_; //Index count for config_out_ stream
   mutable int eig_count_; //Index count for critical_eig_out_
   mutable int bloch_count_; //Index count for bloch_wave_out_ stream
   fstream plot_out_; //Prints out in separate files plotting information. columns are Lambda, U11, U22, U12, E0
   mutable int plot_count_; //Index count for plot_out_ stream

   static const int cachesize = 4;
   mutable int Cached_[cachesize];
   mutable double E0CachedValue_;
   mutable Vector E1CachedValue_;
   mutable Vector W1CachedValue_;
   mutable Vector DispE1CachedValue_; // store conjugate to disp. control
   mutable Vector E1CachedValue_F_; //E1 Cached value for FEAP
   mutable Vector E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable Matrix W2CachedValue_;
   mutable Matrix E2CachedValue_F_;
   mutable int EvaluationCount_[2];
   mutable int CallCount_[cachesize];

public:
   // Functions provided by FEAP
   FEAP(PerlInput const& Input, int const& Echo = 1, int const& Width = 20);
   ~FEAP();

   // Virtual Functions required by Lattice
   Vector const& DOF() const
   {
      return DOF_;
   }

   void SetDOF(Vector const& dof)
   {
      DOF_ = dof;
      switch (LoadingType_)
      {
        case DEAD_LOAD:
        case PRESSURE_LOAD:
          U_[0][0] = DOF_[0];
          U_[1][1] = DOF_[1];
          U_[0][1] = DOF_[2]/sqrt(2.0);
          U_[1][0] = U_[0][1];
          F_ = U_;
          break;
        case DISPLACEMENT_CONTROL:
          U_[0][0] = Lambda_;
          U_[1][1] = StretchRatio_ * Lambda_;
          U_[0][1] = 0.0;
          U_[1][0] = 0.0;
          F_ = U_;
          break;
      }
      for (int i = 0; i < cachesize; ++i)
      {
         Cached_[i] = 0;
      }
   }

   double Lambda() const
   {
      return Lambda_;
   }

   void SetLambda(double const& lambda)
   {
      Lambda_ = lambda;
      if (LoadingType_ == DISPLACEMENT_CONTROL)
      {
        U_[0][0] = Lambda_;
        U_[1][1] = StretchRatio_ * Lambda_;
        U_[0][1] = 0.0;
        U_[1][0] = 0.0;
        F_ = U_;
      }
      for (int i = 0; i < cachesize; ++i)
      {
         Cached_[i] = 0;
      }
   }

   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Vector const& E1DLoad() const;
   virtual Vector const& StressDL() const
   {
      return E1DLoad();
   }

   virtual Matrix const& E2() const;
   virtual Matrix const& StiffnessDL() const;
   virtual char const* const Type() const
   {
      return "FEAP";
   }

   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);
   virtual int CriticalPointInfo(int* const CPCrossingNum, int const& TFIndex,
                                 Vector const& DrDt, int const& CPorBif,
                                 int const& NumZeroEigenVals, double const& Tolerance,
                                 int const& Width, PerlInput const& Input, ostream& out);
   virtual void ExtraTestFunctions(Vector& TF) const;


   friend ostream& operator<<(ostream& out, FEAP& A);

   // ignore these
   double Entropy() const
   {
      return 0.0;
   }

   double HeatCapacity() const
   {
      return 0.0;
   }

   Vector const& StressDT() const
   {
      return EmptyV_;
   }

   Matrix const& StiffnessDT() const
   {
      return EmptyM_;
   }

   double Temp() const
   {
      return 0.0;
   }

   void SetTemp(double const& Ntemp)
   {
   }

   virtual Matrix const& E3() const
   {
      cerr << "FEAP::E3() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual Matrix const& E4() const
   {
      cerr << "FEAP::E4() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
   {
   }

   virtual void SetGridSize(int const& Grid)
   {
   }

private:
   // statice for StiffnessDL
   mutable Matrix stiffdl_static;

   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;
};

#endif
