#ifndef RSE__NeoHookean2D
#define RSE__NeoHookean2D

#include <string>
#include <sstream>
#include <cstdlib>
#include "Lattice.h"

class NeoHookean2D : public Lattice
{
    private:
        mutable int DOFS_; // Number of DOFs in BFB
        mutable int DOFS_D_; // Number of DOFs in deal.II
        mutable int dofs_vertical_;
        double factor_penalty_ = 1;
        
        mutable Vector DOF_; // DOFs in BFB
        mutable Vector DOF_D_; // DOFs in deal.II
        mutable Vector dofs_properties_; // Description of the deal.II DOFs
        mutable Vector constraint_properties_; // Description of the periodic properties (which DOF is linked with which)
        mutable Vector links_from_constrained_to_unconstrained_;
        mutable Vector RHS_;
        mutable Vector RHS_D_;
        mutable Matrix Stiff_;
        mutable Matrix Stiff_D_;
        mutable Vector E1DLoad_;
        mutable Vector E1DLoad_D_;
        mutable Vector dof_vertical_;
        mutable double Lambda_;
        
//        mutable Matrix Jacobian_; // d(DOF_D)/d(DOF_)
//        mutable Matrix FJacobian_; // d(DOF_D)/d(DOF_+1) //Bastien : what's that?
//        mutable Matrix DispJacobian_; // d(DOF_D)/d(U) //Bastien : what's that?
//        int** Map_ ; // Represents the sparse 3D array d2(DOF_F)/d(DOF_)2 //Bastien : I don't like the type
//        int** FMap_ ; // Represents the sparse 3D array d2(DOF_F)/d(DOF_+1)2 //Bastien : I don't like the type
        
        int Width_;
        //std::size_t system_size_;
        
        int Caching_;
        static const int cachesize = 4;
        mutable int Cached_[cachesize];
        mutable double E0CachedValue_;
        mutable Vector E1CachedValue_;
        mutable Vector E1DLoadCachedValue_;
        mutable Matrix E2CachedValue_;
        mutable int CallCount_[cachesize];
        
        public:
            // Functions provided by NeoHookean2D
            NeoHookean2D(PerlInput const& Input, int const& Echo = 1,
                    int const& Width = 20);
            ~NeoHookean2D();
            
            // Virtual Functions required by Lattice
            Vector const& DOF() const
            {
                return DOF_;
            }
            
            void fill_links_from_constrained_to_unconstrained();
            
            void SetDOF(Vector const& dof);
            
            double Lambda() const
            {
                return Lambda_;
            }
            
            void SetLambda(double const& lambda);
            
            virtual double E0() const;
            virtual Vector const& E1() const;
            virtual Vector const& E1DLoad() const;
            virtual Vector const& StressDL() const
            {
                return E1DLoad();
            }
            
            virtual Matrix const& E2() const;
            virtual Matrix const& StiffnessDL() const
            {
                return EmptyM_;
            }
            virtual void ExtraTestFunctions(Vector& TF) const;
            virtual char const* const Type() const
            {
                return "NeoHookean2D";
            }
            
            virtual void Print(ostream& out, PrintDetail const& flag,
                    PrintPathSolutionType const& SolType = RegularPt);
            
            virtual void PrintPath(ostream& out, ostream& pathout, int const& width);
            
            friend ostream& operator<<(ostream& out, NeoHookean2D& A);
            
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
                cerr << "NeoHookean2D::E3() Not Programmed\n"; exit(-1); return EmptyM_;
            }
            
            virtual Matrix const& E4() const
            {
                cerr << "NeoHookean2D::E4() Not Programmed\n"; exit(-1); return EmptyM_;
            }
            
            virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
            {
            }
            
            virtual void SetGridSize(int const& Grid)
            {
            }
            
            private:
                // statice for StiffnessDL and E3
                mutable Matrix stiffdl_static;
                mutable Matrix E3_static;
                
                // place holder
                Vector EmptyV_;
                Matrix EmptyM_;
};

#endif
