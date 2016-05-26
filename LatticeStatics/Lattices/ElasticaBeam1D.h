#ifndef RSE__ElasticaBeam1D
#define RSE__ElasticaBeam1D

#include <string>
#include <sstream>
#include <cstdlib>
#include "Lattice.h"

using namespace std;

class ElasticaBeam1D : public Lattice
{
    private:
        mutable int DOFS_;
        
        mutable Vector DOF_;
        mutable Vector RHS_;
        mutable Matrix Stiff_;
        mutable Vector E1DLoad_;
        mutable double Lambda_;
        
        int Width_;
        std::size_t system_size_;
        unsigned int unconstrained_system_size_;
        
        int Caching_;
        static const int cachesize = 4;
        mutable int Cached_[cachesize];
        mutable double E0CachedValue_;
        mutable Vector E1CachedValue_;
        mutable Vector E1DLoadCachedValue_;
        mutable Matrix E2CachedValue_;
        mutable int CallCount_[cachesize];
        
        public:
            // Functions provided by ElasticaBeam1D
            ElasticaBeam1D(PerlInput const& Input, int const& Echo = 1,
                    int const& Width = 20);
            ~ElasticaBeam1D();
            
            // Virtual Functions required by Lattice
            Vector const& DOF() const
            {
                return DOF_;
            }
            
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
                return "ElasticaBeam1D";
            }
            
            virtual void Print(ostream& out, PrintDetail const& flag,
                    PrintPathSolutionType const& SolType = RegularPt);
            
            virtual void PrintPath(ostream& out, ostream& pathout, int const& width);
            
            friend ostream& operator<<(ostream& out, ElasticaBeam1D& A);
            
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
                cerr << "ElasticaBeam1D::E3() Not Programmed\n"; exit(-1); return EmptyM_;
            }
            
            virtual Matrix const& E4() const
            {
                cerr << "ElasticaBeam1D::E4() Not Programmed\n"; exit(-1); return EmptyM_;
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
