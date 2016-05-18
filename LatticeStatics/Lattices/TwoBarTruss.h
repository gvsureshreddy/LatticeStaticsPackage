#ifndef RSE__TwoBarTruss
#define RSE__TwoBarTruss

#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class TwoBarTruss : public Lattice
{
private:
    int DOFS_;
    
    // DOF[i] = [u v]
    Vector DOF_;
    enum LDeriv {L0, DL};
    double Lambda_;  // applied load
    double Gamma_;   // imperfection in modulus
    double Theta_;
    double COSTheta_;
    double SINTheta_;
    
    int Width_;
    
    int Caching_;
    static const int cachesize = 6;
    mutable int Cached_[cachesize];
    mutable double E0CachedValue_;
    mutable Vector E1CachedValue_;
    mutable Vector E1DLoadCachedValue_;
    mutable Matrix E2CachedValue_;
    mutable Matrix E3CachedValue_;
    mutable Matrix E4CachedValue_;
    mutable Vector ExtraTestFunctions_;
    mutable int CallCount_[cachesize];
    
public:
    // Functions provided by TwoBarTruss
    TwoBarTruss(PerlInput const& Input, int const& Echo = 1, int const& Width = 20);
    ~TwoBarTruss();
    
    // Virtual Functions required by Lattice
    Vector const& DOF() const
    {
        return DOF_;
    }
    
    void SetDOF(Vector const& dof)
    {
        DOF_ = dof; if (Caching_)
        {
            for (int i = 0; i < cachesize; ++i)
            {
                Cached_[i] = 0;
            }
        }
    }
    
    double Lambda() const
    {
        return Lambda_;
    }
    
    void SetLambda(double const& lambda)
    {
        Lambda_ = lambda; if (Caching_)
        {
            for (int i = 0; i < cachesize; ++i)
            {
                Cached_[i] = 0;
            }
        }
    }
    
    virtual double E0() const;
    virtual Vector const& E1() const;
    virtual Vector const& E1DLoad() const;
    virtual Matrix const& E2() const;
    virtual Matrix const& E3() const;
    virtual Matrix const& E4() const;
    virtual void ExtraTestFunctions(Vector& TF) const;
    virtual char const* const Type() const
    {
        return "TwoBarTruss";
    }
    
    virtual void Print(ostream& out, PrintDetail const& flag,
        PrintPathSolutionType const& SolType = RegularPt);
    
    virtual void PrintPath(ostream& out, ostream& pathout, int const& width);
    
    friend ostream& operator<<(ostream& out, TwoBarTruss& A);
    
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
    
    Vector const& StressDL() const
    {
        return E1DLoad();
    }
    
    Matrix const& StiffnessDL() const
    {
        return EmptyM_;
    }
    
    virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
    {
    }
    
    virtual void SetGridSize(int const& Grid)
    {
    }
    
private:
    // temp storage space
    mutable double eps1_;
    mutable double eps2_;
    mutable double eps1u_;
    mutable double eps2u_;
    mutable double eps1v_;
    mutable double eps2v_;
    mutable double eps1uu_;
    mutable double eps1vv_;
    mutable double eps1uv_;
    mutable double eps2uu_;
    mutable double eps2vv_;
    mutable double eps2uv_;
    
    // place holder
    Vector EmptyV_;
    Matrix EmptyM_;
};

#endif
