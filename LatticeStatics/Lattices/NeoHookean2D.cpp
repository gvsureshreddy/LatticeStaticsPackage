#include "NeoHookean2D.h"

namespace neo_hookean
{
    void run();
    std::size_t get_system_size();
    void set_solution(double const* const solution);
    void get_rhs_and_tangent(double* const rhs, double* const tm);
    float get_energy();
}


NeoHookean2D::~NeoHookean2D()
{
    cout << "TwoBarTruss Function Calls:\n"
            << "\tE0 calls - " << CallCount_[0] << "\n"
            << "\tE1 calls - " << CallCount_[1] << "\n"
            << "\tE1DLoad calls - " << CallCount_[2] << "\n"
            << "\tE2 calls - " << CallCount_[3] << "\n";
}

NeoHookean2D::NeoHookean2D(PerlInput const& Input, int const& Echo, int const& Width) :
Lattice(Input, Echo),
        Lambda_(0.0),
        Width_(Width)
{
    LoadParameter_ = Load;
    for (int i = 0; i < cachesize; ++i)
    {
        Cached_[i] = 0;
        CallCount_[i] = 0;
    }
    Caching_ = 0;
    
    
    system_size_ = neo_hookean::get_system_size();
    std::cout << "NeoHookean2D size is " << system_size_ << std::endl;
    DOF_.Resize(system_size_,0.0);
    RHS_.Resize(system_size_,0.0);
    Stiff_.Resize(system_size_,system_size_,0.0);
    neo_hookean::set_solution(&(DOF_[0]));
    neo_hookean::get_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]));
    //std::cout << setw(20) << RHS_;
    //std::cout << setw(20) << Stiff_;
}

double NeoHookean2D::E0() const
{
    if ((!Caching_) || (!Cached_[0]))
    {
        E0CachedValue_ = neo_hookean::get_energy();
        Cached_[0] = 1;
        CallCount_[0]++;
    }
    
    return E0CachedValue_;
}

Vector const& NeoHookean2D::E1() const
{
    neo_hookean::get_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]));
    return RHS_[0];
}

Vector const& NeoHookean2D::E1DLoad() const
{
    return EmptyV_;
}

Matrix const& NeoHookean2D::E2() const
{
    neo_hookean::get_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]));
    return Stiff_[0][0];
}

void NeoHookean2D::ExtraTestFunctions(Vector& TF) const
{
    return;
}



void NeoHookean2D::Print(ostream& out, PrintDetail const& flag,
        PrintPathSolutionType const& SolType)
{
    return;
}

ostream& operator<<(ostream& out, NeoHookean2D& A)
{
    A.Print(out, Lattice::PrintShort);
    return out;
}
