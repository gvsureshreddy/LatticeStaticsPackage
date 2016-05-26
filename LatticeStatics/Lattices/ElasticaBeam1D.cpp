#include "ElasticaBeam1D.h"
#include <iostream>
#include <fstream>

namespace elastica_beam
{
  void run();
  std::size_t get_system_size();
  unsigned int get_unconstrained_system_size();
  void set_solution(double const* const solution);
  void get_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_unconstrained_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_E1DLoad(double* const E1Dload);
  double get_energy();
  void set_P(const double value_P);
}


ElasticaBeam1D::~ElasticaBeam1D()
{
  cout << "TwoBarTruss Function Calls:\n"
       << "\tE0 calls - " << CallCount_[0] << "\n"
       << "\tE1 calls - " << CallCount_[1] << "\n"
       << "\tE1DLoad calls - " << CallCount_[2] << "\n"
       << "\tE2 calls - " << CallCount_[3] << "\n";
}

ElasticaBeam1D::ElasticaBeam1D(PerlInput const& Input, int const& Echo, int const& Width) :
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


  system_size_ = elastica_beam::get_system_size();
  unconstrained_system_size_ = elastica_beam::get_unconstrained_system_size();
//  std::cout << "ElasticaBeam1D size is " << system_size_ << std::endl;
//  std::cout << "ElasticaBeam1D unconstrained size is "
//          << unconstrained_system_size_ << std::endl;
  DOF_.Resize(unconstrained_system_size_,0.0);
  RHS_.Resize(unconstrained_system_size_,0.0);
  E1DLoad_.Resize(unconstrained_system_size_,0.0);
  Stiff_.Resize(unconstrained_system_size_,unconstrained_system_size_,0.0);
  elastica_beam::set_solution(&(DOF_[0]));
  elastica_beam::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),0);
  //std::cout << setw(20) << RHS_;
  //std::cout << std::endl << std::endl << "Tangent matrix :\n\n\n" << std::endl;
  //std::cout << setw(20) << Stiff_;
}

void ElasticaBeam1D::SetLambda(double const& lambda)
{
    Lambda_ = lambda;
    elastica_beam::set_P(lambda);
    for (int i = 0; i < cachesize; ++i)
      Cached_[i] = 0;
}

void ElasticaBeam1D::SetDOF(Vector const& dof)
{
    DOF_ = dof;
    elastica_beam::set_solution(&(DOF_[0]));
}

double ElasticaBeam1D::E0() const
{
  if ((!Caching_) || (!Cached_[0]))
    {
      E0CachedValue_ = elastica_beam::get_energy();
      Cached_[0] = 1;
      CallCount_[0]++;
    }

  return E0CachedValue_;
}

Vector const& ElasticaBeam1D::E1() const
{
  elastica_beam::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),0);
  return RHS_;
}

Vector const& ElasticaBeam1D::E1DLoad() const
{
  elastica_beam::get_E1DLoad(&(E1DLoad_[0]));
  return E1DLoad_;
}

Matrix const& ElasticaBeam1D::E2() const
{
  elastica_beam::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),0);
  return Stiff_;
}

void ElasticaBeam1D::ExtraTestFunctions(Vector& TF) const
{
  return;
}



void ElasticaBeam1D::Print(ostream& out, PrintDetail const& flag,
			 PrintPathSolutionType const& SolType)
{
  return;
}

ostream& operator<<(ostream& out, ElasticaBeam1D& A)
{
  A.Print(out, Lattice::PrintShort);
  return out;
}
