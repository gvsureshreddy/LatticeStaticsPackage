#include "ElasticaBeam1D.h"
#include <iostream>
#include <fstream>

namespace elastica_beam
{
  void createObject();
  void deleteObject();
  void run();
  std::size_t get_system_size();
  unsigned int get_unconstrained_system_size();
  void set_solution(double const* const solution);
  void get_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_unconstrained_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_E1DLoad(double* const E1Dload);
  double get_energy();
  void set_P(const double value_P);
  int NoNegTestFunctions;
}


ElasticaBeam1D::~ElasticaBeam1D()
{
  elastica_beam::deleteObject();
  cout << "ElasticaBeam Function Calls:\n"
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
  elastica_beam::createObject();
  LoadParameter_ = Load;
  for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
      CallCount_[i] = 0;
    }
  Caching_ = 0;


  system_size_ = elastica_beam::get_system_size();
  unconstrained_system_size_ = elastica_beam::get_unconstrained_system_size();
  DOFS_ = unconstrained_system_size_;
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
  return -RHS_;
}

Vector const& ElasticaBeam1D::E1DLoad() const
{
  elastica_beam::get_E1DLoad(&(E1DLoad_[0]));
  return -E1DLoad_;
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
    int W;
    elastica_beam::NoNegTestFunctions = 0;
    double engy;
    double mintestfunct;
    out << "\nDOFS_ = " << DOFS_  << "\n";
    Matrix stiff(DOFS_, DOFS_);
    Vector str(DOFS_);
    Vector TestFunctVals(NumTestFunctions());
    W = out.width();
    out.width(0);
    if (Echo_)
    {
        cout.width(0);
    }
    engy = E0();
    str = E1();
    stiff = E2();
    TestFunctions(TestFunctVals, LHS);
    mintestfunct = TestFunctVals[0];
    for (int i = 0; i < NumTestFunctions(); ++i)
    {
        if ((TestFunctVals[i] < 0.0) && (i < DOFS_))
        {
            ++elastica_beam::NoNegTestFunctions;
        }
        if (mintestfunct > TestFunctVals[i])
        {
            mintestfunct = TestFunctVals[i];
        }
    }

    switch (flag)
    {
        case PrintLong:
            out << "ElasticaBeam:" << "\n" << "\n";
            out << "I don't know what to print in print long!\n";

            if (Echo_)
            {
                cout << "ElasticaBeam:" << "\n" << "\n";
                cout << "I don't know what to print in print long!\n";
            }
            // passthrough to short
        case PrintShort:
            out << "\n__________________________________________\n\n"
                    << "Lambda: " << setw(W) << Lambda_ << "\n"
                    << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                    << "Potential Value:" << setw(W) << engy << "\n";

            out << "Stress:" << "\n" << setw(W) << str << "\n\n" //To be deleted
                    << "Stiffness:" << setw(W) << stiff //To be deleted
                    << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
                    << "Bifurcation Info:" << setw(W) << mintestfunct
                    << setw(W) << elastica_beam::NoNegTestFunctions << "\n";
            out << setw(W) << Lambda_ << " " << setw(W) << DOF_ << "\n";
            // send to cout also
            if (Echo_)
            {
                cout << "__________________________________________\n\n"
                        << "Lambda: " << setw(W) << Lambda_ << "\n"
                        << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                        << "Potential Value:" << setw(W) << engy << "\n";

                cout << "Stress:" << "\n" << setw(W) << str << "\n\n"
                        << "Stiffness:" << setw(W) << stiff
                        << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
                        << "Bifurcation Info:" << setw(W) << mintestfunct
                        << setw(W) << elastica_beam::NoNegTestFunctions << "\n";
            }
            break;
    }
}

ostream& operator<<(ostream& out, ElasticaBeam1D& A)
{
  A.Print(out, Lattice::PrintShort);
  return out;
}

void ElasticaBeam1D::PrintPath(ostream& out, ostream& pathout, int const& width)
{
    pathout << setw(width) << elastica_beam::NoNegTestFunctions << setw(width) << Lambda_ << " " << setw(width) << DOF_ << "\n";
}
