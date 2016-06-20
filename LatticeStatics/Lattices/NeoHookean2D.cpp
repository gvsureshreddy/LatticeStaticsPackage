#include "NeoHookean2D.h"
#include <iostream>
#include <fstream>

namespace neo_hookean
{
  void createObject();
  void deleteObject();
  void run();
  std::size_t get_system_size();
  unsigned int get_unconstrained_system_size();
  void set_solution(double const* const solution, double* const solution_deal);
  void set_lambda(double const lambda);
  void get_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_unconstrained_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_unconstrained_E1DLoad(double* const E1DLoad, unsigned int iter_value);
  void output_results_BFB(double const lambda);
  double get_energy();
  int NoNegTestFunctions;
}


NeoHookean2D::~NeoHookean2D()
{
  neo_hookean::deleteObject();
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
  neo_hookean::createObject();
  LoadParameter_ = Load;
  for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
      CallCount_[i] = 0;
    }
  Caching_ = 0;


  //system_size_ = neo_hookean::get_system_size();
  DOFS_D_ = neo_hookean::get_unconstrained_system_size();
//  std::cout << "NeoHookean2D size is " << system_size_ << std::endl;
//  std::cout << "NeoHookean2D unconstrained size is "
//          << DOFS_D_ << std::endl;
  DOFS_ = DOFS_D_ + 0;
  DOF_.Resize(DOFS_, 0.0);
  DOF_D_.Resize(DOFS_D_, 0.0);
  RHS_.Resize(DOFS_,0.0);
  E1DLoad_.Resize(DOFS_,0.0);
  Stiff_.Resize(DOFS_,DOFS_,0.0);
  neo_hookean::set_solution(&(DOF_[0]), &(DOF_D_[0]));
  neo_hookean::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),1);
  //std::cout << setw(20) << RHS_;
  //std::cout << std::endl << std::endl << "Tangent matrix :\n\n\n" << std::endl;
  //std::cout << setw(20) << Stiff_;
}

void NeoHookean2D::SetLambda(double const& lambda)
{
    Lambda_ = lambda;
    neo_hookean::set_lambda(lambda);
}

void NeoHookean2D::SetDOF(Vector const& dof)
{
    DOF_ = dof;
    neo_hookean::set_solution(&(DOF_[0]), &(DOF_D_[0]));
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
  neo_hookean::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),1);
  return RHS_;
}

Vector const& NeoHookean2D::E1DLoad() const
{
  neo_hookean::get_unconstrained_E1DLoad(&(E1DLoad_[0]),1);
  return E1DLoad_;
}

Matrix const& NeoHookean2D::E2() const
{
  neo_hookean::get_unconstrained_rhs_and_tangent(&(RHS_[0]),&(Stiff_[0][0]),1);
  return Stiff_;
}

void NeoHookean2D::ExtraTestFunctions(Vector& TF) const
{
  return;
}



void NeoHookean2D::Print(ostream& out, PrintDetail const& flag,
			 PrintPathSolutionType const& SolType)
{
    int W;
    neo_hookean::NoNegTestFunctions = 0;
    double engy;
    double mintestfunct;
    out << "\nDOFS_ = " << DOFS_  << "\n";
    Matrix stiff(DOFS_, DOFS_);
    Vector str(DOFS_);
    Vector TestFunctVals(NumTestFunctions());
    W = out.width();
    W = 16;
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
            ++neo_hookean::NoNegTestFunctions;
        }
        if (mintestfunct > TestFunctVals[i])
        {
            mintestfunct = TestFunctVals[i];
        }
    }

    neo_hookean::output_results_BFB(Lambda_);

    switch (flag)
    {
        case PrintLong:
            out << "\nNeoHookean:\n";

            if (Echo_)
            {
                cout << "\nNeoHookean:\n";
            }
            // passthrough to short
        case PrintShort:
            out << "\n__________________________________________\n\n"
                    << "Lambda: " << setw(W) << Lambda_ << "\n"
                    << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                    << "Potential Value:" << setw(W) << engy << "\n";

            out /*<< "Stress:" << "\n" << setw(W) << str << "\n\n" //To be deleted
                    << "Stiffness:" << setw(W) << stiff //To be deleted*/
                    << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
                    << "Bifurcation Info:" << setw(W) << mintestfunct
                    << setw(W) << neo_hookean::NoNegTestFunctions << "\n";
            out << setw(W) << Lambda_ << " " << setw(W) << DOF_ << "\n";
            // send to cout also
            if (Echo_)
            {
                cout << "__________________________________________\n\n"
                        << "Lambda: " << setw(W) << Lambda_ << "\n"
                        << "DOF's :" << "\n" << setw(W) << DOF_ << "\n"
                        << "Potential Value:" << setw(W) << engy << "\n";

                cout /*<< "Stress:" << "\n" << setw(W) << str << "\n\n"
                        << "Stiffness:" << setw(W) << stiff*/
                        << "Eigenvalue Info:" << "\n" << setw(W) << TestFunctVals << "\n"
                        << "Bifurcation Info:" << setw(W) << mintestfunct
                        << setw(W) << neo_hookean::NoNegTestFunctions << "\n";
            }
            break;
    }
}

ostream& operator<<(ostream& out, NeoHookean2D& A)
{
  A.Print(out, Lattice::PrintShort);
  return out;
}

void NeoHookean2D::PrintPath(ostream& out, ostream& pathout, int const& width)
{
    Vector TestFunctVals(NumTestFunctions());
    TestFunctions(TestFunctVals, LHS);
    pathout << setw(width) << neo_hookean::NoNegTestFunctions << setw(width) << Lambda_ << " " << setw(width) << TestFunctVals /*<< " " << setw(width) << DOF_ */<< "\n";
}
