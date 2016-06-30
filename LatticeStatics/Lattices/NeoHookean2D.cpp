#include "NeoHookean2D.h"
#include <iostream>
#include <fstream>

namespace neo_hookean
{
  void createObject();
  void deleteObject();
  void run();
  std::size_t get_system_size();
  std::size_t get_unconstrained_system_size();
  std::size_t get_system_size_with_periodic();
  void get_dofs_properties(double* const dofs_properties);
  void get_constraint_properties(double* const constraint_properties);
  void set_solution(double const* const solution);
  void set_lambda(double const lambda);
  void get_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_unconstrained_rhs_and_tangent(double* const rhs, double* const tm, unsigned int iter_value);
  void get_free_rhs(double* const rhs);
  void get_free_tangent(double* const tm);
  void get_unconstrained_E1DLoad(double* const E1DLoad, unsigned int iter_value);
  void output_results_BFB();
  void output_results_BFB(bool isBifurcationTangent, unsigned int numBifurcationPoint, unsigned int indexLocal);
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


  DOFS_D_ = neo_hookean::get_system_size();
  DOFS_   = neo_hookean::get_unconstrained_system_size() + 1;
//  std::cout << "NeoHookean2D size with periodic is "
//          << DOFS_D_ << std::endl;
  DOF_.Resize(DOFS_, 0.0);
  for(unsigned int i = DOFS_-4; i < DOFS_; ++i)
      DOF_[i] = 1.0;
  DOF_D_.Resize(DOFS_D_, 0.0);
  dofs_properties_.Resize(3*DOFS_D_, 0.0);
  neo_hookean::get_dofs_properties(&(dofs_properties_[0]));
//  std::cout << "\nDOFs properties :\n";
//  for(unsigned int i = 0; i < DOFS_D_; ++i)
//  {
//        switch((int) dofs_properties_[3*i]){
//            case 0 :
//                std::cout << "\nDOF " << i << " is horizontal, at position : (" << dofs_properties_[3*i+1] << "," << dofs_properties_[3*i+2] << ")";
//                break;
//            case 1 :
//                std::cout << "\nDOF " << i << " is vertical, at position : (" << dofs_properties_[3*i+1] << "," << dofs_properties_[3*i+2] << ")";
//                break;
//            case 2 :
//                std::cout << "\nDOF " << i << " is a pressure. The position shall be zero : (" << dofs_properties_[3*i+1] << "," << dofs_properties_[3*i+2] << ")";
//                break;
//            default :
//                break;
//        }
//  }
  constraint_properties_.Resize(DOFS_D_, 0.0);
  neo_hookean::get_constraint_properties(&(constraint_properties_[0]));
  //std::cout << "\nPrint constraint_properties :\n";
//  for(unsigned int i = 0; i < DOFS_D_; ++i)
//  {
//        switch((int) constraint_properties_[i]){
//            case -3 :
//                std::cout << "\nDOF " << i << " is free.";
//                break;
//            case -2 :
//                std::cout << "\nDOF " << i << " is the node at the bottom right : zero horizontal displacement.";
//                break;
//            case -1 :
//                std::cout << "\nDOF " << i << " is a node at the bottom : zero vertical displacement.";
//                break;
//            default :
//                std::cout << "\nDOF " << i << " is periodically constrained with the node : " << constraint_properties_[i];
//                break;
//        }
//  }
  dofs_vertical_ = 0;
  links_from_constrained_to_unconstrained_.Resize(DOFS_D_, 0.0);
  fill_links_from_constrained_to_unconstrained();
//  std::cout << "\nPrint links_from_constrained_to_unconstrained :\n";
//  for(unsigned int i = 0; i < DOFS_D_; ++i)
//  {
//        if(links_from_constrained_to_unconstrained_[i] < 0){
//                std::cout << "\nDOF " << i << " is periodically constrained with the absolute node : " << -links_from_constrained_to_unconstrained_[i]-1;
//        } else if(links_from_constrained_to_unconstrained_[i] > 0)
//                std::cout << "\nDOF " << i << " is free and its unconstrained number is : " << links_from_constrained_to_unconstrained_[i];
//        else
//            std::cout << "\nDOF " << i << " blocked to a zero displacement.";
//  }
  RHS_.Resize(DOFS_,0.0);
  RHS_D_.Resize(DOFS_D_,0.0);
  E1DLoad_.Resize(DOFS_,0.0);
  E1DLoad_D_.Resize(DOFS_D_,0.0);
  Stiff_.Resize(DOFS_,DOFS_,0.0);
  Stiff_D_.Resize(DOFS_D_,DOFS_D_,0.0);
  SetDOF(DOF_);
}

void NeoHookean2D::fill_links_from_constrained_to_unconstrained()
{
    unsigned int i_unconstrained = 1;
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        if(constraint_properties_[i] > -1)
            // Periodically constrained DOFs. We set the indice to the related node minus 2.
            links_from_constrained_to_unconstrained_[i] = -constraint_properties_[i]-1;
        else if(constraint_properties_[i] == -3)
        {
            links_from_constrained_to_unconstrained_[i] = i_unconstrained++; // Free DOFs
            if(dofs_properties_[3*i] == 1.0)
                dofs_vertical_++;
        }
        else
            links_from_constrained_to_unconstrained_[i] = 0; // Blocked to zero DOFs
    }
    std::cout << "\ndofs_vertical_ = " << dofs_vertical_ << "\n";
    dof_vertical_.Resize(dofs_vertical_, (unsigned int) 0);
    unsigned int indice = 0;
    for(unsigned int i = 0; i < DOFS_D_; ++i)
        if(constraint_properties_[i] == -3 && dofs_properties_[3*i] == 1.0)
            dof_vertical_[indice++] = i;;
    std::cout << "dof_vertical_ = " << dof_vertical_ << "\n";
}

void NeoHookean2D::SetLambda(double const& lambda)
{
    for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
    }
    Lambda_ = lambda;
    neo_hookean::set_lambda(lambda);
    SetDOF(DOF_);
}

void NeoHookean2D::SetDOF(Vector const& dof)
{
    for (int i = 0; i < cachesize; ++i)
    {
      Cached_[i] = 0;
    }
    DOF_ = dof;
    //std::cout << "\n\nWe have : lambda = " << Lambda_ << " and eta = " << DOF_[0] << "\n";
    for(int i = 0; i < DOFS_D_; ++i)
    {
        switch((int) constraint_properties_[i])
        {
            case -3 : // Free DOF
                switch((int) dofs_properties_[3*i])
                {
                    case 0 :
                        DOF_D_[i] = (1 - dofs_properties_[3*i+1]) * Lambda_ + DOF_[links_from_constrained_to_unconstrained_[i]];
                        //std::cout << "\nThis is a H free dof. Added : " << (1 - dofs_properties_[3*i+1]) * Lambda_ << " to the displacement : " << DOF_[links_from_constrained_to_unconstrained_[i]];
                        break;
                    case 1 :
                        DOF_D_[i] = dofs_properties_[3*i+2] * DOF_[0] + DOF_[links_from_constrained_to_unconstrained_[i]];
                        //std::cout << "\nThis is a V free dof. Added : " << dofs_properties_[3*i+2] * DOF_[0] << " to the displacement : " << DOF_[links_from_constrained_to_unconstrained_[i]];
                        break;
                    case 2 :
                        DOF_D_[i] = DOF_[links_from_constrained_to_unconstrained_[i]];
                        //std::cout << "\nThis is a P dof. Added : " << 0 << " to the displacement : " << DOF_[links_from_constrained_to_unconstrained_[i]];
                        break;
                    default :
                        std::cerr << "Internal problem in Set_DOF, type 1." << std::endl;
                        break;
                }
                break;
            case -2 : // Bottom right horizontal DOF
                DOF_D_[i] = 0.0;
                break;
            case -1 : // Bottom vertical DOFs
                DOF_D_[i] = 0.0;
                break;
            default : // Periodically constrained DOFs (positive integer)
                 if(constraint_properties_[constraint_properties_[i]] != -3)
                 { // Handle the bottom left DOF
                     DOF_D_[i] = Lambda_;
                     break;
                 }
                switch((int) dofs_properties_[3*i])
                {
                    case 0 :
                        DOF_D_[i] = (1 - dofs_properties_[3*i+1]) * Lambda_ + DOF_[links_from_constrained_to_unconstrained_[constraint_properties_[i]]];
                        //std::cout << "\nThis is a H free dof. Added : " << (1 - dofs_properties_[3*i+1]) * Lambda_ << " to the displacement : " << DOF_[links_from_constrained_to_unconstrained_[constraint_properties_[i]]];
                        break;
                    case 1 :
                        DOF_D_[i] = dofs_properties_[3*i+2] * DOF_[0] + DOF_[links_from_constrained_to_unconstrained_[constraint_properties_[i]]];
                        //std::cout << "\nThis is a V free dof. Added : " << dofs_properties_[3*i+2] * DOF_[0] << " to the displacement : " << DOF_[links_from_constrained_to_unconstrained_[constraint_properties_[i]]];
                        break;
                    default :
                        std::cerr << "Internal problem in Set_DOF, type 2." << std::endl;
                        break;
                }
                break;
        }
    }
    //std::cout << "\nSetting solution to :\n" << DOF_D_ << "\n\nDOF was equal to : " << DOF_;
    neo_hookean::set_solution(&(DOF_D_[0]));
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
  if ((!Caching_) || (!Cached_[1]))
  {
    RHS_.Resize(DOFS_,0.0);
    RHS_D_.Resize(DOFS_D_,0.0);
    neo_hookean::get_free_rhs(&(RHS_D_[0]));
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        const bool condition = (constraint_properties_[i] == -3 || constraint_properties_[i] > -1) && dofs_properties_[3*i] == 1.0;
        RHS_[0] -= condition ? dofs_properties_[3*i+2] * RHS_D_[i] : 0.0;
        if(constraint_properties_[i] == -3)
        {
            RHS_[links_from_constrained_to_unconstrained_[i]] -= RHS_D_[i];
        } else if(constraint_properties_[i] > -1 && constraint_properties_[constraint_properties_[i]] == -3)
        {
            RHS_[links_from_constrained_to_unconstrained_[constraint_properties_[i]]] -= RHS_D_[i];
        }
    }
    double sum_vertical_displacements = 0.0;
    for(unsigned int i = 0; i < dofs_vertical_; ++i)
    {
        sum_vertical_displacements += DOF_[links_from_constrained_to_unconstrained_[dof_vertical_[i]]];
    }
    for(unsigned int i = 0; i < dofs_vertical_; ++i)
    {
        RHS_[links_from_constrained_to_unconstrained_[dof_vertical_[i]]] += factor_penalty_ * sum_vertical_displacements;
    }
    Cached_[1] = 1;
    CallCount_[1]++;
  }
  //std::cout << "\n\nPrint RHS_D :\n\n" << setw(16) << RHS_D_;
  //std::cout << "\n\nPrint RHS :\n\n" << setw(16) << RHS_;
  return RHS_;
}

Vector const& NeoHookean2D::E1DLoad() const
{
  if ((!Caching_) || (!Cached_[2]))
  {
    E1DLoad_.Resize(DOFS_,0.0);
    Stiff_D_.Resize(DOFS_D_,DOFS_D_,0.0);
    neo_hookean::get_free_tangent(&(Stiff_D_[0][0]));
    Matrix temp;
    temp.Resize(DOFS_,DOFS_D_, 0.0);
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        const bool condition = (constraint_properties_[i] == -3 || constraint_properties_[i] > -1) && dofs_properties_[3*i] == 1.0;
        for(unsigned int j = 0; j < DOFS_D_; ++j)
        {
            temp[0][j] += condition ? dofs_properties_[3*i+2] * Stiff_D_[i][j] : 0.0;
            if(constraint_properties_[i] == -3)
            {
                temp[(int) links_from_constrained_to_unconstrained_[i]][j] += Stiff_D_[i][j];
            } else if(constraint_properties_[i] > -1 && constraint_properties_[constraint_properties_[i]] == -3)
            {
                temp[(int) links_from_constrained_to_unconstrained_[(int) constraint_properties_[i]]][j] += Stiff_D_[i][j];
            }
        }
    }
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        const bool condition = (constraint_properties_[i] == -3 || constraint_properties_[i] > -1) && dofs_properties_[3*i] == 0.0;
        if(condition)
            for(unsigned int j = 0; j < DOFS_; ++j)
            {
                E1DLoad_[j] += temp[j][i] * (1 - dofs_properties_[3*i+1]);
            }
    }
    Cached_[2] = 1;
    CallCount_[2]++;
  }
//  std::cout << "\n\nPrint Stiff_D :\n\n" << setw(17) << Stiff_D_;
//  std::cout << "\n\nPrint E1DLoad :\n\n" << setw(17) << E1DLoad_;
  return E1DLoad_;
}

Matrix const& NeoHookean2D::E2() const
{
  if ((!Caching_) || (!Cached_[3]))
  {
    Stiff_.Resize(DOFS_,DOFS_,0.0);
    Stiff_D_.Resize(DOFS_D_,DOFS_D_,0.0);
    neo_hookean::get_free_tangent(&(Stiff_D_[0][0]));
    Matrix temp;
    temp.Resize(DOFS_,DOFS_D_, 0.0);
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        const bool condition = (constraint_properties_[i] == -3 || constraint_properties_[i] > -1) && dofs_properties_[3*i] == 1.0;
        for(unsigned int j = 0; j < DOFS_D_; ++j)
        {
            temp[0][j] += condition ? dofs_properties_[3*i+2] * Stiff_D_[i][j] : 0.0;
            if(constraint_properties_[i] == -3)
            {
                temp[(int) links_from_constrained_to_unconstrained_[i]][j] += Stiff_D_[i][j];
            } else if(constraint_properties_[i] > -1 && constraint_properties_[constraint_properties_[i]] == -3)
            {
                temp[(int) links_from_constrained_to_unconstrained_[(int) constraint_properties_[i]]][j] += Stiff_D_[i][j];
            }
        }
    }
    for(unsigned int i = 0; i < DOFS_D_; ++i)
    {
        const bool condition = (constraint_properties_[i] == -3 || constraint_properties_[i] > -1) && dofs_properties_[3*i] == 1.0;
        Stiff_[0][0] += condition ? temp[0][i] * dofs_properties_[3*i+2] : 0.0;
        if(constraint_properties_[i] == -3)
        {
            Stiff_[0][(int) links_from_constrained_to_unconstrained_[i]] += temp[0][i];
            Stiff_[(int) links_from_constrained_to_unconstrained_[i]][0] += temp[0][i];
        } else if(constraint_properties_[i] > -1 && constraint_properties_[constraint_properties_[i]] == -3)
        {
            Stiff_[0][(int) links_from_constrained_to_unconstrained_[(int) constraint_properties_[i]]] += temp[0][i];
            Stiff_[(int) links_from_constrained_to_unconstrained_[(int) constraint_properties_[i]]][0] += temp[0][i];
        }
        for(unsigned int j = 1; j < DOFS_; ++j)
        {
            if(constraint_properties_[i] == -3)
            {
                Stiff_[j][(int) links_from_constrained_to_unconstrained_[i]] += temp[j][i];
            } else if(constraint_properties_[i] > -1 && constraint_properties_[constraint_properties_[i]] == -3)
            {
                Stiff_[j][(int) links_from_constrained_to_unconstrained_[(int) constraint_properties_[i]]] += temp[j][i];
            }
        }
    }
    for(unsigned int i = 0; i < dofs_vertical_; ++i)
    {
        Stiff_[(int) links_from_constrained_to_unconstrained_[dof_vertical_[i]]][(int) links_from_constrained_to_unconstrained_[dof_vertical_[i]]] += factor_penalty_;
    }
    Cached_[3] = 1;
    CallCount_[3]++;
  }
//  std::cout << "\n\nPrint Stiff_D :\n\n" << setw(17) << Stiff_D_;
//  std::cout << "\n\nPrint Stiff :\n\n" << setw(17) << Stiff_ << "\n\n";
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

    neo_hookean::output_results_BFB();

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

void NeoHookean2D::DrawBifurcatedPath(Vector const& tangent, unsigned int numBifurcationPoint, unsigned int indexLocal)
{
    Vector backupSolution(DOF_);
    std::cout << "\n\n\n_________________________________________________________________________\n\nOk, here is backupSolution"     << backupSolution << "\n\n";
    Vector tempSolution(DOF_);
    for(unsigned int i = 0; i < DOFS_; i++)
    {
        tempSolution[i] += tangent[i]*sqrt(DOFS_)/50;
        // This is just a factor to make it nice to display, and to have the
        // same average displacement with respect to the refinement level
    }
    SetDOF(tempSolution);
    neo_hookean::output_results_BFB(true, numBifurcationPoint, indexLocal);
    SetDOF(backupSolution);
}
