/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2015 by the deal.II authors and
 *                              & Jean-Paul Pelteret and Andrew McBride
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Jean-Paul Pelteret, University of Cape Town,
 *          Andrew McBride, University of Erlangen-Nuremberg, 2010
 *          Bastien Lauras, University of Minnesota, 2016
 */


// We start by including all the necessary deal.II header files and some C++
// related ones.
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>


// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace elastica_beam
{
  using namespace dealii;


  // @sect3{Run-time parameters}
  //
  // There are several parameters that can be set in the code so we set up a
  // ParameterHandler object to read in the choices at run-time.
  namespace Parameters
  {
    // @sect4{Finite Element system}

    struct FESystem
    {
      unsigned int poly_degree;
      unsigned int quad_order;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };


    void FESystem::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
	prm.declare_entry("Polynomial degree", "2",
			  Patterns::Integer(0),
			  "Displacement system polynomial order");

	prm.declare_entry("Quadrature order", "3",
			  Patterns::Integer(0),
			  "Gauss quadrature order");
      }
      prm.leave_subsection();
    }

    void FESystem::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
	poly_degree = prm.get_integer("Polynomial degree");
	quad_order = prm.get_integer("Quadrature order");
      }
      prm.leave_subsection();
    }

    // @sect4{Geometry}

    // Make adjustments to the problem geometry and the applied load.  Since the
    // problem modelled here is quite specific, the load scale can be altered to
    // specific values to compare with the results given in the literature.
    struct Geometry
    {
      unsigned int global_refinement;
      double       scale;
      double       p_p0;
      double       elongation;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Geometry::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
	prm.declare_entry("Global refinement", "0",
			  Patterns::Integer(0),
			  "Global refinement level");

	prm.declare_entry("Grid scale", "1e0",
			  Patterns::Double(0.0),
			  "Global grid scaling factor");

	prm.declare_entry("Pressure ratio p/p0", "100",
			  Patterns::Selection("20|40|60|80|100"),
			  "Ratio of applied pressure to reference pressure");

	prm.declare_entry("Elongation", "0.99",
			  Patterns::Double(0.0),
			  "Elongation");
      }
      prm.leave_subsection();
    }

    void Geometry::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
	global_refinement = prm.get_integer("Global refinement");
	scale = prm.get_double("Grid scale");
	p_p0 = prm.get_double("Pressure ratio p/p0");
	elongation = prm.get_double("Elongation");
      }
      prm.leave_subsection();
    }

    // @sect4{Materials}

    struct Materials
    {
      double EI;

      double k;

      double alpha;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Materials::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
	prm.declare_entry("Bending stiffness", "1",
			  Patterns::Double(),
			  "Bending stiffness");
        prm.declare_entry("Spring stiffness", "1",
			  Patterns::Double(),
			  "Spring stiffness");
        prm.declare_entry("Spring nonlinear stiffness", "0",
			  Patterns::Double(),
			  "Spring nonlinear stiffness");
      }
      prm.leave_subsection();
    }

    void Materials::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
	EI = prm.get_double("Bending stiffness");
        k = prm.get_double("Spring stiffness");
        alpha = prm.get_double("Spring nonlinear stiffness");
      }
      prm.leave_subsection();
    }

    // @sect4{Nonlinear solver}

    // A Newton-Raphson scheme is used to solve the nonlinear system of governing
    // equations.  We now define the tolerances and the maximum number of
    // iterations for the Newton-Raphson nonlinear solver.
    struct NonlinearSolver
    {
      unsigned int max_iterations_NR;
      double       tol_f;
      double       tol_u;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void NonlinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
	prm.declare_entry("Max iterations Newton-Raphson", "10",
			  Patterns::Integer(0),
			  "Number of Newton-Raphson iterations allowed");

	prm.declare_entry("Tolerance force", "1.0e-10",
			  Patterns::Double(0.0),
			  "Force residual tolerance");

	prm.declare_entry("Tolerance displacement", "1.0e-12",
			  Patterns::Double(0.0),
			  "Displacement error tolerance");
      }
      prm.leave_subsection();
    }

    void NonlinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
	max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
	tol_f = prm.get_double("Tolerance force");
	tol_u = prm.get_double("Tolerance displacement");
      }
      prm.leave_subsection();
    }

    // @sect4{Time}

    // Set the timestep size $ \varDelta t $ and the simulation end-time.
    struct Time
    {
      double delta_t;
      double end_time;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Time::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
	prm.declare_entry("End time", "1",
			  Patterns::Double(),
			  "End time");

	prm.declare_entry("Time step size", "0.1",
			  Patterns::Double(),
			  "Time step size");
      }
      prm.leave_subsection();
    }

    void Time::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
	end_time = prm.get_double("End time");
	delta_t = prm.get_double("Time step size");
      }
      prm.leave_subsection();
    }

    // @sect4{All parameters}

    // Finally we consolidate all of the above structures into a single container
    // that holds all of our run-time selections.
    struct AllParameters : public FESystem,
			   public Geometry,
			   public Materials,
			   public NonlinearSolver,
			   public Time

    {
      AllParameters(const std::string &input_file);

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    AllParameters::AllParameters(const std::string &input_file)
    {
      ParameterHandler prm;
      declare_parameters(prm);
      prm.read_input(input_file);
      parse_parameters(prm);
    }

    void AllParameters::declare_parameters(ParameterHandler &prm)
    {
      FESystem::declare_parameters(prm);
      Geometry::declare_parameters(prm);
      Materials::declare_parameters(prm);
      NonlinearSolver::declare_parameters(prm);
      Time::declare_parameters(prm);
    }

    void AllParameters::parse_parameters(ParameterHandler &prm)
    {
      FESystem::parse_parameters(prm);
      Geometry::parse_parameters(prm);
      Materials::parse_parameters(prm);
      NonlinearSolver::parse_parameters(prm);
      Time::parse_parameters(prm);
    }
  }

  // @sect3{Time class}

  // A simple class to store time data. Its functioning is transparent so no
  // discussion is necessary. For simplicity we assume a constant time step
  // size.
  class Time
  {
  public:
    Time (const double time_end,
	  const double delta_t)
      :
      timestep(0),
      time_current(0.0),
      time_end(time_end),
      delta_t(delta_t)
    {}

    virtual ~Time()
    {}

    double current() const
    {
      return time_current;
    }
    double end() const
    {
      return time_end;
    }
    double get_delta_t() const
    {
      return delta_t;
    }
    unsigned int get_timestep() const
    {
      return timestep;
    }
    void increment()
    {
      time_current += delta_t;
      ++timestep;
    }

  private:
    unsigned int timestep;
    double       time_current;
    const double time_end;
    const double delta_t;
  };


  template <int dim>
  class Elastica_Beam_On_Spring_Fundation
  {
  public:
    Elastica_Beam_On_Spring_Fundation(const double bending_stiffness,
            const double spring_stiffness, const double spring_nonlinear_stiffness)
      :
      EI(bending_stiffness),
      k(spring_stiffness),
      alpha(spring_nonlinear_stiffness)
    {}

    ~Elastica_Beam_On_Spring_Fundation()
    {}

    // This function is useless but we keep it in case we want to change the
    // material property during the iterations
    void update_material_data()
    {}

    // The next few functions return various data that we choose to store with
    // the material:

    double get_EI() const
    {
      return EI;
    }

    double get_k() const
    {
      return k;
    }

    double get_alpha() const
    {
      return alpha;
    }

  protected:
    // Define constitutive model parameters:
    const double EI;
    const double k;
    const double alpha;

  };

  //___________________________________________________________________________
  // @sect3{Quadrature point history}

  // As seen in step-18, the <code> PointHistory </code> class offers a method
  // for storing data at the quadrature points.  Here each quadrature point
  // holds a pointer to a material description.  Thus, different material models
  // can be used in different regions of the domain.
  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
      :
      material(NULL),
      y(0.0),
      d_y(Tensor<1, dim>()),
      dd_y(Tensor<2, dim>())
    {}

    virtual ~PointHistory()
    {
      delete material;
      material = NULL;
    }

    // The first function is used to create a material object and to
    // initialize all tensors correctly: The second one updates the stored
    // values and stresses based on the current deformation measure
    // $\textrm{Grad}\mathbf{u}_{\textrm{n}}$ values.
    void setup_lqp (const Parameters::AllParameters &parameters,
                const double value_y, const Tensor<1, dim> value_d_y)
    {
      material = new Elastica_Beam_On_Spring_Fundation<dim>(parameters.EI,
                                            parameters.k,parameters.alpha);
      update_values(value_y, value_d_y, Tensor<2, dim>());
    }

    // To this end, we calculate the deformation gradient $\mathbf{F}$ from
    // the displacement gradient $\textrm{Grad}\ \mathbf{u}$, i.e.
    // $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad}\ \mathbf{u}$ and
    // then let the material model associated with this quadrature point
    // update itself. When computing the deformation gradient, we have to take
    // care with which data types we compare the sum $\mathbf{I} +
    // \textrm{Grad}\ \mathbf{u}$: Since $I$ has data type SymmetricTensor,
    // just writing <code>I + Grad_u_n</code> would convert the second
    // argument to a symmetric tensor, perform the sum, and then cast the
    // result to a Tensor (i.e., the type of a possibly nonsymmetric
    // tensor). However, since <code>Grad_u_n</code> is nonsymmetric in
    // general, the conversion to SymmetricTensor will fail. We can avoid this
    // back and forth by converting $I$ to Tensor first, and then performing
    // the addition as between nonsymmetric tensors:
    void update_values (const double value_y, const Tensor<1, dim> value_d_y,
                                        const Tensor<2, dim> value_dd_y)
    {
        y = value_y;
        d_y = value_d_y;
        dd_y = value_dd_y;
    }

    // We offer an interface to retrieve certain data.  Here are the kinematic
    // variables:

    double get_y() const
    {
      return y;
    }

    const double &get_d_y() const
    {
      return d_y[0];
    }

    const double &get_dd_y() const
    {
      return dd_y[0][0];
    }

    // ...and the kinetic variables.  These are used in the material and
    // global tangent matrix and residual assembly operations:
    double get_EI() const
    {
      return material->get_EI();
    }

    double get_k() const
    {
      return material->get_k();
    }

    double get_alpha() const
    {
      return material->get_alpha();
    }

    // In terms of member functions, this class stores for the quadrature
    // point it represents a copy of a material type in case different
    // materials are used in different regions of the domain, as well as the
    // inverse of the deformation gradient...
  private:
    Elastica_Beam_On_Spring_Fundation<dim> *material;
    double y;
    Tensor<1, dim> d_y;
    Tensor<2, dim> dd_y;
  };



  // @sect3{Quasi-static quasi-incompressible finite-strain solid}

  // The Solid class is the central class in that it represents the problem at
  // hand. It follows the usual scheme in that all it really has is a
  // constructor, destructor and a <code>run()</code> function that dispatches
  // all the work to private functions of this class:
  template <int dim>
  class Solid
  {
  public:
    Solid(const std::string &input_file);

    virtual
    ~Solid();

    void
    run();

    std::size_t
    get_system_size();

    unsigned int
    get_unconstrained_system_size();

    double
    get_energy();

    void
    set_solution(Vector<double> const &solution);

    Vector<double> const&
    get_solution()
    {
      return solution_n;
    }

    void
    get_rhs_and_tangent(Vector<double> const* &sys_rhs,
			SparseMatrix<double> const* &tm, unsigned int iter_value);

    void
    get_E1DLoad(Vector<double> const* &sys_E1DLoad);

    void
    get_constraints_matrix(ConstraintMatrix const* &constraints_matrix);

    void
    set_P(const double value_P);

  private:

    // In the private section of this class, we first forward declare a number
    // of objects that are used in parallelizing work using the WorkStream
    // object (see the @ref threads module for more information on this).
    //
    // We declare such structures for the computation of tangent (stiffness)
    // matrix, right hand side, and for updating quadrature points:
    struct PerTaskData_K;
    struct ScratchData_K;

    struct PerTaskData_RHS;
    struct ScratchData_RHS;

    struct PerTaskData_E1DLoad;
    struct ScratchData_E1DLoad;

    struct PerTaskData_UQPH;
    struct ScratchData_UQPH;

    struct PerTaskData_Energy;
    struct ScratchData_Energy;

    // We start the collection of member functions with one that builds the
    // grid:
    void
    make_grid();

    // Set up the finite element system to be solved:
    void
    system_setup();

    // Several functions to assemble the system and right hand side matrices
    // using multithreading. Each of them comes as a wrapper function, one
    // that is executed to do the work in the WorkStream model on one cell,
    // and one that copies the work done on this one cell into the global
    // object that represents it:
    void
    assemble_system_tangent();

    void
    assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
				     ScratchData_K &scratch,
				     PerTaskData_K &data);

    void
    copy_local_to_global_K(const PerTaskData_K &data);

    void
    assemble_system_rhs();

    void
    assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
				 ScratchData_RHS &scratch,
				 PerTaskData_RHS &data);

    void
    copy_local_to_global_rhs(const PerTaskData_RHS &data);

    void
    assemble_system_E1DLoad();

    void
    assemble_system_E1DLoad_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
				 ScratchData_E1DLoad &scratch,
				 PerTaskData_E1DLoad &data);

    void
    copy_local_to_global_E1DLoad(const PerTaskData_E1DLoad &data);

    void
    assemble_system_energy();

    void
    assemble_system_energy_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
				    ScratchData_Energy &scratch,
				    PerTaskData_Energy &data);

    void
    copy_local_to_global_energy(const PerTaskData_Energy &/*data*/)
    {}

    // Apply Dirichlet boundary conditions on the displacement field
    void
    make_constraints(const int &it_nr);

    // Create and update the quadrature points. Here, no data needs to be
    // copied into a global object, so the copy_local_to_global function is
    // empty:
    void
    setup_qph();

    void
    update_qph_incremental(const Vector<double> &solution_delta);

    void
    update_qph_incremental_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
				    ScratchData_UQPH &scratch,
				    PerTaskData_UQPH &data);

    void
    copy_local_to_global_UQPH(const PerTaskData_UQPH &/*data*/)
    {}

    // Solve for the displacement using a Newton-Raphson method. We break this
    // function into the nonlinear loop and the function that solves the
    // linearized Newton-Raphson step:
    void
    solve_nonlinear_timestep(Vector<double> &solution_delta);

    void
    solve_linear_system(Vector<double> &newton_update, const int &it_nr);

    // Solution retrieval as well as post-processing and writing data to file:
    Vector<double>
    get_total_solution(const Vector<double> &solution_delta) const;

    void
    output_results() const;

    // Finally, some member variables that describe the current state: A
    // collection of the parameters used to describe the problem setup...
    Parameters::AllParameters        parameters;

    // ...and description of the geometry on which the problem is solved:
    Triangulation<dim>               triangulation;

    // Also, keep track of the current time and the time spent evaluating
    // certain functions
    Time                             time;
    TimerOutput                      timer;

    // A storage object for quadrature point information.  See step-18 for
    // more on this:
    std::vector<PointHistory<dim> >  quadrature_point_history;

    // A description of the finite-element system including the displacement
    // polynomial degree, the degree-of-freedom handler, number of DoFs per
    // cell and the extractor objects used to retrieve information from the
    // solution vectors:
    const unsigned int               degree;
    const FESystem<dim>              fe;
    DoFHandler<dim>                  dof_handler_ref;
    const unsigned int               dofs_per_cell;
    unsigned int                     n_dofs;


    // Rules for Gauss-quadrature on the cells. The number of
    // quadrature points on cells is recorded.
    const QGauss<dim>                qf_cell;
    const unsigned int               n_q_points;

    // Objects that store the converged solution and right-hand side vectors,
    // as well as the tangent matrix. There is a ConstraintMatrix object used
    // to keep track of constraints.  We make use of a sparsity pattern
    // designed for a block system.
    ConstraintMatrix            constraints;
    SparsityPattern             sparsity_pattern;
    SparseMatrix<double>        tangent_matrix;
    Vector<double>              system_rhs;
    Vector<double>              system_E1DLoad;
    Vector<double>              solution_n;
    double                      system_energy;
    double                       P;

    //Some boolean to decide if we want to display and save the tangent matrix and RHS
    const bool                       print_tangent_matrix = false;
    const bool                       print_tangent_matrix_constrained = false;
    const bool                       print_RHS = false;
    const bool                       print_steps_computation = true;
    const bool                       apply_non_zero_displacement_beginning = false;
    //This is just to know the size of the tangent matrix "by hand", for debugging:

    // Then define a number of variables to store norms and update norms and
    // normalisation factors.
    struct Errors
    {
      Errors()
	:
	norm(1.0)
      {}

      void reset()
      {
	norm = 1.0;
      }
      void normalise(const Errors &rhs)
      {
	if (rhs.norm != 0.0)
	  norm /= rhs.norm;
      }

      double norm;
    };

    Errors error_residual, error_residual_0, error_residual_norm, error_update,
      error_update_0, error_update_norm;

    // Methods to calculate error measures
    void get_error_residual(Errors &error_residual);

    void get_error_update(const Vector<double> &newton_update,
			  Errors &error_update);

    // Print information to screen in a pleasing way...
    static
    void
    print_conv_header();

    void
    print_conv_footer();
  };

  // @sect3{Implementation of the <code>Solid</code> class}

  // @sect4{Public interface}

  // We initialise the Solid class using data extracted from the parameter file.
  template <int dim>
  Solid<dim>::Solid(const std::string &input_file)
    :
    parameters(input_file),
    triangulation(Triangulation<dim>::maximum_smoothing),
    time(parameters.end_time, parameters.delta_t),
    timer(std::cout,
	  TimerOutput::summary,
	  TimerOutput::wall_times),
    degree(parameters.poly_degree),
    fe(FE_Q<dim>(parameters.poly_degree), dim), //Bastien : adjust
    dof_handler_ref(triangulation),
    dofs_per_cell (fe.dofs_per_cell),
    n_dofs (dof_handler_ref.n_dofs()),
    qf_cell(parameters.quad_order),
    n_q_points (qf_cell.size()),
    P (0.0)
  {
    make_grid();
    system_setup();
    output_results();
    time.increment();
  }

  // The class destructor simply clears the data held by the DOFHandler
  template <int dim>
  Solid<dim>::~Solid()
  {
    dof_handler_ref.clear();
  }


  // In solving the quasi-static problem, the time becomes a loading parameter,
  // i.e. we increasing the loading linearly with time, making the two concepts
  // interchangeable. We choose to increment time linearly using a constant time
  // step size.

  template <int dim>
  void Solid<dim>::run()
  {
    // Declare the incremental solution update $\varDelta
    // \mathbf{\Xi}:= \{\varDelta \mathbf{u}}$ and
    // start the loop over the time domain.
    //
    // At the beginning, we reset the solution update for this time step...
    Vector<double> solution_delta(n_dofs);
    while (time.current() < time.end())
      {
        if(apply_non_zero_displacement_beginning)
            P = 100*time.get_timestep();
	solution_delta = 0.0;

	// ...solve the current time step and update total solution vector
	// $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} +
	// \varDelta \mathbf{\Xi}$...
	solve_nonlinear_timestep(solution_delta);
	solution_n += solution_delta;

	// ...and plot the results before moving on happily to the next time
	// step:
	output_results();
	time.increment();
      }
  }


  // @sect3{Private interface}

  // @sect4{Threading-building-blocks structures}

  // The first group of private member functions is related to parallization.
  // We use the Threading Building Blocks library (TBB) to perform as many
  // computationally intensive distributed tasks as possible. In particular, we
  // assemble the tangent matrix and right hand side vector, the static
  // condensation contributions, and update data stored at the quadrature points
  // using TBB. Our main tool for this is the WorkStream class (see the @ref
  // threads module for more information).

  // Firstly we deal with the tangent matrix assembly structures.  The
  // PerTaskData object stores local contributions.
  template <int dim>
  struct Solid<dim>::PerTaskData_K
  {
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_K(const unsigned int dofs_per_cell)
      :
      cell_matrix(dofs_per_cell, dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
    }
  };


  // On the other hand, the ScratchData object stores the larger objects such as
  // the shape-function values array (<code>Nx</code>) and a shape function
  // gradient and symmetric gradient vector which we will use during the
  // assembly.
  template <int dim>
  struct Solid<dim>::ScratchData_K
  {
    FEValues<dim> fe_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<Tensor<1, dim> > >          d_Nx;
    std::vector<std::vector<Tensor<2, dim> > >          dd_Nx;


    ScratchData_K(const FiniteElement<dim> &fe_cell,
		  const QGauss<dim> &qf_cell,
		  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      Nx(qf_cell.size(),
	 std::vector<double>(fe_cell.dofs_per_cell)),
      d_Nx(qf_cell.size(),
	      std::vector<Tensor<1, dim> >(fe_cell.dofs_per_cell)),
      dd_Nx(qf_cell.size(),
	      std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell))
    {}

    ScratchData_K(const ScratchData_K &k)
      :
      fe_values_ref(k.fe_values_ref.get_fe(),
		    k.fe_values_ref.get_quadrature(),
		    k.fe_values_ref.get_update_flags()),
      Nx(k.Nx),
      d_Nx(k.d_Nx),
      dd_Nx(k.dd_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	  Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
	  Assert( d_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( dd_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
	  for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
	    {
	      Nx[q_point][k] = 0.0;
	      d_Nx[q_point][k] = 0.0;
              dd_Nx[q_point][k] = 0.0;
	    }
	}
    }

  };

  // Next, the same approach is used for the right-hand side assembly.  The
  // PerTaskData object again stores local contributions and the ScratchData
  // object the shape function object and precomputed values vector:
  template <int dim>
  struct Solid<dim>::PerTaskData_RHS
  {
    Vector<double>            cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_RHS(const unsigned int dofs_per_cell)
      :
      cell_rhs(dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_rhs = 0.0;
    }
  };

  template <int dim>
  struct Solid<dim>::ScratchData_RHS
  {
    FEValues<dim> fe_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<Tensor<1, dim> > >          d_Nx;
    std::vector<std::vector<Tensor<2, dim> > >          dd_Nx;


    ScratchData_RHS(const FiniteElement<dim> &fe_cell,
		  const QGauss<dim> &qf_cell,
		  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      Nx(qf_cell.size(),
	 std::vector<double>(fe_cell.dofs_per_cell)),
      d_Nx(qf_cell.size(),
	      std::vector<Tensor<1, dim> >(fe_cell.dofs_per_cell)),
      dd_Nx(qf_cell.size(),
	      std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell))
    {}

    ScratchData_RHS(const ScratchData_RHS &rhs)
      :
      fe_values_ref(rhs.fe_values_ref.get_fe(),
		    rhs.fe_values_ref.get_quadrature(),
		    rhs.fe_values_ref.get_update_flags()),
      Nx(rhs.Nx),
      d_Nx(rhs.d_Nx),
      dd_Nx(rhs.dd_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	  Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
	  Assert( d_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( dd_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
	  for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
	    {
	      Nx[q_point][k] = 0.0;
	      d_Nx[q_point][k] = 0.0;
              dd_Nx[q_point][k] = 0.0;
	    }
	}
    }

  };

  // Next, the same approach is used for the E1DLoad assembly.  The
  // PerTaskData object again stores local contributions and the ScratchData
  // object the shape function object and precomputed values vector:
  template <int dim>
  struct Solid<dim>::PerTaskData_E1DLoad
  {
    Vector<double>            cell_E1DLoad;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_E1DLoad(const unsigned int dofs_per_cell)
      :
      cell_E1DLoad(dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_E1DLoad = 0.0;
    }
  };

  template <int dim>
  struct Solid<dim>::ScratchData_E1DLoad
  {
    FEValues<dim> fe_values_ref;

    std::vector<std::vector<Tensor<1, dim> > >          d_Nx;

    ScratchData_E1DLoad(const FiniteElement<dim> &fe_cell,
		  const QGauss<dim> &qf_cell,
		  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      d_Nx(qf_cell.size(),
	      std::vector<Tensor<1, dim> >(fe_cell.dofs_per_cell))
    {}

    ScratchData_E1DLoad(const ScratchData_E1DLoad &E1DLoad)
      :
      fe_values_ref(E1DLoad.fe_values_ref.get_fe(),
		    E1DLoad.fe_values_ref.get_quadrature(),
		    E1DLoad.fe_values_ref.get_update_flags()),
      d_Nx(E1DLoad.d_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points = d_Nx.size();
      const unsigned int n_dofs_per_cell = d_Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	  Assert( d_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
	  for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
	    {
	      d_Nx[q_point][k] = 0.0;
	    }
	}
    }

  };


  // And finally we define the structures to assist with updating the quadrature
  // point information. We do not need the
  // PerTaskData object (since there is nothing to store here) but must define
  // one nonetheless. Note that this is because for the operation that we have
  // here -- updating the data on quadrature points -- the operation is purely
  // local: the things we do on every cell get consumed on every cell, without
  // any global aggregation operation as is usually the case when using the
  // WorkStream class. The fact that we still have to define a per-task data
  // structure points to the fact that the WorkStream class may be ill-suited to
  // this operation (we could, in principle simply create a new task using
  // Threads::new_task for each cell) but there is not much harm done to doing
  // it this way anyway.
  // Furthermore, should there be different material models associated with a
  // quadrature point, requiring varying levels of computational expense, then
  // the method used here could be advantageous.
  template <int dim>
  struct Solid<dim>::PerTaskData_UQPH
  {
    void reset()
    {}
  };


  // The ScratchData object will be used to store an alias for the solution
  // vector so that we don't have to copy this large data structure. We then
  // define a number of vectors to extract the solution values and gradients at
  // the quadrature points.
  template <int dim>
  struct Solid<dim>::ScratchData_UQPH
  {
    const Vector<double>   &solution_total;

    std::vector<double>          solution_y_total;
    std::vector<Tensor<1, dim> > solution_d_y_total;
    std::vector<Tensor<2, dim> > solution_dd_y_total;

    FEValues<dim>                fe_values_ref;

    ScratchData_UQPH(const FiniteElement<dim> &fe_cell,
		     const QGauss<dim> &qf_cell,
		     const UpdateFlags uf_cell,
		     const Vector<double> &solution_total)
      :
      solution_total(solution_total),
      solution_y_total(qf_cell.size()),
      solution_d_y_total(qf_cell.size()),
      solution_dd_y_total(qf_cell.size()),
      fe_values_ref(fe_cell, qf_cell, uf_cell)
    {}

    ScratchData_UQPH(const ScratchData_UQPH &uqph)
      :
      solution_total(uqph.solution_total),
      solution_y_total(uqph.solution_y_total),
      solution_d_y_total(uqph.solution_d_y_total),
      solution_dd_y_total(uqph.solution_dd_y_total),
      fe_values_ref(uqph.fe_values_ref.get_fe(),
		    uqph.fe_values_ref.get_quadrature(),
		    uqph.fe_values_ref.get_update_flags())
    {}

    void reset()
    {
      const unsigned int n_q_points = solution_y_total.size();
      for (unsigned int q = 0; q < n_q_points; ++q)
	{
	  solution_y_total[q] = 0.0;
	  solution_d_y_total[q] = 0.0;
          solution_dd_y_total[q] = 0.0;
	}
    }
  };

  template <int dim>
  struct Solid<dim>::PerTaskData_Energy
  {
    void reset()
    {}
  };


  // The ScratchData object will be used to store an alias for the solution
  // vector so that we don't have to copy this large data structure. We then
  // define a number of vectors to extract the solution values and gradients at
  // the quadrature points.


   template <int dim>
  struct Solid<dim>::ScratchData_Energy
  {
    FEValues<dim> fe_values_ref;



    ScratchData_Energy(const FiniteElement<dim> &fe_cell,
		  const QGauss<dim> &qf_cell,
		  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell)
    {}

    ScratchData_Energy(const ScratchData_Energy &energy)
      :
      fe_values_ref(energy.fe_values_ref.get_fe(),
		    energy.fe_values_ref.get_quadrature(),
		    energy.fe_values_ref.get_update_flags())
    {}

    void reset()
    {}

  };


  // @sect4{Solid::make_grid}

  // On to the first of the private member functions. Here we create the
  // triangulation of the domain, for which we choose the scaled cube with each
  // face given a boundary ID number.  The grid must be refined at least once
  // for the indentation problem.
  template <int dim>
  void Solid<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    //GridTools::scale(parameters.scale, triangulation);

    // We refine our mesh globally, at least once
    triangulation.refine_global(std::max (1U, parameters.global_refinement));

    const double length = GridTools::volume(triangulation);
    std::cout << "Grid:\n\t Length: " << length << std::endl;

    // We mark the surfaces in order to apply the boundary conditions after
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation.begin_active(), endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
	if (cell->face(f)->at_boundary())
	  {
	    const Point<dim> face_center = cell->face(f)->center();
            if (face_center[0] == 0) //Left
	      cell->face(f)->set_boundary_id (1);
	    else if (face_center[0] == 1) //Right
	      cell->face(f)->set_boundary_id (2);
	    else                          //No way
	      cell->face(f)->set_boundary_id (3);
            //Bastien : is this working in 1D ?
	  }
  }


  // @sect4{Solid::system_setup}

  // Next we describe how the FE system is setup.  We first determine the number
  // of components per block. Since the displacement is a vector component, the
  // first dim components belong to it, while the next one describe scalar
  // pressure.
  template <int dim>
  void Solid<dim>::system_setup()
  {
    timer.enter_subsection("Setup system");

    // The DOF handler is then initialised and we renumber the grid in an
    // efficient manner. We also record the number of DOFs per block.
    dof_handler_ref.distribute_dofs(fe);
    n_dofs = dof_handler_ref.n_dofs();
    //DoFRenumbering::Cuthill_McKee(dof_handler_ref);
    //DoFRenumbering::hierarchical(dof_handler_ref);
    //Point<dim> 	direction(0);
    //DoFRenumbering::downstream	(dof_handler_ref,direction,false);

      Triangulation<1>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
      std::cout << "\nPrinting coordinates of vertices...\n";
      for (; cell!=endc; ++cell)
        for (unsigned int v=0;
             v < GeometryInfo<1>::vertices_per_cell;
             ++v)
          {
            std::cout << cell->vertex(v) << "  ";
          }
      std::cout << std::endl;

    // Setup the sparsity pattern and tangent matrix
    tangent_matrix.clear();
    {
      //We print the number of cells and dofs:
      std::cout << std::endl;
      std::cout << "Triangulation:"
		<< "\n\t Number of active cells: " << triangulation.n_active_cells()
		<< "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
		<< "\n\t Number of dofs per cell: " << fe.dofs_per_cell
		<< std::endl;


      DynamicSparsityPattern dsp(dof_handler_ref.n_dofs());

      // This block creates the periodic constraints on our beam. The function
      // make_constraints() could maybe be deleted.
      const FEValuesExtractors::Scalar y_displacement(0);
      IndexSet selected_dofs;
      std::set< types::boundary_id > boundary_ids= std::set<types::boundary_id>();
      boundary_ids.insert(1);
      boundary_ids.insert(2);
      DoFTools::extract_boundary_dofs(dof_handler_ref,
                fe.component_mask(y_displacement), selected_dofs, boundary_ids);
      ConstraintMatrix::size_type dof_left;
      ConstraintMatrix::size_type dof_right;
      IndexSet::ElementIterator el = selected_dofs.begin();
      dof_left = *el;
      el++;
      dof_right = *el;
      constraints.add_line(dof_left);
      constraints.add_entry(dof_left,dof_right,-1.0);
      constraints.close();

      DoFTools::make_sparsity_pattern(dof_handler_ref,
				      dsp,
				      constraints,
				      true);
      sparsity_pattern.copy_from(dsp);
    }

    tangent_matrix.reinit(sparsity_pattern);

    // We then set up storage vectors
    system_rhs.reinit(dof_handler_ref.n_dofs());
    system_E1DLoad.reinit(dof_handler_ref.n_dofs());

    solution_n.reinit(dof_handler_ref.n_dofs());

    // ...and finally set up the quadrature point history:
    setup_qph();

    if(apply_non_zero_displacement_beginning)
        for(unsigned int i = 0; i < n_dofs; i++)
        {
            solution_n(i) = i * (n_dofs - i - 1) * 1.0 / (n_dofs * n_dofs);
        }

    timer.leave_subsection();
  }

  // @sect4{Solid::setup_qph}
  // The method used to store quadrature information is already described in
  // step-18. Here we implement a similar setup for a SMP machine.
  //
  // Firstly the actual QPH data objects are created. This must be done only
  // once the grid is refined to its finest level.
  template <int dim>
  void Solid<dim>::setup_qph()
  {
    if(print_steps_computation)
        std::cout << "    Setting up quadrature point data..." << std::endl;

    {
      triangulation.clear_user_data();
      {
	std::vector<PointHistory<dim> > tmp;
	tmp.swap(quadrature_point_history);
      }

      quadrature_point_history
	.resize(triangulation.n_active_cells() * n_q_points);

      unsigned int history_index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
	     triangulation.begin_active(); cell != triangulation.end();
	   ++cell)
	{
	  cell->set_user_pointer(&quadrature_point_history[history_index]);
	  history_index += n_q_points;
	}

      Assert(history_index == quadrature_point_history.size(),
	     ExcInternalError());
    }

    // Usefull only if we plan to apply a non-zero displacement at the beginning,
    // or use a graded displacement for instance. This allow us to obtain the
    // coordinates of the quadrature points in the real space.
    const UpdateFlags uf_cell(update_quadrature_points);
    FEValues<dim> fe_values(fe, qf_cell, uf_cell);

    // Next we setup the initial quadrature
    // point data:
    for (typename Triangulation<dim>::active_cell_iterator cell =
	   triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        fe_values.reinit(cell);
	PointHistory<dim> *lqph =
	  reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
	Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          double value_y = 0;
          Tensor<1, dim> value_d_y;
          if(apply_non_zero_displacement_beginning)
          {
            const Point<dim> local_point = fe_values.quadrature_point(q_point);
            //value_y = local_point[0] * (1 - local_point[0]) * 1.0;
            value_y = local_point[0] < 0.5 ? local_point[0] : 1 - local_point[0];
            //value_d_y[0] = 1 - 2 * local_point[0];
            value_d_y[0] = local_point[0] < 0.5 ? 1 : -1;
          }
	  lqph[q_point].setup_lqp(parameters, value_y, value_d_y);
        }
      }
  }

  // @sect4{Solid::update_qph_incremental}
  // As the update of QP information occurs frequently and involves a number of
  // expensive operations, we define a multithreaded approach to distributing
  // the task across a number of CPU cores.
  //
  // To start this, we first need to obtain the total solution as it stands at
  // this Newton increment and then create the initial copy of the scratch and
  // copy data objects:
  template <int dim>
  void Solid<dim>::update_qph_incremental(const Vector<double> &solution_delta)
  {
    timer.enter_subsection("Update QPH data");
    if(print_steps_computation)
        std::cout << " UQPH " << std::flush;

    const Vector<double> solution_total(get_total_solution(solution_delta));

    const UpdateFlags uf_UQPH(update_values | update_gradients | update_hessians);
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);

    // We then pass them and the one-cell update function to the WorkStream to
    // be processed:
    WorkStream::run(dof_handler_ref.begin_active(),
		    dof_handler_ref.end(),
		    *this,
		    &Solid::update_qph_incremental_one_cell,
		    &Solid::copy_local_to_global_UQPH,
		    scratch_data_UQPH,
		    per_task_data_UQPH);

    timer.leave_subsection();
  }


  // Now we describe how we extract data from the solution vector and pass it
  // along to each QP storage object for processing.
  template <int dim>
  void
  Solid<dim>::update_qph_incremental_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
					      ScratchData_UQPH &scratch,
					      PerTaskData_UQPH &/*data*/)
  {
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

    Assert(scratch.solution_y_total.size() == n_q_points,
	   ExcInternalError());
    Assert(scratch.solution_d_y_total.size() == n_q_points,
	   ExcInternalError());
    Assert(scratch.solution_dd_y_total.size() == n_q_points,
	   ExcInternalError());

    scratch.reset();

    // We first need to find the values and gradients at quadrature points
    // inside the current cell and then we update each local QP using the
    // displacement gradient and total pressure and dilatation solution
    // values:
    scratch.fe_values_ref.reinit(cell);

    scratch.fe_values_ref.get_function_values(scratch.solution_total,
						    scratch.solution_y_total);
    scratch.fe_values_ref.get_function_gradients(scratch.solution_total,
						    scratch.solution_d_y_total);
    scratch.fe_values_ref.get_function_hessians(scratch.solution_total,
						    scratch.solution_dd_y_total);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      lqph[q_point].update_values(scratch.solution_y_total[q_point],
                                        scratch.solution_d_y_total[q_point],
                                        scratch.solution_dd_y_total[q_point]);
  }


  // @sect4{Solid::solve_nonlinear_timestep}

  // The next function is the driver method for the Newton-Raphson scheme. At
  // its top we create a new vector to store the current Newton update step,
  // reset the error storage objects and print solver header.
  template <int dim>
  void
  Solid<dim>::solve_nonlinear_timestep(Vector<double> &solution_delta)
  {
    std::cout << std::endl << "Timestep " << time.get_timestep() << " @ "
	      << time.current() << "s" << std::endl;

    Vector<double> newton_update(n_dofs);

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();

    print_conv_header();

    // We now perform a number of Newton iterations to iteratively solve the
    // nonlinear problem.  Since the problem is fully nonlinear and we are
    // using a full Newton method, the data stored in the tangent matrix and
    // right-hand side vector is not reusable and must be cleared at each
    // Newton step.  We then initially build the right-hand side vector to
    // check for convergence (and store this value in the first iteration).
    // The unconstrained DOFs of the rhs vector hold the out-of-balance
    // forces. The building is done before assembling the system matrix as the
    // latter is an expensive operation and we can potentially avoid an extra
    // assembly process by not assembling the tangent matrix when convergence
    // is attained.
    unsigned int newton_iteration = 0;
    for (; newton_iteration < parameters.max_iterations_NR;
	 ++newton_iteration)
      {
	std::cout << " " << std::setw(2) << newton_iteration << " " << std::flush;


	tangent_matrix = 0.0;
	system_rhs = 0.0;

	// assemble rhs
	assemble_system_rhs();
	// assemble the tangent, make and impose the Dirichlet constraints,
	// and do the solve of the linearised system:
	assemble_system_tangent();
	make_constraints(newton_iteration);
	constraints.condense(tangent_matrix, system_rhs);

	get_error_residual(error_residual);

	if (newton_iteration <= 1)
	  error_residual_0 = error_residual;

	// We can now determine the normalised residual error and check for
	// solution convergence:
	error_residual_norm = error_residual;
	error_residual_norm.normalise(error_residual_0);

	if (newton_iteration > 0 && error_update_norm.norm <= parameters.tol_u
	    && error_residual_norm.norm <= parameters.tol_f)
	  {
	    std::cout << " CONVERGED! " << std::endl;
	    print_conv_footer();
	    break;
	  }

	solve_linear_system(newton_update,newton_iteration);

	get_error_update(newton_update, error_update);
	if (newton_iteration == 0)
	  error_update_0 = error_update;

	// We can now determine the normalised Newton update error, and
	// perform the actual update of the solution increment for the current
	// time step, update all quadrature point information pertaining to
	// this new displacement and stress state and continue iterating:
	error_update_norm = error_update;
	error_update_norm.normalise(error_update_0);

	solution_delta += newton_update;
	update_qph_incremental(solution_delta);

	std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
		  << std::scientific <<  error_residual_norm.norm
		  << "  " << error_update_norm.norm << std::endl;
      }

    // At the end, if it turns out that we have in fact done more iterations
    // than the parameter file allowed, we raise an exception that can be
    // caught in the main() function. The call <code>AssertThrow(condition,
    // exc_object)</code> is in essence equivalent to <code>if (!cond) throw
    // exc_object;</code> but the former form fills certain fields in the
    // exception object that identify the location (filename and line number)
    // where the exception was raised to make it simpler to identify where the
    // problem happened.

    AssertThrow (newton_iteration < parameters.max_iterations_NR,
		 ExcMessage("No convergence in nonlinear solver!"));
  }


  // @sect4{Solid::print_conv_header and Solid::print_conv_footer}

  // This program prints out data in a nice table that is updated
  // on a per-iteration basis. The next two functions set up the table
  // header and footer:
  template <int dim>
  void Solid<dim>::print_conv_header()
  {
    static const unsigned int l_width = 80;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "          SOLVER STEP            "
	      << " | RES_NORM     NU_NORM" << std::endl;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
  }



  template <int dim>
  void Solid<dim>::print_conv_footer()
  {
    static const unsigned int l_width = 80;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "Relative errors:" << std::endl
	      << "Displacement:\t" << error_update.norm / error_update_0.norm << std::endl
	      << "Force: \t\t" << error_residual.norm / error_residual_0.norm << std::endl;
  }


  // @sect4{Solid::get_error_residual}

  // Determine the true residual error for the problem.  That is, determine the
  // error in the residual for the unconstrained degrees of freedom.  Note that to
  // do so, we need to ignore constrained DOFs by setting the residual in these
  // vector components to zero.
  template <int dim>
  void Solid<dim>::get_error_residual(Errors &error_residual)
  {
    Vector<double> error_res(n_dofs);
    for (unsigned int i = 0; i < n_dofs; ++i)
      error_res(i) = (constraints.is_constrained(i)) ? 0.0 : system_rhs(i);
    error_residual.norm = error_res.l2_norm();
  }


  // @sect4{Solid::get_error_udpate}

  // Determine the true Newton update error for the problem
  template <int dim>
  void Solid<dim>::get_error_update(const Vector<double> &newton_update,
				    Errors &error_update)
  {
    Vector<double> error_ud(n_dofs);
    for (unsigned int i = 0; i < n_dofs; ++i)
      if (!constraints.is_constrained(i))
	error_ud(i) = newton_update(i);

    error_update.norm = error_ud.l2_norm();
  }



  // @sect4{Solid::get_total_solution}

  // This function provides the total solution, which is valid at any Newton step.
  // This is required as, to reduce computational error, the total solution is
  // only updated at the end of the timestep.
  template <int dim>
  Vector<double>
  Solid<dim>::get_total_solution(const Vector<double> &solution_delta) const
  {
    Vector<double> solution_total(solution_n);
    solution_total += solution_delta;
    return solution_total;
  }


  // @sect4{Solid::assemble_system_tangent}

  // Since we use TBB for assembly, we simply setup a copy of the
  // data structures required for the process and pass them, along
  // with the memory addresses of the assembly functions to the
  // WorkStream object for processing. Note that we must ensure that
  // the matrix is reset before any assembly operations can occur.
  template <int dim>
  void Solid<dim>::assemble_system_tangent()
  {
    timer.enter_subsection("Assemble tangent matrix");
    if(print_steps_computation)
        std::cout << " ASM_K" << std::flush;

    const UpdateFlags uf_cell(update_values    |  update_gradients |
			      update_hessians  |  update_JxW_values);

    PerTaskData_K per_task_data(dofs_per_cell);
    ScratchData_K scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
		    dof_handler_ref.end(),
		    *this,
		    &Solid::assemble_system_tangent_one_cell,
		    &Solid::copy_local_to_global_K,
		    scratch_data,
		    per_task_data);

    timer.leave_subsection();
  }

  // This function adds the local contribution to the system matrix.
  // Note that we choose not to use the constraint matrix to do the
  // job for us because the tangent matrix and residual processes have
  // been split up into two separate functions.
  template <int dim>
  void Solid<dim>::copy_local_to_global_K(const PerTaskData_K &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
	tangent_matrix.add(data.local_dof_indices[i],
			   data.local_dof_indices[j],
			   data.cell_matrix(i, j));
  }

  // Of course, we still have to define how we assemble the tangent matrix
  // contribution for a single cell.  We first need to reset and initialise some
  // of the scratch data structures and retrieve some basic information
  // regarding the DOF numbering on this cell.  We can precalculate the cell
  // shape function values and gradients. Note that the shape function gradients
  // are defined with regard to the intial configuration.
  template <int dim>
  void
  Solid<dim>::assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
					       ScratchData_K &scratch,
					       PerTaskData_K &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	for (unsigned int k = 0; k < dofs_per_cell; ++k)
	  {
	    scratch.Nx[q_point][k] = scratch.fe_values_ref.shape_value_component(k, q_point, 0);
            scratch.d_Nx[q_point][k] = scratch.fe_values_ref.shape_grad(k, q_point);
            scratch.dd_Nx[q_point][k] = scratch.fe_values_ref.shape_hessian(k, q_point);
	  }
      }

    // Now we build the local cell stiffness matrix. Since the global and
    // local system matrices are symmetric, we can exploit this property by
    // building only the lower half of the local matrix and copying the values
    // to the upper half.  So we only assemble half of the
    // $\mathsf{\mathbf{k}}_{uu}$ block, while complementary parts of the
    // $\mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{p}}$ and
    // $\mathsf{\mathbf{k}}_{\mathbf{p} \widetilde{u}}$ blocks are built.
    //
    // In doing so, we first extract some configuration dependent variable
    // from our QPH history objects for the current quadrature point.

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	const double EI                      = lqph[q_point].get_EI();
        const double k                       = lqph[q_point].get_k();
        const double alpha                   = lqph[q_point].get_alpha();
        const double y                       = lqph[q_point].get_y();
        const double d_y                     = lqph[q_point].get_d_y();
	const double dd_y                    = lqph[q_point].get_dd_y();

	const std::vector<double>
	  &vec_Nx = scratch.Nx[q_point];
	const std::vector<Tensor<1, dim> >
	  &vec_d_Nx = scratch.d_Nx[q_point];
        const std::vector<Tensor<2, dim> >
	  &vec_dd_Nx = scratch.dd_Nx[q_point];
	const double JxW = scratch.fe_values_ref.JxW(q_point);

        const double par = 1 - d_y * d_y;
	// We first compute the contributions
	// from the internal forces.  Note, by
	// definition of the rhs as the negative
	// of the residual, these contributions
	// are subtracted.

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	  {
            const double Nx_i = vec_Nx[i];
            const double d_Nx_i = vec_d_Nx[i][0];
            const double dd_Nx_i = vec_dd_Nx[i][0][0];
	    for (unsigned int j = 0; j <= i; ++j)
	      {
                const double Nx_j = vec_Nx[j];
                const double d_Nx_j = vec_d_Nx[j][0];
                const double dd_Nx_j = vec_dd_Nx[j][0][0];
		data.cell_matrix(i, j) += EI * (dd_Nx_j * dd_Nx_i / par
                                            + 2 * dd_y * d_y * dd_Nx_i * d_Nx_j / (par * par)
                                            + 2 * dd_y * d_y * dd_Nx_j * d_Nx_i / (par * par)
                                            + dd_y * dd_y * d_Nx_j * d_Nx_i / (par * par)
                                            + 4 * dd_y * dd_y * d_y * d_y * d_Nx_i * d_Nx_j / (par * par * par)
                                                ) * JxW;
                data.cell_matrix(i, j) -= P * (d_Nx_j * d_Nx_i / sqrt(par)
                                            + d_y * d_y * d_Nx_i * d_Nx_j / (par * sqrt(par))) * JxW;
                data.cell_matrix(i, j) += (k * Nx_j * Nx_i + 2 * alpha * y * Nx_j * Nx_i) * JxW;
	      }
	  }
      }

    //Finally, we need to copy the lower half of the local matrix into the
    //upper half:
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
	data.cell_matrix(i, j) = data.cell_matrix(j, i);

  }

  // @sect4{Solid::assemble_system_rhs}
  // The assembly of the right-hand side process is similar to the
  // tangent matrix, so we will not describe it in too much detail.
  // Note that since we want, maybe, to describe a problem with
  // Neumann BCs, we will need the face normals and so must specify
  // this in the update flags.
  template <int dim>
  void Solid<dim>::assemble_system_rhs()
  {
    timer.enter_subsection("Assemble system right-hand side");
    if(print_steps_computation)
        std::cout << " ASM_R" << std::flush;

    const UpdateFlags uf_cell(update_values    |  update_gradients |
			      update_hessians  |  update_JxW_values);

    PerTaskData_RHS per_task_data(dofs_per_cell);
    ScratchData_RHS scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
		    dof_handler_ref.end(),
		    *this,
		    &Solid::assemble_system_rhs_one_cell,
		    &Solid::copy_local_to_global_rhs,
		    scratch_data,
		    per_task_data);

    timer.leave_subsection();
  }



  template <int dim>
  void Solid<dim>::copy_local_to_global_rhs(const PerTaskData_RHS &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      system_rhs(data.local_dof_indices[i]) += data.cell_rhs(i);
  }



  template <int dim>
  void
  Solid<dim>::assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
					   ScratchData_RHS &scratch,
					   PerTaskData_RHS &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	for (unsigned int k = 0; k < dofs_per_cell; ++k)
	  {
	    scratch.Nx[q_point][k] = scratch.fe_values_ref.shape_value_component(k, q_point, 0);
            scratch.d_Nx[q_point][k] = scratch.fe_values_ref.shape_grad(k, q_point);
            scratch.dd_Nx[q_point][k] = scratch.fe_values_ref.shape_hessian(k, q_point);
	  }
      }

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	const double EI                      = lqph[q_point].get_EI();
        const double k                       = lqph[q_point].get_k();
        const double alpha                   = lqph[q_point].get_alpha();
        const double y                       = lqph[q_point].get_y();
        const double d_y                     = lqph[q_point].get_d_y();
	const double dd_y                    = lqph[q_point].get_dd_y();

	const std::vector<double>
	  &vec_Nx = scratch.Nx[q_point];
	const std::vector<Tensor<1, dim> >
	  &vec_d_Nx = scratch.d_Nx[q_point];
        const std::vector<Tensor<2, dim> >
	  &vec_dd_Nx = scratch.dd_Nx[q_point];
	const double JxW = scratch.fe_values_ref.JxW(q_point);

        const double par = 1 - d_y * d_y;
	// We first compute the contributions
	// from the internal forces.  Note, by
	// definition of the rhs as the negative
	// of the residual, these contributions
	// are subtracted.

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	  {
            const double Nx = vec_Nx[i];
            const double d_Nx = vec_d_Nx[i][0];
            const double dd_Nx = vec_dd_Nx[i][0][0];
	    data.cell_rhs(i) -= ( EI * (dd_y * dd_Nx / par
                                + dd_y * dd_y * d_y * d_Nx /(par * par))
                                - P * (d_y * d_Nx / sqrt(par))
                                + k * y * Nx + alpha * y * y * Nx) * JxW;
	  }
      }
  }

    // @sect4{Solid::assemble_system_E1DLoad}
  // The assembly of the right-hand side process is similar to the
  // tangent matrix, so we will not describe it in too much detail.
  // Note that since we want, maybe, to describe a problem with
  // Neumann BCs, we will need the face normals and so must specify
  // this in the update flags.
  template <int dim>
  void Solid<dim>::assemble_system_E1DLoad()
  {
    timer.enter_subsection("Assemble system E1DLoad");
    if(print_steps_computation)
        std::cout << " ASM_E1D" << std::flush;

    const UpdateFlags uf_cell(update_gradients |  update_JxW_values);
    // We use the same thing as the RHS cause it's almost the same assembly process
    PerTaskData_E1DLoad per_task_data(dofs_per_cell);
    ScratchData_E1DLoad scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
		    dof_handler_ref.end(),
		    *this,
		    &Solid::assemble_system_E1DLoad_one_cell,
		    &Solid::copy_local_to_global_E1DLoad,
		    scratch_data,
		    per_task_data);

    timer.leave_subsection();
  }



  template <int dim>
  void Solid<dim>::copy_local_to_global_E1DLoad(const PerTaskData_E1DLoad &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      system_E1DLoad(data.local_dof_indices[i]) += data.cell_E1DLoad(i);
  }



  template <int dim>
  void
  Solid<dim>::assemble_system_E1DLoad_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
					   ScratchData_E1DLoad &scratch,
					   PerTaskData_E1DLoad &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	for (unsigned int k = 0; k < dofs_per_cell; ++k)
	  {
            scratch.d_Nx[q_point][k] = scratch.fe_values_ref.shape_grad(k, q_point);
	  }
      }

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const double d_y                     = lqph[q_point].get_d_y();

	const std::vector<Tensor<1, dim> >
	  &vec_d_Nx = scratch.d_Nx[q_point];
	const double JxW = scratch.fe_values_ref.JxW(q_point);

        const double par = 1 - d_y * d_y;
	// We first compute the contributions
	// from the internal forces.  Note, by
	// definition of the rhs as the negative
	// of the residual, these contributions
	// are subtracted.

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	  {
            const double d_Nx = vec_d_Nx[i][0];
	    data.cell_E1DLoad(i) += (d_y * d_Nx / sqrt(par)) * JxW;
	  }
      }
  }

  // @sect4{Solid::make_constraints}
  // The constraints for this problem are simple to describe.
  // However, since we are dealing with an iterative Newton method,
  // it should be noted that any displacement constraints should only
  // be specified at the zeroth iteration and subsequently no
  // additional contributions are to be made since the constraints
  // are already exactly satisfied.
  template <int dim>
  void Solid<dim>::make_constraints(const int &it_nr)
  {//This is a useless function since constraints are done at the beginning also.
    if(print_steps_computation)
        std::cout << " CST  " << std::flush;

    // Since the constraints are different at different Newton iterations, we
    // need to clear the constraints matrix and completely rebuild
    // it. However, after the first iteration, the constraints remain the same
    // and we can simply skip the rebuilding step if we do not clear it.
    if (it_nr > 1)
      return;

    constraints.clear();

    const FEValuesExtractors::Scalar y_displacement(0);
    IndexSet selected_dofs;
    std::set< types::boundary_id > 	boundary_ids= std::set<types::boundary_id>();
    boundary_ids.insert(1);
    boundary_ids.insert(2);
    DoFTools::extract_boundary_dofs(dof_handler_ref,
                fe.component_mask(y_displacement), selected_dofs, boundary_ids);
    ConstraintMatrix::size_type dof_left;
    ConstraintMatrix::size_type dof_right;
    IndexSet::ElementIterator el = selected_dofs.begin();
    dof_left = *el;
    el++;
    dof_right = *el;
    constraints.add_line(dof_left);
    constraints.add_entry(dof_left,dof_right,-1.0);

    constraints.close();
  }

  // @sect4{Solid::solve_linear_system}

  template <int dim>
  void Solid<dim>::solve_linear_system(Vector<double> &newton_update, const int &it_nr)
  {
    {

      timer.enter_subsection("Linear solver");
      if(print_steps_computation)
          std::cout << " SLV" << std::flush;


      if(print_tangent_matrix)
	{
	  std::cout << std::endl << "Print tangent matrix "
                    << "after applying boundary conditions..." << std::endl;

	  FullMatrix<double> tangent_matrix_f(n_dofs);
	  tangent_matrix_f.copy_from(tangent_matrix);

	  std::cout << std::endl;

	  for(int x=0;x<n_dofs;x++)
	    {
	      for(int y=0;y<n_dofs;y++)
		{
		  if(tangent_matrix_f(x,y) == 0)
		    std::cout << "   0     |" << std::flush;
		  else if(tangent_matrix_f(x,y) < 1e-10 && tangent_matrix_f(x,y) > -1e-10)
		    std::cout << "   0     |" << std::flush;
		  else
		    {
		      if(tangent_matrix_f(x,y) >= 0)
			std::cout << " " << std::flush;
		      std::cout << std::scientific << std::setprecision(1) << tangent_matrix_f(x,y) << " |" << std::flush;  // display the current element out of the array
		    }
		}
	      std::cout << "  ///// Row number " << x+1 << std::endl;
	    }
	}


      if(print_tangent_matrix_constrained)
	{

	  FullMatrix<double> tangent_matrix_f(n_dofs);
	  tangent_matrix_f.copy_from(tangent_matrix);

	  std::cout << std::endl << "Extracting the intersting part..." << std::endl << std::endl;

	  //We open a file to write the tangent matrix in it
	  std::string tangent_matrix_file_name = "Matrices_RHS/tangent_matrix_timestep="
	    +Utilities::int_to_string(time.get_timestep())+"_nr_it="+Utilities::int_to_string(it_nr)+".txt";
	  std::ofstream tangent_matrix_file(tangent_matrix_file_name, std::ios::out | std::ios::trunc);

	  if(tangent_matrix_file)  // If opening the file was done well
	    {
	      for(int x=0;x<n_dofs;x++)
		{
		  if (!constraints.is_constrained(x))
		    {
		      for(int y=0;y<n_dofs;y++)
			if (!constraints.is_constrained(y))
			  {
			    if(tangent_matrix_f(x,y) == 0/* || (tangent_matrix_f(x,y) < 1e-10 && tangent_matrix_f(x,y) > -1e-10)*/){
			      std::cout << "   0      " << std::flush;
			      tangent_matrix_file << "   0      ";
			    }
			    else
			      {
				if(tangent_matrix_f(x,y) >= 0){
				  std::cout << " " << std::flush;
				  tangent_matrix_file << " ";
				}
				std::cout << std::scientific << std::setprecision(1)
					  << tangent_matrix_f(x,y) << "  " << std::flush;
				tangent_matrix_file << std::scientific << std::setprecision(1)
						    << tangent_matrix_f(x,y) << "  ";
			      }
			  }
		      std::cout << "  ///// Row number " << x+1 << std::endl;
		      tangent_matrix_file << std::endl;
		    }
		}
	      tangent_matrix_file.close();  // We close the file
	    }
	  else  // If there was an error opening the file
	    std::cerr << "Error while opening file tangent_matrix!" << std::endl;

	}
      if(print_RHS)
	{
	  //We open a file to write the RHS in it
	  std::string RHS_file_name = "Matrices_RHS/RHS_timestep="
	    +Utilities::int_to_string(time.get_timestep())+"_it_nr="+Utilities::int_to_string(it_nr)+".txt";
	  std::ofstream RHS_file(RHS_file_name, std::ios::out | std::ios::trunc);

	  if(RHS_file)  // If opening the file was done well
	    {
	      std::cout << std::endl << "Print right hand side..." << std::endl << std::endl;

	      for(int x=0;x<n_dofs;x++)
		if (!constraints.is_constrained(x))
		  {
		    if(system_rhs(x) == 0/* || (system_rhs(x) < 1e-10 && system_rhs(x) > -1e-10)*/){
		      std::cout << "   0     |  row " << x << std::endl;
		      RHS_file << "   0" << (constraints.is_constrained(x) ? 1 : 0) << std::endl;
		    }
		    else
		      {
			if(system_rhs(x) >= 0){
			  std::cout << " " << std::flush;
			  RHS_file << " ";
			}
			std::cout << std::scientific << std::setprecision(1)
				  << system_rhs(x) << " |  row " << x << std::endl;
			RHS_file << std::scientific << std::setprecision(1)
				 << system_rhs(x) << std::endl;
		      }
		  }
	      RHS_file.close();  // We close the file
	    }
	  else  // If there was an error opening the file
	    std::cerr << "Error while opening file RHS!" << std::endl;
	}

        SparseDirectUMFPACK A_direct;
        A_direct.initialize(tangent_matrix);
        A_direct.vmult(newton_update, system_rhs);

      timer.leave_subsection();
    }

    // Now that we have the displacement update, distribute the constraints
    // back to the Newton update:
    constraints.distribute(newton_update);
  }



  // @sect4{Solid::output_results}
  // Here we present how the results are written to file to be viewed
  // using ParaView or Visit.
  template <int dim>
  void Solid<dim>::output_results() const
  {
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(dim,
				    DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name(dim, "displacement");

    data_out.attach_dof_handler(dof_handler_ref);
    data_out.add_data_vector(solution_n,
			     solution_name,
			     DataOut<dim>::type_dof_data,
			     data_component_interpretation);

    // Since we are dealing with a large deformation problem, it would be nice
    // to display the result on a displaced grid!  The MappingQEulerian class
    // linked with the DataOut class provides an interface through which this
    // can be achieved without physically moving the grid points in the
    // Triangulation object ourselves.  We first need to copy the solution to
    // a temporary vector and then create the Eulerian mapping. We also
    // specify the polynomial degree to the DataOut object in order to produce
    // a more refined output data set when higher order polynomials are used.
    Vector<double> soln(solution_n.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping(degree, dof_handler_ref, soln);
    data_out.build_patches(q_mapping, degree);

    std::string strPath = "Results";
    int systemCmd = 0;
    if ( access( strPath.c_str(), 0 ) == -1 )
      {
	systemCmd = system("mkdir \"Results\"");
      }
    if(systemCmd != 0)
      std::cerr << std::endl << "Warning : problem while trying to create "
		<< "Results directory, using mkdir command" << std::endl;


    std::ostringstream filename;
    filename << "Results/timestep-" << time.get_timestep() << ".vtk";

    std::ofstream output(filename.str().c_str());
    data_out.write_vtk(output);

  }

  template <int dim>
  std::size_t
  Solid<dim>::get_system_size()
  {
    return n_dofs;
  }

  template <int dim>
  unsigned int
  Solid<dim>::get_unconstrained_system_size()
  {
    make_constraints(0);
    unsigned int nb_unconstrained_dofs = 0;
    for (unsigned int i = 0; i < n_dofs; ++i)
      if (!constraints.is_constrained(i))
	nb_unconstrained_dofs++;
    return nb_unconstrained_dofs;
  }

  template <int dim>
  double
  Solid<dim>::get_energy()
  {
    assemble_system_energy();
    return system_energy;
  }

  template <int dim>
  void Solid<dim>::assemble_system_energy()
  {
    timer.enter_subsection("Assemble system energy");

    const UpdateFlags uf_cell(update_values    |  update_gradients |
			      update_hessians  |  update_JxW_values);

    PerTaskData_Energy per_task_data;
    ScratchData_Energy scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
		    dof_handler_ref.end(),
		    *this,
		    &Solid::assemble_system_energy_one_cell,
		    &Solid::copy_local_to_global_energy,
		    scratch_data,
		    per_task_data);

    timer.leave_subsection();
  }

  template <int dim>
  void
  Solid<dim>::assemble_system_energy_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
					      ScratchData_Energy &scratch,
					      PerTaskData_Energy &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    //cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

   for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	const double EI                      = lqph[q_point].get_EI();
        const double k                       = lqph[q_point].get_k();
        const double alpha                   = lqph[q_point].get_alpha();
        const double y                       = lqph[q_point].get_y();
        const double d_y                     = lqph[q_point].get_d_y();
	const double dd_y                    = lqph[q_point].get_dd_y();

        const double JxW                     = scratch.fe_values_ref.JxW(q_point);
        const double par                     = 1 - d_y * d_y;
	// We first compute the contributions
	// from the internal forces.  Note, by
	// definition of the rhs as the negative
	// of the residual, these contributions
	// are subtracted.

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	  {
	    system_energy += ( 0.5 * EI * (dd_y * dd_y / par )
                                - P * (1 - sqrt(par))
                                + 0.5 * k * y * y + alpha * y * y * y / 3) * JxW;
	  }
      }
  }

  template <int dim>
  void
  Solid<dim>::set_solution(Vector<double> const &solution)
  {
    Vector<double> solution_delta(n_dofs);
    solution_delta = 0.0;

    solution_n = solution;
    update_qph_incremental(solution_delta);
  }

  template <int dim>
  void
  Solid<dim>::get_rhs_and_tangent(Vector<double> const* &sys_rhs,
				  SparseMatrix<double> const* &tm, unsigned int iter_value)
  {
    tangent_matrix = 0.0;
    system_rhs = 0.0;
    assemble_system_rhs();
    assemble_system_tangent();
    make_constraints(iter_value);
    constraints.condense(tangent_matrix, system_rhs);

    sys_rhs = &system_rhs;
    tm = &tangent_matrix;
  }

  template <int dim>
  void
  Solid<dim>::get_E1DLoad(Vector<double> const* &sys_E1DLoad)
  {
    //tangent_matrix = 0.0;
    system_E1DLoad = 0.0;
    assemble_system_E1DLoad();
    //assemble_system_tangent();
    //make_constraints(0);
    //constraints.condense(tangent_matrix, system_rhs);

    sys_E1DLoad = &system_E1DLoad;
  }

  template <int dim>
  void
  Solid<dim>::get_constraints_matrix(ConstraintMatrix const* &constraints_matrix)
  {
    constraints_matrix = &constraints;
  }

  template <int dim>
  void
  Solid<dim>::set_P(const double value_P)
  {
    P = value_P;
  }


#ifdef CREATE_LIBRARY
  template class Elastica_Beam_On_Spring_Fundation<2>;
  template class PointHistory<1>;
  Solid<1>* MyElasticaBeam;

  void createObject()
  {
      MyElasticaBeam = new Solid<1>("parameters.prm");
  }

  void deleteObject()
  {
      delete MyElasticaBeam;
  }

  void
  set_solution(double const* const solution)
  {
    //We assume here that the constraint matrix has already been created before.
    //This should be wise since this function is called after having computed
    //the size of the unconstrained system, in which we assemble the constraints.
    ConstraintMatrix const* constraints_matrix;
    MyElasticaBeam->get_constraints_matrix(constraints_matrix);
    std::size_t size(MyElasticaBeam->get_system_size());
    Vector<double> sol(size);
    unsigned int i_unconstrained = 0;
    for (unsigned int i = 0; i < size; ++i)
    {
        if (!constraints_matrix->is_constrained(i)){
            sol[i] = solution[i_unconstrained++];
        }
    }
    MyElasticaBeam->set_solution(sol);
  }

  void
  get_rhs_and_tangent(double* const sys_rhs,
		      double* const tm, unsigned int iter_value)
  {
    Vector<double> const* rhs;
    SparseMatrix<double> const* tangent;
    MyElasticaBeam->get_rhs_and_tangent(rhs,tangent,iter_value);

    std::size_t size(MyElasticaBeam->get_system_size());
    for (unsigned i = 0; i < size; ++i)
      {
	sys_rhs[i] = (*rhs)[i];
      }
    memset(tm, 0, size*size*sizeof(double));
    for (SparseMatrix<double>::const_iterator itr = tangent->begin();
	 itr != tangent->end(); ++itr)
      {
	tm[size*(itr->row()) + itr->column()] = itr->value();
      }
  }

  void
  get_unconstrained_rhs_and_tangent(double* const sys_rhs,
		      double* const tm, unsigned int iter_value)
  {
    Vector<double> const* rhs;
    SparseMatrix<double> const* tangent;
    MyElasticaBeam->get_rhs_and_tangent(rhs,tangent,iter_value);
    ConstraintMatrix const* constraints_matrix;
    MyElasticaBeam->get_constraints_matrix(constraints_matrix);
    std::size_t size(MyElasticaBeam->get_system_size());
    unsigned int unconstrained_size = 0;

    Vector<int> indices_unconstrained(size);

    unsigned int i_unconstrained = 0;
    for (unsigned int i = 0; i < size; ++i)
        if (!constraints_matrix->is_constrained(i)){
          sys_rhs[i_unconstrained] = (*rhs)[i];
          indices_unconstrained[i] = i_unconstrained++;
        }
    unconstrained_size = i_unconstrained;

    memset(tm, 0, unconstrained_size*unconstrained_size*sizeof(double));
    for (SparseMatrix<double>::const_iterator itr = tangent->begin();
	 itr != tangent->end(); ++itr)
    {
        if (!constraints_matrix->is_constrained(itr->row()))
        {
            if (!constraints_matrix->is_constrained(itr->column()))
            {
                tm[unconstrained_size*(indices_unconstrained[itr->row()]) + indices_unconstrained[itr->column()]] = itr->value();
            }
        }
    }
  }

  void
  get_E1DLoad(double* const sys_E1DLoad)
  {
    Vector<double> const* E1DLoad;
    MyElasticaBeam->get_E1DLoad(E1DLoad);

    std::size_t size(MyElasticaBeam->get_system_size());
    for (unsigned i = 0; i < size; ++i)
      {
	sys_E1DLoad[i] = (*E1DLoad)[i];
      }
  }

  double get_energy()
  {
    return MyElasticaBeam->get_energy();
  }

  std::size_t get_system_size()
  {
    return MyElasticaBeam->get_system_size();
  }

  unsigned int get_unconstrained_system_size()
  {
    return MyElasticaBeam->get_unconstrained_system_size();
  }

  void set_P(const double value_P)
  {
    MyElasticaBeam->set_P(value_P);
  }

  void run()
  {
    MyElasticaBeam->run();
  }
#endif
}


#ifndef CREATE_LIBRARY
// @sect3{Main function}
// Lastly we provide the main driver function which appears
// no different to the other tutorials.
int main ()
{
  using namespace dealii;
  using namespace elastica_beam;

  try
    {
      Solid<1> solid_1d("parameters.prm");
      solid_1d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
                << std::endl << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return -5;
    }

  return 0;
}
#endif
