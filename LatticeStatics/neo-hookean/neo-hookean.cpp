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

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
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
namespace neo_hookean
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
                prm.declare_entry("Polynomial degree", "1",
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
            unsigned int local_refinement_cycles;
            double       scale;
            double       p_p0;
            double       displacement_side_1;
            
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
                
                prm.declare_entry("Local refinement cycles", "2",
                        Patterns::Integer(0),
                        "Local refinement cycles");
                
                prm.declare_entry("Grid scale", "1e0",
                        Patterns::Double(0.0),
                        "Global grid scaling factor");
                
                prm.declare_entry("Pressure ratio p/p0", "100",
                        Patterns::Selection("20|40|60|80|100"),
                        "Ratio of applied pressure to reference pressure");
                
                prm.declare_entry("Displacement side 1", "0.001",
                        Patterns::Double(0.0),
                        "Displacement side 1");
            }
            prm.leave_subsection();
        }
        
        void Geometry::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Geometry");
            {
                global_refinement = prm.get_integer("Global refinement");
                local_refinement_cycles = prm.get_integer("Local refinement cycles");
                scale = prm.get_double("Grid scale");
                p_p0 = prm.get_double("Pressure ratio p/p0");
                displacement_side_1 = prm.get_double("Displacement side 1");
            }
            prm.leave_subsection();
        }
        
        // @sect4{Materials}
        
        // We also need the shear modulus $ \mu $ for the incompressible
        // neo-Hookean material.
        struct Materials
        {
            double mu;
            
            static void
            declare_parameters(ParameterHandler &prm);
            
            void
            parse_parameters(ParameterHandler &prm);
        };
        
        void Materials::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Material properties");
            {
                prm.declare_entry("Shear modulus", "2"/*"80.194e6"*/,
                        Patterns::Double(),
                        "Shear modulus");
            }
            prm.leave_subsection();
        }
        
        void Materials::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Material properties");
            {
                mu = prm.get_double("Shear modulus");
            }
            prm.leave_subsection();
        }
        
        // @sect4{Linear solver}
        
        // Next, we choose both solver and preconditioner settings.  The use of an
        // effective preconditioner is critical to ensure convergence when a large
        // nonlinear motion occurs within a Newton increment.
        struct LinearSolver
        {
            std::string type_lin;
            double      tol_lin;
            double      max_iterations_lin;
            std::string preconditioner_type;
            double      preconditioner_relaxation;
            
            static void
            declare_parameters(ParameterHandler &prm);
            
            void
            parse_parameters(ParameterHandler &prm);
        };
        
        //Appart the first entry, "Linear solver", everything's useless
        //if we use a direct solver
        void LinearSolver::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Linear solver");
            {
                prm.declare_entry("Solver type", "Direct",
                        Patterns::Selection("CG|Direct"),
                        "Type of solver used to solve the linear system");
                
                prm.declare_entry("Residual", "1e-6",
                        Patterns::Double(0.0),
                        "Linear solver residual (scaled by residual norm)");
                
                prm.declare_entry("Max iteration multiplier", "1",
                        Patterns::Double(0.0),
                        "Linear solver iterations (multiples of the system matrix size)");
                
                prm.declare_entry("Preconditioner type", "ssor",
                        Patterns::Selection("jacobi|ssor"),
                        "Type of preconditioner");
                
                prm.declare_entry("Preconditioner relaxation", "0.65",
                        Patterns::Double(0.0),
                        "Preconditioner relaxation value");
            }
            prm.leave_subsection();
        }
        
        void LinearSolver::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Linear solver");
            {
                type_lin = prm.get("Solver type");
                tol_lin = prm.get_double("Residual");
                max_iterations_lin = prm.get_double("Max iteration multiplier");
                preconditioner_type = prm.get("Preconditioner type");
                preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
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
                public LinearSolver,
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
            LinearSolver::declare_parameters(prm);
            NonlinearSolver::declare_parameters(prm);
            Time::declare_parameters(prm);
        }
        
        void AllParameters::parse_parameters(ParameterHandler &prm)
        {
            FESystem::parse_parameters(prm);
            Geometry::parse_parameters(prm);
            Materials::parse_parameters(prm);
            LinearSolver::parse_parameters(prm);
            NonlinearSolver::parse_parameters(prm);
            Time::parse_parameters(prm);
        }
    }
    
    // @sect3{Some standard tensors}
    
    // Now we define one frequently used second-order tensor:
    template <int dim>
    class StandardTensors
    {
        public:
            // $\mathbf{I}$
            static const SymmetricTensor<2, dim> I;
    };
    
    template <int dim>
    const SymmetricTensor<2, dim>
    StandardTensors<dim>::I = unit_symmetric_tensor<dim>();
    
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
    
    // @sect3{Incompressible neo-Hookean material within a two-field formulation}
    
    // As discussed in the Introduction, Neo-Hookean materials are a type of
    // hyperelastic materials.  The entire domain is assumed to be composed of an
    // incompressible neo-Hookean material.  This class defines the behaviour of
    // this material within a two-field formulation.
    //
    // The following class will be used to characterize the material we work with,
    // and provides a central point that one would need to modify if one were to
    // implement a different material model. For it to work, we will store one
    // object of this type per quadrature point, and in each of these objects
    // store the current state (characterized by the values or measures  of the two fields)
    // so that we can compute the elastic coefficients linearized around the
    // current state.
    template <int dim>
    class Material_Incompressible_Neo_Hook_Two_Field
    {
    public:
        Material_Incompressible_Neo_Hook_Two_Field(const double mu)
        :
        c_1(mu / 2.0),
                det_F(1.0),
                p(mu / 2.0)
        {}
        
        ~Material_Incompressible_Neo_Hook_Two_Field()
        {}
        
        // We update the material model with various deformation dependent data
        // based on $F$ and the pressure p, and at the end of the function include a physical
        // check for internal consistency:
        void update_material_data(const Tensor<2, dim> &F,
        const double p_in)
        {
            det_F = determinant(F);
            p = p_in;
            
            Assert(det_F > 0, ExcInternalError());
        }
        
        // The next few functions return various data that we choose to store with
        // the material:
        double get_det_F() const
        {
            return det_F;
        }
        
        double get_p() const
        {
            return p;
        }
        
        double get_c_1() const
        {
            return c_1;
        }
        
    protected:
        // Define constitutive neo-Hookean model parameter $c_1$:
        const double c_1;
        
        // Model specific data that is convenient to store with the material:
        double det_F;
        double p;
        
    };
    
    //________________________________________________________________________________________________
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
                C(Tensor<2, dim>()),
                C_inv(Tensor<2, dim>()),
                Grad_U(Tensor<2, dim>())
        {}
        
        virtual ~PointHistory()
        {
            delete material;
            material = NULL;
        }
        
        // The first function is used to create a material object and to
        // initialize all tensors correctly: The second one updates the stored
        // values and stresses based on the current deformation measure
        // $\textrm{Grad}\mathbf{u}_{\textrm{n}}$ and pressure $\widetilde{p}$ field values.
        void setup_lqp (const Parameters::AllParameters &parameters)
        {
            material = new Material_Incompressible_Neo_Hook_Two_Field<dim>(parameters.mu);
            update_values(Tensor<2, dim>(), parameters.mu / 2.0);
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
        void update_values (const Tensor<2, dim> &Grad_u_n,
        const double p)
        {
            const Tensor<2, dim> F
            = (Tensor<2, dim>(StandardTensors<dim>::I) +
            Grad_u_n);
            material->update_material_data(F, p);
            C = transpose(F) * F;
            C_inv = invert(C);
            Grad_U = Grad_u_n;
            
        }
        
        // We offer an interface to retrieve certain data.  Here are the kinematic
        // variables:
        
        double get_det_F() const
        {
            return material->get_det_F();
        }
        
        
        const Tensor<2, dim> &get_Grad_U() const
        {
            return Grad_U;
        }
        
        // ...and the kinetic variables.  These are used in the material and
        // global tangent matrix and residual assembly operations:
        double get_p() const
        {
            return material->get_p();
        }
        
        const Tensor<2, dim> &get_C() const
        {
            return C;
        }
        
        const Tensor<2, dim> &get_C_inv() const
        {
            return C_inv;
        }
        
        double get_c_1() const
        {
            return material->get_c_1();
        }
        
        // In terms of member functions, this class stores for the quadrature
        // point it represents a copy of a material type in case different
        // materials are used in different regions of the domain, as well as the
        // inverse of the deformation gradient...
    private:
        Material_Incompressible_Neo_Hook_Two_Field<dim> *material;
        Tensor<2, dim> C;
        Tensor<2, dim> C_inv;
        Tensor<2, dim> Grad_U;
    };
    
    // The next three classes are used for the CG solver, and can be skipped if
    // one uses a direct solver. See step 20 for details.
    template <class MatrixType>
    class InverseMatrix : public Subscriptor
    {
    public:
        InverseMatrix(const MatrixType &m);
        void vmult(Vector<double>       &dst,
        const Vector<double> &src) const;
    private:
        const SmartPointer<const MatrixType> matrix;
    };
    template <class MatrixType>
    InverseMatrix<MatrixType>::InverseMatrix (const MatrixType &m)
    :
    matrix (&m)
    {}
    template <class MatrixType>
    void InverseMatrix<MatrixType>::vmult (Vector<double>       &dst,
            const Vector<double> &src) const
    {
        SolverControl solver_control (std::max(src.size(), static_cast<std::size_t> (200)),
                1e-8*src.l2_norm());
        SolverCG<>    cg (solver_control);
        dst = 0;
        std::cout << std::endl << "Try to inverse mass matrix..." << std::endl;
        cg.solve (*matrix, dst, src, PreconditionIdentity());
        std::cout << "Mass matrix inversed." << std::endl;
    }
    class SchurComplement : public Subscriptor
    {
    public:
        SchurComplement (const BlockSparseMatrix<double>            &A,
                const InverseMatrix<SparseMatrix<double> > &Minv);
        void vmult (Vector<double>       &dst,
        const Vector<double> &src) const;
    private:
        const SmartPointer<const BlockSparseMatrix<double> > tangent_matrix;
        const SmartPointer<const InverseMatrix<SparseMatrix<double> > > m_inverse;
        mutable Vector<double> tmp1, tmp2;
    };
    SchurComplement
            ::SchurComplement (const BlockSparseMatrix<double>            &A,
            const InverseMatrix<SparseMatrix<double> > &Minv)
    :
    tangent_matrix (&A),
            m_inverse (&Minv),
            tmp1 (A.block(0,0).m()),
            tmp2 (A.block(0,0).m())
    {}
    void SchurComplement::vmult (Vector<double>       &dst,
            const Vector<double> &src) const
    {
        tangent_matrix->block(0,1).vmult (tmp1, src);
        m_inverse->vmult (tmp2, tmp1);
        tangent_matrix->block(1,0).vmult (dst, tmp2);
    }
    class ApproximateSchurComplement : public Subscriptor
    {
    public:
        ApproximateSchurComplement (const BlockSparseMatrix<double> &A);
        void vmult (Vector<double>       &dst,
        const Vector<double> &src) const;
    private:
        const SmartPointer<const BlockSparseMatrix<double> > tangent_matrix;
        mutable Vector<double> tmp1, tmp2;
    };
    ApproximateSchurComplement::ApproximateSchurComplement
    (const BlockSparseMatrix<double> &A) :
    tangent_matrix (&A),
            tmp1 (A.block(0,0).m()),
            tmp2 (A.block(0,0).m())
    {}
    void
    ApproximateSchurComplement::vmult
    (Vector<double>       &dst,
            const Vector<double> &src) const
    {
        tangent_matrix->block(0,1).vmult (tmp1, src);
        tangent_matrix->block(0,0).precondition_Jacobi (tmp2, tmp1);
        tangent_matrix->block(1,0).vmult (dst, tmp2);
    }
    
    
    
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
        
        std::vector<types::global_dof_index> const&
        get_dofs_per_block()
        {
            return dofs_per_block;
        }
        
        void
        set_solution(BlockVector<double> const &solution);
        
        BlockVector<double> const&
        get_solution()
        {
            return solution_n;
        }
        
        void
        get_rhs_and_tangent(BlockVector<double> const* &sys_rhs,
        BlockSparseMatrix<double> const* &tm);
        
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
        
        struct PerTaskData_UQPH;
        struct ScratchData_UQPH;
        
        // We start the collection of member functions with one that builds the
        // grid:
        void
        make_grid();
        
        // Set up the finite element system to be solved:
        void
        system_setup();
        
        void
        determine_component_extractors();
        
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
        
        // Apply Dirichlet boundary conditions on the displacement field
        void
        make_constraints(const int &it_nr);
        
        // Create and update the quadrature points. Here, no data needs to be
        // copied into a global object, so the copy_local_to_global function is
        // empty:
        void
        setup_qph();
        
        void
        update_qph_incremental(const BlockVector<double> &solution_delta);
        
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
        solve_nonlinear_timestep(BlockVector<double> &solution_delta);
        
        void
        solve_linear_system(BlockVector<double> &newton_update, const int &it_nr);
        
        // Solution retrieval as well as post-processing and writing data to file:
        BlockVector<double>
        get_total_solution(const BlockVector<double> &solution_delta) const;
        
        void
        output_results() const;
        
        // Finally, some member variables that describe the current state: A
        // collection of the parameters used to describe the problem setup...
        Parameters::AllParameters        parameters;
        
        // ...the volume of the reference and current configurations...
        double                           vol_reference;
        double                           vol_current;
        
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
        const FEValuesExtractors::Vector u_fe;
        const FEValuesExtractors::Scalar p_fe;
        
        // Description of how the block-system is arranged. There are 2 blocks,
        // the first contains a vector DOF $\mathbf{u}$ while the other one
        // describes the scalar DOF, p.
        static const unsigned int        n_blocks = 2;
        static const unsigned int        n_components = dim + 1;
        static const unsigned int        first_u_component = 0;
        static const unsigned int        p_component = dim;
        
        enum
        {
            u_dof = 0,
            p_dof = 1,
        };
        
        std::vector<types::global_dof_index>        dofs_per_block;
        std::vector<types::global_dof_index>        element_indices_u;
        std::vector<types::global_dof_index>        element_indices_p;
        
        // Rules for Gauss-quadrature on both the cell and faces. The number of
        // quadrature points on both cells and faces is recorded.
        const QGauss<dim>                qf_cell;
        const QGauss<dim - 1>            qf_face;
        const unsigned int               n_q_points;
        const unsigned int               n_q_points_f;
        
        // Objects that store the converged solution and right-hand side vectors,
        // as well as the tangent matrix. There is a ConstraintMatrix object used
        // to keep track of constraints.  We make use of a sparsity pattern
        // designed for a block system.
        ConstraintMatrix                 constraints;
        BlockSparsityPattern             sparsity_pattern;
        BlockSparseMatrix<double>        tangent_matrix;
        BlockVector<double>              system_rhs;
        BlockVector<double>              solution_n;
        
        //Some boolean to decide if we want to display and save the tangent matrix and RHS
        const bool                       print_tangent_matrix = false;
        const bool                       print_tangent_matrix_constrained = false;
        const bool                       print_RHS = false;
        //This is just to know the size of the tangent matrix "by hand", for debugging:
    	const int                        dim_matrix = 150;
        
        // Then define a number of variables to store norms and update norms and
        // normalisation factors.
        struct Errors
        {
            Errors()
            :
            norm(1.0), u(1.0), p(1.0)
            {}
            
            void reset()
            {
                norm = 1.0;
                u = 1.0;
                p = 1.0;
            }
            void normalise(const Errors &rhs)
            {
                if (rhs.norm != 0.0)
                    norm /= rhs.norm;
                if (rhs.u != 0.0)
                    u /= rhs.u;
                if (rhs.p != 0.0)
                    p /= rhs.p;
            }
            
            double norm, u, p;
        };
        
        Errors error_residual, error_residual_0, error_residual_norm, error_update,
        error_update_0, error_update_norm;
        
        // Methods to calculate error measures
        void get_error_residual(Errors &error_residual);
        
        void get_error_update(const BlockVector<double> &newton_update,
        Errors &error_update);
        
        std::pair<double, double>
        get_error_dilation();
        
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
            // The Finite Element System is composed of dim continuous displacement
            // DOFs, and discontinuous pressure DOF. In an attempt to
            // satisfy the Babuska-Brezzi or LBB stability conditions (see Hughes
            // (2000)), we setup a $Q_n \times DGP_{n-1}$
            // system. $Q_2 \times DGPM_1$ elements satisfy this
            // condition, while $Q_1 \times DGP_0$ elements do
            // not. However, it has been shown that the latter demonstrate good
            // convergence characteristics nonetheless.
            
            fe(FE_Q<dim>(parameters.poly_degree), dim, // displacement
            FE_DGP<dim>(parameters.poly_degree-1),1), //pressure
            dof_handler_ref(triangulation),
            dofs_per_cell (fe.dofs_per_cell),
            u_fe(first_u_component),
            p_fe(p_component),
            dofs_per_block(n_blocks),
            qf_cell(parameters.quad_order),
            qf_face(parameters.quad_order),
            n_q_points (qf_cell.size()),
            n_q_points_f (qf_face.size())
    {
        determine_component_extractors();
        
        make_grid();
        system_setup();
        //        {
        //            ConstraintMatrix constraints;
        //            constraints.close();
        //        }
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
    //
    // We start the function with preprocessing, setting the initial dilatation
    // values, and then output the initial grid before starting the simulation
    //  proper with the first time (and loading)
    // increment.
    
    template <int dim>
    void Solid<dim>::run()
    {
        // Declare the incremental solution update $\varDelta
        // \mathbf{\Xi}:= \{\varDelta \mathbf{u},\varDelta p\}$ and
        // start the loop over the time domain.
        //
        // At the beginning, we reset the solution update for this time step...
        BlockVector<double> solution_delta(dofs_per_block);
        while (time.current() < time.end())
        {
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
        std::vector<std::vector<Tensor<2, dim> > >          grad_Nx;
        
        ScratchData_K(const FiniteElement<dim> &fe_cell,
        const QGauss<dim> &qf_cell,
        const UpdateFlags uf_cell)
        :
        fe_values_ref(fe_cell, qf_cell, uf_cell),
        Nx(qf_cell.size(),
        std::vector<double>(fe_cell.dofs_per_cell)),
        grad_Nx(qf_cell.size(),
        std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell))
        {}
        
        ScratchData_K(const ScratchData_K &rhs)
        :
        fe_values_ref(rhs.fe_values_ref.get_fe(),
        rhs.fe_values_ref.get_quadrature(),
        rhs.fe_values_ref.get_update_flags()),
        Nx(rhs.Nx),
        grad_Nx(rhs.grad_Nx)
        {}
        
        void reset()
        {
            const unsigned int n_q_points = Nx.size();
            const unsigned int n_dofs_per_cell = Nx[0].size();
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
                Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
                Assert( grad_Nx[q_point].size() == n_dofs_per_cell,
                        ExcInternalError());
                for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
                {
                    Nx[q_point][k] = 0.0;
                    grad_Nx[q_point][k] = 0.0;
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
        FEValues<dim>     fe_values_ref;
        FEFaceValues<dim> fe_face_values_ref;
        
        std::vector<std::vector<double> >          Nx;
        std::vector<std::vector<Tensor<2, dim> > > grad_Nx;
        
        ScratchData_RHS(const FiniteElement<dim> &fe_cell,
        const QGauss<dim> &qf_cell, const UpdateFlags uf_cell,
        const QGauss<dim - 1> & qf_face, const UpdateFlags uf_face)
        :
        fe_values_ref(fe_cell, qf_cell, uf_cell),
        fe_face_values_ref(fe_cell, qf_face, uf_face),
        Nx(qf_cell.size(),
        std::vector<double>(fe_cell.dofs_per_cell)),
        grad_Nx(qf_cell.size(),
        std::vector<Tensor<2, dim> >
        (fe_cell.dofs_per_cell))
        {}
        
        ScratchData_RHS(const ScratchData_RHS &rhs)
        :
        fe_values_ref(rhs.fe_values_ref.get_fe(),
        rhs.fe_values_ref.get_quadrature(),
        rhs.fe_values_ref.get_update_flags()),
        fe_face_values_ref(rhs.fe_face_values_ref.get_fe(),
        rhs.fe_face_values_ref.get_quadrature(),
        rhs.fe_face_values_ref.get_update_flags()),
        Nx(rhs.Nx),
        grad_Nx(rhs.grad_Nx)
        {}
        
        void reset()
        {
            const unsigned int n_q_points      = Nx.size();
            const unsigned int n_dofs_per_cell = Nx[0].size();
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
                Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
                Assert( grad_Nx[q_point].size() == n_dofs_per_cell,
                        ExcInternalError());
                for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
                {
                    Nx[q_point][k] = 0.0;
                    grad_Nx[q_point][k] = 0.0;
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
        const BlockVector<double>   &solution_total;
        
        std::vector<Tensor<2, dim> > solution_grads_u_total;
        std::vector<double>          solution_values_p_total;
        
        FEValues<dim>                fe_values_ref;
        
        ScratchData_UQPH(const FiniteElement<dim> &fe_cell,
        const QGauss<dim> &qf_cell,
        const UpdateFlags uf_cell,
        const BlockVector<double> &solution_total)
        :
        solution_total(solution_total),
        solution_grads_u_total(qf_cell.size()),
        solution_values_p_total(qf_cell.size()),
        fe_values_ref(fe_cell, qf_cell, uf_cell)
        {}
        
        ScratchData_UQPH(const ScratchData_UQPH &rhs)
        :
        solution_total(rhs.solution_total),
        solution_grads_u_total(rhs.solution_grads_u_total),
        solution_values_p_total(rhs.solution_values_p_total),
        fe_values_ref(rhs.fe_values_ref.get_fe(),
        rhs.fe_values_ref.get_quadrature(),
        rhs.fe_values_ref.get_update_flags())
        {}
        
        void reset()
        {
            const unsigned int n_q_points = solution_grads_u_total.size();
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                solution_grads_u_total[q] = 0.0;
                solution_values_p_total[q] = 0.0;
            }
        }
    };
    
    
    // @sect4{Solid::make_grid}
    
    // On to the first of the private member functions. Here we create the
    // triangulation of the domain, for which we choose the scaled cube with each
    // face given a boundary ID number.  The grid must be refined at least once
    // for the indentation problem.
    //
    // We then determine the volume of the reference configuration and print it
    // for comparison:
    template <int dim>
    void Solid<dim>::make_grid()
    {
        GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
        //GridTools::scale(parameters.scale, triangulation);
        
        // We will refine the grid in five steps towards the inner circle of
        // the domain. See step 1 for details.
        for (unsigned int step=0; step<parameters.local_refinement_cycles; ++step)
        {
            typename Triangulation<dim>::active_cell_iterator
            cell_ref = triangulation.begin_active(),
                    endc = triangulation.end();
            for (; cell_ref!=endc; ++cell_ref)
            {
                for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
                    if (cell_ref->face(f)->at_boundary())
                    {
                        const Point<dim> face_center = cell_ref->face(f)->center();
                        if (face_center[1] == 1) //Top
                            cell_ref->set_refine_flag ();
                    }
            }
            
            // Now that we have marked all the cells that we want refined, we let
            // the triangulation actually do this refinement.
            triangulation.execute_coarsening_and_refinement ();
            triangulation.clear_user_flags();
        }
        
        // We refine our mesh globally, at least once, to not too have a poor
        // refinement at the bottom of our square (like two cells)
        triangulation.refine_global(std::max (1U, parameters.global_refinement));
        
        vol_reference = GridTools::volume(triangulation);
        vol_current = vol_reference;
        std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
        
        // We mark the surfaces in order to apply the boundary conditions after
        typename Triangulation<dim>::active_cell_iterator cell =
        triangulation.begin_active(), endc = triangulation.end();
        for (; cell != endc; ++cell)
            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->face(f)->at_boundary())
                {
                    const Point<dim> face_center = cell->face(f)->center();
                    if (face_center[1] == 0) //Bottom
                        cell->face(f)->set_boundary_id (0);
                    else if (face_center[1] == 1) //Top
                        cell->face(f)->set_boundary_id (2);
                    else if (face_center[0] == 0) //Left
                        cell->face(f)->set_boundary_id (1);
                    else if (face_center[0] == 1) //Right
                        cell->face(f)->set_boundary_id (3);
                    else                          //No way
                        cell->face(f)->set_boundary_id (4);
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
        
        std::vector<unsigned int> block_component(n_components, u_dof); // Displacement
        block_component[p_component] = p_dof; // Pressure
        
        // The DOF handler is then initialised and we renumber the grid in an
        // efficient manner. We also record the number of DOFs per block.
        dof_handler_ref.distribute_dofs(fe);
        DoFRenumbering::Cuthill_McKee(dof_handler_ref);
        DoFRenumbering::component_wise(dof_handler_ref, block_component);
        DoFTools::count_dofs_per_block(dof_handler_ref, dofs_per_block,
        block_component);
        
        
        // Setup the sparsity pattern and tangent matrix
        tangent_matrix.clear();
        {
            const types::global_dof_index n_dofs_u = dofs_per_block[u_dof];
            const types::global_dof_index n_dofs_p = dofs_per_block[p_dof];
            
            //We print the number of cells and dofs:
            std::cout << std::endl;
            std::cout << "Triangulation:"
                    << "\n\t Number of active cells: " << triangulation.n_active_cells()
                    << "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
                    << " (" << n_dofs_u << " for u, "
                    << n_dofs_p << " for p)"
                    << "\n\t Number of dofs per cell: " << fe.dofs_per_cell
                    << std::endl;
            
            
            BlockDynamicSparsityPattern dsp(n_blocks, n_blocks);
            
            dsp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
            dsp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
            
            dsp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
            dsp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
            
            dsp.collect_sizes();
            
            // We optimise the sparsity pattern of the global system matrix
            // to reflect its structure and prevent unnecessary data creation
            // for the right-diagonal block components.
            Table<2, DoFTools::Coupling> coupling(n_components, n_components);
            for (unsigned int ii = 0; ii < n_components; ++ii)
                for (unsigned int jj = 0; jj < n_components; ++jj)
                    if  ((ii == p_component) && (jj == p_component))
                        coupling[ii][jj] = DoFTools::none;
                    else
                        coupling[ii][jj] = DoFTools::always;
                
                DoFTools::make_hanging_node_constraints (dof_handler_ref, constraints);
                
                DoFTools::make_sparsity_pattern(dof_handler_ref,
                        coupling,
                        dsp,
                        constraints,
                        true);
                sparsity_pattern.copy_from(dsp);
        }
        
        tangent_matrix.reinit(sparsity_pattern);
        
        // We then set up storage vectors
        system_rhs.reinit(dofs_per_block);
        system_rhs.collect_sizes();
        
        solution_n.reinit(dofs_per_block);
        solution_n.collect_sizes();
        
        //We initialize the pressure to c_1, i.e. mu/2
        for(unsigned int i = solution_n.block(u_dof).size(); i<solution_n.size(); ++i)
        {
            solution_n(i) = parameters.mu / 2.0;
        }
        
        
        // ...and finally set up the quadrature point history:
        setup_qph();
        
        timer.leave_subsection();
    }
    
    
    // @sect4{Solid::determine_component_extractors}
    // Next we compute some information from the FE system that describes which local
    // element DOFs are attached to which block component.  This is used later to
    // extract sub-blocks from the global matrix.
    //
    // In essence, all we need is for the FESystem object to indicate to which
    // block component a DOF on the reference cell is attached to.  Currently, the
    // interpolation fields are setup such that 0 indicates a displacement DOF and 1
    // a pressure DOF.
    template <int dim>
    void
    Solid<dim>::determine_component_extractors()
    {
        element_indices_u.clear();
        element_indices_p.clear();
        
        for (unsigned int k = 0; k < fe.dofs_per_cell; ++k)
        {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;
            if (k_group == u_dof)
                element_indices_u.push_back(k);
            else if (k_group == p_dof)
                element_indices_p.push_back(k);
            else
            {
                Assert(k_group <= p_dof, ExcInternalError());
            }
        }
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
        
        // Next we setup the initial quadrature
        // point data:
        for (typename Triangulation<dim>::active_cell_iterator cell =
                triangulation.begin_active(); cell != triangulation.end(); ++cell)
        {
            PointHistory<dim> *lqph =
                    reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());
            
            Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
            Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());
            
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                lqph[q_point].setup_lqp(parameters);
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
    void Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta)
    {
        timer.enter_subsection("Update QPH data");
        std::cout << " UQPH " << std::flush;
        
        const BlockVector<double> solution_total(get_total_solution(solution_delta));
        
        const UpdateFlags uf_UQPH(update_values | update_gradients);
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
        
        Assert(scratch.solution_grads_u_total.size() == n_q_points,
        ExcInternalError());
        Assert(scratch.solution_values_p_total.size() == n_q_points,
        ExcInternalError());
        
        scratch.reset();
        
        // We first need to find the values and gradients at quadrature points
        // inside the current cell and then we update each local QP using the
        // displacement gradient and total pressure and dilatation solution
        // values:
        scratch.fe_values_ref.reinit(cell);
        scratch.fe_values_ref[u_fe].get_function_gradients(scratch.solution_total,
        scratch.solution_grads_u_total);
        scratch.fe_values_ref[p_fe].get_function_values(scratch.solution_total,
        scratch.solution_values_p_total);
        
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            lqph[q_point].update_values(scratch.solution_grads_u_total[q_point],
                    scratch.solution_values_p_total[q_point]);
    }
    
    
    // @sect4{Solid::solve_nonlinear_timestep}
    
    // The next function is the driver method for the Newton-Raphson scheme. At
    // its top we create a new vector to store the current Newton update step,
    // reset the error storage objects and print solver header.
    template <int dim>
    void
    Solid<dim>::solve_nonlinear_timestep(BlockVector<double> &solution_delta)
    {
        std::cout << std::endl << "Timestep " << time.get_timestep() << " @ "
                << time.current() << "s" << std::endl;
        
        BlockVector<double> newton_update(dofs_per_block);
        
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
            
            if (newton_iteration > 0 && error_update_norm.u <= parameters.tol_u
                    && error_residual_norm.u <= parameters.tol_f)
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
                    << "  " << error_residual_norm.u << "  "
                    << error_residual_norm.p
                    << "  " << error_update_norm.norm << "  " << error_update_norm.u
                    << "  " << error_update_norm.p << "  " << std::endl;
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
        static const unsigned int l_width = 103;
        
        for (unsigned int i = 0; i < l_width; ++i)
            std::cout << "_";
        std::cout << std::endl;
        
        std::cout << "          SOLVER STEP            "
                << " | RES_NORM   "
                << " RES_U      RES_P     NU_NORM    "
                << " NU_U       NU_P " << std::endl;
        
        for (unsigned int i = 0; i < l_width; ++i)
            std::cout << "_";
        std::cout << std::endl;
    }
    
    
    
    template <int dim>
    void Solid<dim>::print_conv_footer()
    {
        static const unsigned int l_width = 103;
        
        for (unsigned int i = 0; i < l_width; ++i)
            std::cout << "_";
        std::cout << std::endl;
        
        const std::pair <double,double> error_dil = get_error_dilation();
        
        std::cout << "Relative errors:" << std::endl
                << "Displacement:\t" << error_update.u / error_update_0.u << std::endl
                << "Force: \t\t" << error_residual.u / error_residual_0.u << std::endl
                << "Dilatation:\t" << error_dil.first << std::endl
                << "v - V_0:\t" << vol_current << " - " << vol_reference
                << " = " << error_dil.second << std::endl;
    }
    
    
    // @sect4{Solid::get_error_dilation}
    
    // Calculate how well the dilatation $J$ agrees with $J = 1$ from the $L^2$
    // error $ \bigl[ \int_{\Omega_0} {[ J - 1 ]}^{2}\textrm{d}V \bigr]^{1/2}$.
    // We also return the ratio of the current volume of the
    // domain to the reference volume. This is of interest for incompressible
    // media where we want to check how well the isochoric constraint has been
    // enforced.
    template <int dim>
    std::pair<double, double>
    Solid<dim>::get_error_dilation()
    {
        double dil_L2_error = 0.0;
        vol_current = 0.0;
        
        FEValues<dim> fe_values_ref(fe, qf_cell, update_JxW_values);
        
        for (typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
                cell != triangulation.end(); ++cell)
        {
            fe_values_ref.reinit(cell);
            
            PointHistory<dim> *lqph =
                    reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());
            
            Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
            Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());
            
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
                const double det_F_qp = lqph[q_point].get_det_F();
                const double the_error_qp_squared = std::pow((det_F_qp - 1),
                        2);
                const double JxW = fe_values_ref.JxW(q_point);
                
                dil_L2_error += the_error_qp_squared * JxW;
                vol_current += det_F_qp * JxW;
            }
            Assert(vol_current > 0, ExcInternalError());
        }
        
        std::pair<double, double> error_dil;
        error_dil.first = std::sqrt(dil_L2_error);
        error_dil.second = vol_current - vol_reference;
        
        return error_dil;
    }
    
    
    // @sect4{Solid::get_error_residual}
    
    // Determine the true residual error for the problem.  That is, determine the
    // error in the residual for the unconstrained degrees of freedom.  Note that to
    // do so, we need to ignore constrained DOFs by setting the residual in these
    // vector components to zero.
    template <int dim>
    void Solid<dim>::get_error_residual(Errors &error_residual)
    {
        BlockVector<double> error_res(dofs_per_block);
        for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
            error_res(i) = (constraints.is_constrained(i)) ? 0.0 : system_rhs(i);
        
        error_residual.norm = error_res.l2_norm();
        error_residual.u = error_res.block(u_dof).l2_norm();
        error_residual.p = error_res.block(p_dof).l2_norm();
    }
    
    
    // @sect4{Solid::get_error_udpate}
    
    // Determine the true Newton update error for the problem
    template <int dim>
    void Solid<dim>::get_error_update(const BlockVector<double> &newton_update,
            Errors &error_update)
    {
        BlockVector<double> error_ud(dofs_per_block);
        for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
            if (!constraints.is_constrained(i))
                error_ud(i) = newton_update(i);
        
        error_update.norm = error_ud.l2_norm();
        error_update.u = error_ud.block(u_dof).l2_norm();
        error_update.p = error_ud.block(p_dof).l2_norm();
    }
    
    
    
    // @sect4{Solid::get_total_solution}
    
    // This function provides the total solution, which is valid at any Newton step.
    // This is required as, to reduce computational error, the total solution is
    // only updated at the end of the timestep.
    template <int dim>
    BlockVector<double>
    Solid<dim>::get_total_solution(const BlockVector<double> &solution_delta) const
    {
        BlockVector<double> solution_total(solution_n);
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
        std::cout << " ASM_K" << std::flush;
        
        const UpdateFlags uf_cell(update_values    |
        update_gradients |
        update_JxW_values);
        
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
                const unsigned int k_group = fe.system_to_base_index(k).first.first;
                
                if (k_group == u_dof)
                {
                    scratch.grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point);
                }
                else if (k_group == p_dof)
                    scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k, q_point);
                else
                    Assert(k_group <= p_dof, ExcInternalError());
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
            const double det_F                   = lqph[q_point].get_det_F();
            const double p                       = lqph[q_point].get_p();
            const Tensor<2, dim> C_inv           = lqph[q_point].get_C_inv();
            const double c_1                     = lqph[q_point].get_c_1();
            const Tensor<2, dim> Grad_U          = lqph[q_point].get_Grad_U();
            
            //_________Computing of d_rond_C_inv/d_rond_C :_____________
            Tensor<4, dim> d_C_inv_d_C;
	    for (unsigned int i=0; i<dim; ++i){
                for (unsigned int j=0; j<dim; ++j){
                    for (unsigned int k=0; k<dim; ++k){
                        for (unsigned int l=0; l<dim; ++l){
                            d_C_inv_d_C[i][j][k][l] = -( C_inv[i][k] * C_inv[l][j] + C_inv[i][l]
                                    * C_inv[k][j] ) / 2;
                        }}}}
            
            // Next we define some aliases to make the assembly process easier to
            // follow
            const std::vector<double>
            &N = scratch.Nx[q_point];
            const std::vector<Tensor<2, dim> >
            &grad_Nx = scratch.grad_Nx[q_point];
            const double JxW = scratch.fe_values_ref.JxW(q_point);
            
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const unsigned int component_i = fe.system_to_component_index(i).first;
                const unsigned int i_group     = fe.system_to_base_index(i).first.first;
                
                for (unsigned int j = 0; j <= i; ++j)
                {
                    const unsigned int component_j = fe.system_to_component_index(j).first;
                    const unsigned int j_group     = fe.system_to_base_index(j).first.first;
                    
                    // This is the Kuu contribution
                    
                    if ((i_group == j_group) && (i_group == u_dof))
                    {
                        //_______________Term 1______________________________
                        if (component_i == component_j)
                            data.cell_matrix(i, j) += 2 * c_1 * grad_Nx[i][component_i]
                                    * grad_Nx[j][component_j] * JxW;
                        
                        
                        //_______________Term 4 (or 3)______________________________
                        for(unsigned int k = 0; k < dim; ++k)
                        {
                            data.cell_matrix(i, j) -= 2 * p * det_F * det_F * C_inv[k]
                                    * (transpose(grad_Nx[i]) * grad_Nx[j])[k] * JxW;
                        }
                        
                        
                        //_______________Term 2 (or 1)_______________________
                        {
                            Tensor<2, dim> dot_prod_1;
                            {
                                //______Construction of first dot product_________
                                Tensor<2, dim> tmp;
                                tmp = 2 * (grad_Nx[j] + transpose(grad_Nx[j]) * Grad_U);
                                for(unsigned int k = 0; k < dim; ++k)
                                    for(unsigned int l = 0; l < dim; ++l)
                                        for(unsigned int m = 0; m < dim; ++m) //todo : put symmetric tensors
                                            for(unsigned int n = 0; n < dim ; ++n)
                                                dot_prod_1[k][l] +=  d_C_inv_d_C[k][l][m][n] * tmp[m][n];
                            }
                            {
                                //______Construction of second dot product________
                                Tensor<2, dim> tmp;
                                tmp = 2 * (grad_Nx[i] + transpose(grad_Nx[i]) * Grad_U);
                                for(unsigned int k = 0; k < dim; ++k)
                                    data.cell_matrix(i,j) -= p * det_F * det_F * dot_prod_1[k] * tmp[k] * JxW;
                            }
                        }
                        //_______________Term 3 (or 2)________________________
                        {
                            double dot_prod_1 = 0.0;
                            double dot_prod_2 = 0.0;
                            Tensor<2, dim> tmp1;
                            Tensor<2, dim> tmp2;
                            tmp1 = 2 * (grad_Nx[i] + transpose(grad_Nx[i]) * Grad_U);
                            tmp2 = 2 * (grad_Nx[j] + transpose(grad_Nx[j]) * Grad_U);
                            for(unsigned int k = 0; k < dim; ++k)
                            {
                                dot_prod_1 += C_inv[k] * tmp1[k];
                                dot_prod_2 += C_inv[k] * tmp2[k];
                            }
                            data.cell_matrix(i, j) -= p * det_F * det_F * dot_prod_1 * dot_prod_2 * JxW;
                        }
                    }
                    
                    //_____________________________ Contribution K_pu ____________________________
                    
                    else if ((i_group == p_dof) && (j_group == u_dof))
                    {
                        Tensor<2, dim> tmp;
                        tmp = 2 * (grad_Nx[j] + transpose(grad_Nx[j]) * Grad_U);
                        for(unsigned int k = 0; k < dim; ++k)
                            data.cell_matrix(i, j) -= N[i] *  det_F * det_F
                                    * C_inv[k] * tmp[k]
                                    * JxW;
                    }
                    
                    //_____________________________ Contribution K_up ____________________________
                    
                    // We need to construct both of them because we only
                    // compute that for $j<=i$ and we don't know how the
                    // u and p dofs are ordered in the cell, although they
                    // are ordered globally
                    
                    else if ((i_group == u_dof) && (j_group == p_dof))
                    {
                        Tensor<2, dim> tmp;
                        tmp = 2 * (grad_Nx[i] + transpose(grad_Nx[i]) * Grad_U);
                        for(unsigned int k = 0; k < dim; ++k)
                            data.cell_matrix(i, j) -= N[j] *  det_F * det_F
                                    * C_inv[k] * tmp[k]
                                    * JxW;
                    }
                    else
                        Assert((i_group <= p_dof) && (j_group <= p_dof),
                                ExcInternalError());
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
        std::cout << " ASM_R" << std::flush;
        
        const UpdateFlags uf_cell(update_values |
        update_gradients |
        update_JxW_values);
        const UpdateFlags uf_face(update_values |
        update_normal_vectors |
        update_JxW_values);
        
        PerTaskData_RHS per_task_data(dofs_per_cell);
        ScratchData_RHS scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);
        
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
                const unsigned int k_group = fe.system_to_base_index(k).first.first;
                
                if (k_group == u_dof)
                    scratch.grad_Nx[q_point][k]
                            = scratch.fe_values_ref[u_fe].gradient(k, q_point);
                else if (k_group == p_dof)
                    scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k, q_point);
                else
                    Assert(k_group <= p_dof, ExcInternalError());
            }
        }
        
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
            const double det_F                   = lqph[q_point].get_det_F();
            const double p                       = lqph[q_point].get_p();
            const Tensor<2, dim> C_inv           = lqph[q_point].get_C_inv();
            const double c_1                     = lqph[q_point].get_c_1();
            const Tensor<2, dim> Grad_U          = lqph[q_point].get_Grad_U();
            
            const std::vector<double>
            &N = scratch.Nx[q_point];
            const std::vector<Tensor<2, dim> >
            &Grad_Nx = scratch.grad_Nx[q_point];
            const double JxW = scratch.fe_values_ref.JxW(q_point);
            
            // We first compute the contributions
            // from the internal forces.  Note, by
            // definition of the rhs as the negative
            // of the residual, these contributions
            // are subtracted.
            
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const unsigned int i_group = fe.system_to_base_index(i).first.first;
                
                if (i_group == u_dof){
                    Tensor<2, dim> tmp;
                    tmp = Grad_Nx[i] + transpose(Grad_Nx[i]) + transpose(Grad_Nx[i]) * Grad_U
                            + transpose(Grad_U) * Grad_Nx[i];
                    data.cell_rhs(i) -= c_1 * trace(tmp) * JxW;
                    
                    Tensor<2, dim> tmp2;
                    tmp2 = 2 * (Grad_Nx[i] + transpose(Grad_Nx[i]) * Grad_U);
                    for(unsigned int k = 0; k < dim; ++k)
                    {
                        data.cell_rhs(i) += p * det_F * det_F * C_inv[k] * tmp2[k] * JxW;
                    }
                }
                else if (i_group == p_dof)
                    data.cell_rhs(i) += N[i] * (det_F*det_F - 1) * JxW;
                else
                    Assert(i_group <= p_dof, ExcInternalError());
            }
        }
        
        // Next we assemble the Neumann contribution. We first check to see it the
        // cell face exists on a boundary on which a traction is applied and add
        // the contribution if this is the case.
        // for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
        //      ++face)
        //   if (cell->face(face)->at_boundary() == true
        //       && cell->face(face)->boundary_id() == 1)
        //     {
        //       scratch.fe_face_values_ref.reinit(cell, face);
        
        //       for (unsigned int f_q_point = 0; f_q_point < n_q_points_f;
        //            ++f_q_point)
        //         {
        //           const Tensor<1, dim> &N =
        //             scratch.fe_face_values_ref.normal_vector(f_q_point);
        
        //           // Using the face normal at this quadrature point we specify the
        //           // traction in reference configuration. For this problem, a
        //           // defined pressure is applied in the reference configuration.
        //           // The direction of the applied traction is assumed not to
        //           // evolve with the deformation of the domain. The traction is
        //           // defined using the first Piola-Kirchhoff stress is simply
        //           // $\mathbf{t} = \mathbf{P}\mathbf{N} = [p_0 \mathbf{I}]
        //           // \mathbf{N} = p_0 \mathbf{N}$ We use the time variable to
        //           // linearly ramp up the pressure load.
        //           //
        //           // Note that the contributions to the right hand side vector we
        //           // compute here only exist in the displacement components of the
        //           // vector.
        //           // static const double  p0        = -4.0
        //           //                                  /
        //           //                                  (parameters.scale * parameters.scale);
        //           // const double         time_ramp = (time.current() / time.end());
        //           // const double         pressure  = p0 * parameters.p_p0 * time_ramp;
        //           // const Tensor<1, dim> traction  = pressure * N;
        
        // for (unsigned int i = 0; i < dofs_per_cell; ++i)
        //   {
        //     const unsigned int i_group =
        //       fe.system_to_base_index(i).first.first;
        
        //     if (i_group == u_dof)
        //       {
        //         const unsigned int component_i =
        //           fe.system_to_component_index(i).first;
        //         const double Ni =
        //           scratch.fe_face_values_ref.shape_value(i,
        //                                                  f_q_point);
        //         const double JxW = scratch.fe_face_values_ref.JxW(
        //                              f_q_point);
        
        //         data.cell_rhs(i) += (Ni * traction[component_i])
        //                             * JxW;
        //       }
        //   }
        //}
        //    }
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
    {
        std::cout << " CST" << std::flush;
        
        // Since the constraints are different at different Newton iterations, we
        // need to clear the constraints matrix and completely rebuild
        // it. However, after the first iteration, the constraints remain the same
        // and we can simply skip the rebuilding step if we do not clear it.
        if (it_nr > 1)
            return;
        
        constraints.clear();
        
        DoFTools::make_hanging_node_constraints (dof_handler_ref, constraints);
        
        const bool apply_dirichlet_bc = (it_nr == 0);
        
        // In the following, we will have to tell the function interpolation
        // boundary values which components of the solution vector should be
        // constrained (i.e., whether it's the x-, y-displacements or
        // combinations thereof). This is done using ComponentMask objects (see
        // @ref GlossComponentMask) which we can get from the finite element if we
        // provihde it with an extractor object for the component we wish to
        // select. To this end we first set up such extractor objects and later
        // use it when generating the relevant component masks:
        const FEValuesExtractors::Scalar x_displacement(0);
        const FEValuesExtractors::Scalar y_displacement(1);
        
        {
            const int boundary_id = 0;
            
            VectorTools::interpolate_boundary_values(dof_handler_ref,
                    boundary_id,
                    ZeroFunction<dim>(n_components),
                    constraints,
                    fe.component_mask(y_displacement));
        }
        {
            const int boundary_id = 1;
            
            if(apply_dirichlet_bc)
                VectorTools::interpolate_boundary_values(dof_handler_ref,
                        boundary_id,
                        ConstantFunction<dim>(parameters.displacement_side_1,n_components),
                        constraints,
                        fe.component_mask(x_displacement));
            else
                VectorTools::interpolate_boundary_values(dof_handler_ref,
                        boundary_id,
                        ZeroFunction<dim>(n_components),
                        constraints,
                        fe.component_mask(x_displacement));
            
            VectorTools::interpolate_boundary_values(dof_handler_ref,
                    boundary_id,
                    ZeroFunction<dim>(n_components),
                    constraints,
                    fe.component_mask(y_displacement));
        }
        
        {
            const int boundary_id = 3;
            
            VectorTools::interpolate_boundary_values(dof_handler_ref,
                    boundary_id,
                    ZeroFunction<dim>(n_components),
                    constraints,
                    fe.component_mask(x_displacement));
            
            VectorTools::interpolate_boundary_values(dof_handler_ref,
                    boundary_id,
                    ZeroFunction<dim>(n_components),
                    constraints,
                    fe.component_mask(y_displacement));
        }
        constraints.close();
    }
    
    // @sect4{Solid::solve_linear_system}
    
    template <int dim>
    void Solid<dim>::solve_linear_system(BlockVector<double> &newton_update, const int &it_nr)
    {
        {
            
            timer.enter_subsection("Linear solver");
            std::cout << " SLV" << std::flush;
            
            
            if(print_tangent_matrix)
            {
                std::cout << std::endl << "Print tangent matrix after applying boundary conditions..." << std::endl;
                
                FullMatrix<double> tangent_matrix_f(dim_matrix);
                tangent_matrix_f.copy_from(tangent_matrix);
                
                std::cout << std::endl;
                
                for(int x=0;x<dim_matrix;x++)
                {
                    for(int y=0;y<dim_matrix;y++)
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
                
                FullMatrix<double> tangent_matrix_f(dim_matrix);
                tangent_matrix_f.copy_from(tangent_matrix);
                
                std::cout << std::endl << "Extracting the intersting part..." << std::endl << std::endl;
                
                //We open a file to write the tangent matrix in it
                std::string tangent_matrix_file_name = "Matrices_RHS/tangent_matrix_timestep="
                        +Utilities::int_to_string(time.get_timestep())+"_nr_it="+Utilities::int_to_string(it_nr)+".txt";
                std::ofstream tangent_matrix_file(tangent_matrix_file_name, std::ios::out | std::ios::trunc);
                
                if(tangent_matrix_file)  // If opening the file was done well
                {
                    for(int x=0;x<dim_matrix;x++)
                    {
                        if (!constraints.is_constrained(x))
                        {
                            for(int y=0;y<dim_matrix;y++)
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
                    
                    for(int x=0;x<dim_matrix;x++)
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
            
            
            if (parameters.type_lin == "CG")
            {
                InverseMatrix<SparseMatrix<double> > inverse_mass (tangent_matrix.block(u_dof,u_dof));
                Vector<double> tmp (newton_update.block(u_dof).size());
                {
                    SchurComplement schur_complement (tangent_matrix, inverse_mass);
                    Vector<double> schur_rhs (newton_update.block(p_dof).size());
                    inverse_mass.vmult (tmp, system_rhs.block(u_dof)); //This is not working
                    tangent_matrix.block(p_dof,u_dof).vmult (schur_rhs, tmp);
                    schur_rhs -= system_rhs.block(p_dof);
                    SolverControl solver_control (newton_update.block(p_dof).size(),
                            1e-12*schur_rhs.l2_norm());
                    SolverCG<> cg (solver_control);
                    ApproximateSchurComplement approximate_schur (tangent_matrix);
                    InverseMatrix<ApproximateSchurComplement> approximate_inverse
                    (approximate_schur);
                    cg.solve (schur_complement, newton_update.block(p_dof), schur_rhs,
                            approximate_inverse);
                    std::cout << solver_control.last_step()
                            << " CG Schur complement iterations to obtain convergence."
                            << std::endl;
                }
                {
                    tangent_matrix.block(u_dof,p_dof).vmult (tmp, newton_update.block(p_dof));
                    tmp *= -1;
                    tmp += system_rhs.block(u_dof);
                    inverse_mass.vmult (newton_update.block(u_dof), tmp);
                }
            }
            else
                if (parameters.type_lin == "Direct")
                {
                    // Otherwise if the problem is small enough, a direct solver can be utilised
                    
                    SparseDirectUMFPACK A_direct;
                    A_direct.initialize(tangent_matrix);
                    A_direct.vmult(newton_update, system_rhs);
                }
                else
                    Assert (false, ExcMessage("Linear solver type not implemented"));
            
            timer.leave_subsection();
	}
        
        // Now that we have the displacement update, distribute the constraints
        // back to the Newton update:
        constraints.distribute(newton_update);
        
        timer.enter_subsection("Linear solver postprocessing");
        std::cout << " PP" << std::flush;
        
        
        timer.leave_subsection();
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
        DataComponentInterpretation::component_is_part_of_vector);
        data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
        
        std::vector<std::string> solution_name(dim, "displacement");
        solution_name.push_back("pressure");
        
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
        
        std::ostringstream filename;
        filename << "Results/timestep-" << time.get_timestep() << ".vtk";
        
        std::ofstream output(filename.str().c_str());
        data_out.write_vtk(output);
        
    }
    
    template <int dim>
    std::size_t
    Solid<dim>::get_system_size()
    {
        return solution_n.size();
    }
    
    template <int dim>
    void
    Solid<dim>::set_solution(BlockVector<double> const &solution)
    {
        BlockVector<double> solution_delta(dofs_per_block);
        solution_delta = 0.0;
        
        solution_n = solution;
        update_qph_incremental(solution_delta);
    }
    
    template <int dim>
    void
    Solid<dim>::get_rhs_and_tangent(BlockVector<double> const* &sys_rhs,
            BlockSparseMatrix<double> const* &tm)
    {
        tangent_matrix = 0.0;
        system_rhs = 0.0;
        assemble_system_rhs();
        assemble_system_tangent();
        const int iter_val = 0;  // what value do we want here?
        make_constraints(iter_val);
        constraints.condense(tangent_matrix, system_rhs);
        
        sys_rhs = &system_rhs;
        tm = &tangent_matrix;
    }
    
#ifdef CREATE_LIBRARY
    template class StandardTensors<2>;
    template class Material_Incompressible_Neo_Hook_Two_Field<2>;
    template class PointHistory<2>;
    Solid<2> MyNeoHookean("parameters.prm");
    
    void
    set_solution(double const* const solution)
    {
        BlockVector<double> sol(MyNeoHookean.get_dofs_per_block());
        std::size_t size(MyNeoHookean.get_system_size());
        for (unsigned i = 0; i < size; ++i)
        {
            sol[i] = solution[i];
        }
        MyNeoHookean.set_solution(sol);
    }
    
    void
    get_rhs_and_tangent(double* const sys_rhs,
            double* const tm)
    {
        //MyNeoHookean.get_solution().print(std::cout);
        BlockVector<double> const* rhs;
        BlockSparseMatrix<double> const* tangent;
        MyNeoHookean.get_rhs_and_tangent(rhs,tangent);
        
        std::size_t size(MyNeoHookean.get_system_size());
        for (unsigned i = 0; i < size; ++i)
        {
            sys_rhs[i] = (*rhs)[i];
        }
        memset(tm, 0, size*size*sizeof(double));
        for (BlockSparseMatrix<double>::const_iterator itr = tangent->begin();
                itr != tangent->end(); ++itr)
        {
            tm[size*(itr->row()) + itr->column()] = itr->value();
        }
    }
    
    std::size_t get_system_size()
    {
        return MyNeoHookean.get_system_size();
    }
    
    void run()
    {
        MyNeoHookean.run();
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
    using namespace neo_hookean;
    
    try
    {
        Solid<2> solid_2d("parameters.prm");
        solid_2d.run();
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
