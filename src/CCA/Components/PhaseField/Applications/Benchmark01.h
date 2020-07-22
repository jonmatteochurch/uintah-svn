/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 * @file CCA/Components/PhaseField/Applications/Benchmark01.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Benchmark01_h
#define Packages_Uintah_CCA_Components_PhaseField_Benchmark01_h

#include <sci_defs/hypre_defs.h>

#ifdef HAVE_HYPRE
#   include <CCA/Components/Solvers/HypreSStruct/AdditionalEntriesP.h>
#   include <CCA/Components/Solvers/HypreSStruct/Solver.h>
#endif

#include <CCA/Components/PhaseField/DataTypes/ScalarProblem.h>
#include <CCA/Components/PhaseField/Applications/Application.h>
#include <CCA/Components/PhaseField/Views/View.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWView.h>
#include <CCA/Components/PhaseField/AMR/AMRInterpolator.h>
#include <CCA/Components/PhaseField/AMR/AMRRestrictor.h>

#include <CCA/Ports/Regridder.h>

#include <Core/Util/Factory/Implementation.h>
#include <Core/Grid/SimpleMaterial.h>
#include <Core/Grid/Variables/PerPatchVars.h>

namespace Uintah
{
namespace PhaseField
{

/**
 * @brief Benchmark I application: 2D Allen Cahn
 *
 * \f[
 * \dot u = \epsilon^2 \nabla^2 u - W^\prime(u)
 * \f]
 * where
 * \f[
 *  W (u) = \frac{1}{4} (u^2 − 1)^2
 * \f]
 * on \f$[0, 2\pi]^2\f$ with initial data
 * \f[
 * u_{|t=0} = \tanh \frac{\sqrt{(x-\pi)^2+(y-\pi)^2}-2}{\epsilon\sqrt{2}}
 * \f]
 * with periodic boundary conditions.
 *
 * The model parameters are:
 * - \f$ \epsilon \f$   interface width
 *
 * @tparam VAR type of variable representation
 * @tparam STN finite-difference stencil
 */

/**
 * SemiImplicit0
 * \f[
 * [1 - k\epsilon^2 \nabla^2] u_{k+1} = [1 + k - k u_k^2] u_k
 * \f]
 *
 * SemiImplicit1
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - k] u_{k+1} = [1 - k u_k^2] u_k
 * \f]
 *
 * SemiImplicit2
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - k + k u_k^2] u_{k+1} = u_k
 * \f]
 *
 * SemiImplicit3
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - k + 3k u_k^2] u_{k+1} = [1 + 2k u_k^2] u_k
 * \f]
 *
 * SemiImplicit4
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - \tfrac k2] u_{k+1} = [1 + \tfrac k2 - k u_k^2] u_k
 * \f]
 *
 * SemiImplicit5
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - \tfrac k2 + \tfrac k2 u_k^2] u_{k+1} = [1 + \tfrac k2 - \tfrac k2 u_k^2] u_k
 * \f]
 *
 * SemiImplicit6
 * \f[
 * [1 - k\epsilon^2 \nabla^2 - \tfrac k2 + \tfrac{3k}2 u_k^2] u_{k+1} = [1 + \tfrac k2 + \tfrac k2 u_k^2] u_k
 * \f]
 *
 * SemiImplicit7
 * \f[
 * [1 - \tfrac{k\epsilon^2}2\nabla^2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2\nabla^2 + k - ku_k^2] u_k
 * \f]
 *
 * SemiImplicit8
 * \f[
 * [1 - \tfrac{k\epsilon^2}2 \nabla^2 - k] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2 - k u_k^2] u_k
 * \f]
 *
 * SemiImplicit9
 * \f[
 * [1 - \tfrac{k\epsilon^2}2 \nabla^2 - k + k u_k^2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2] u_k
 * \f]
 *
 * SemiImplicit10
 * \f[
 * [1 - \tfrac{k\epsilon^2}2 \nabla^2 - k + 3k u_k^2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2 + 2ku_k^2] u_k
 * \f]
 *
 * SemiImplicit11
 * \f[
 * [1 - \tfrac{k\epsilon^2}2 \nabla^2 - \tfrac k2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2 + \tfrac k2 - k u_k^2] u_k
 * \f]
 *
 * SemiImplicit12
 * \f[
 * [1 - \tfrac{k\epsilon^2}2\nabla^2 - \tfrac k2 + \tfrac k2 u_k^2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2 + \tfrac k2 - \tfrac k2 u_k^2] u_k
 * \f]
 *
 * SemiImplicit13
 * \f[
 * [1 + \tfrac{k\epsilon^2}2 \nabla^2 - \tfrac k2 + \tfrac{3k}2 u_k^2] u_{k+1} = [1 + \tfrac{k\epsilon^2}2 \nabla^2 + \tfrac k2 + \tfrac k2 ku_k^2] u_k
 * \f]
 */

template<VarType VAR, StnType STN, bool AMR = false>
class Benchmark01
    : public Application< ScalarProblem<VAR, STN>, AMR >
    , public Implementation<Benchmark01<VAR, STN, AMR>, UintahParallelComponent, const ProcessorGroup *, const MaterialManagerP, const std::string &, int >
{
private: // STATIC MEMBERS

    /// Index for solution
    static constexpr size_t U = 0;

    /// ScalarProblem material index (only one SimpleMaterial)
    static constexpr int material = 0;

    /// ScalarProblem dimension
    static constexpr DimType DIM = D2;

    /// Number of ghost elements required by STN (on the same level)
    using Application< ScalarProblem<VAR, STN> >::FGN;

    /// Type of ghost elements required by VAR and STN (on coarser level)
    using Application< ScalarProblem<VAR, STN> >::FGT;

    /// Number of ghost elements required by STN (on coarser level)
    using Application< ScalarProblem<VAR, STN> >::CGN;

    /// Type of ghost elements required by VAR and STN (on the same level)
    using Application< ScalarProblem<VAR, STN> >::CGT;

    /// Interpolation type for refinement
    using Application< ScalarProblem<VAR, STN> >::C2F;

    /// Restriction type for coarsening
    using Application< ScalarProblem<VAR, STN> >::F2C;

#ifdef HAVE_HYPRE
    /// Non-stencil entries type for implicit matrix
    using _AdditionalEntries = PerPatch<HypreSStruct::AdditionalEntriesP>;
#endif

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const FactoryString Name;

protected: // MEMBERS

    /// Label for u field into the DataWarehouse
    const VarLabel * u_label;

    /// Label for solution value at domain center (\f$ \pi, \pi \f$) into the DataWarehouse
    const VarLabel * u0_label;

    /// Label for system energy into the DataWarehouse
    const VarLabel * energy_label;

#ifdef HAVE_HYPRE
    /// Label for the implicit matrix (stencil entries) in the DataWarehouse
    const VarLabel * matrix_label;

    /// Label for the implicit vector in the DataWarehouse
    const VarLabel * rhs_label;

    /// Label for the implicit matrix (non-stencil entries) in the DataWarehouse
    const VarLabel * additional_entries_label;
#endif // HAVE_HYPRE

    /// Time step size
    double delt;

    /// Interface width
    double epsilon;

    /// Threshold for AMR
    double refine_threshold;

#ifdef HAVE_HYPRE
    /// Time advance scheme
    TS time_scheme;
#endif

public: // CONSTRUCTORS/DESTRUCTOR

    /**
     * @brief Constructor
     *
     * Intantiate a Benchmark I application
     *
     * @param myworld data structure to manage mpi processes
     * @param materialManager data structure to manage materials
     * @param verbosity constrols amount of debugging output
     */
    Benchmark01 (
        const ProcessorGroup * myworld,
        const MaterialManagerP materialManager,
        const std::string & uda,
        int verbosity = 0
    );

    /**
     * @brief Destructor
     */
    virtual ~Benchmark01();

    /// Prevent copy (and move) constructor
    Benchmark01 ( const Benchmark01 & ) = delete;

    /// Prevent copy (and move) assignment
    /// @return deleted
    Benchmark01 & operator= ( const Benchmark01 & ) = delete;

protected: // SETUP

    /**
       * @brief Setup
       *
       * Initialize problem parameters with values from problem specifications
       *
       * @param params problem specifications parsed from input file
       * @param restart_prob_spec unused
       * @param grid unused
       */
    virtual void
    problemSetup (
        const ProblemSpecP & params,
        const ProblemSpecP & restart_prob_spec,
        GridP & grid
    ) override;

protected: // SCHEDULINGS

    /**
       * @brief Schedule the initialization tasks
       *
       * Specify all tasks to be performed at initial timestep to initialize
       * variables in the DataWarehouse
       *
       * @param level grid level to be initialized
       * @param sched scheduler to manage the tasks
       */
    virtual void
    scheduleInitialize (
        const LevelP & level,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_initialize_solution (non AMR implementation)
     *
     * Defines the dependencies and output of the task which initializes the
     * solution allowing sched to control its execution order
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleInitialize_solution (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_initialize_solution (AMR implementation)
     *
     * Defines the dependencies and output of the task which initializes the
     * solution allowing sched to control its execution order
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleInitialize_solution (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the initialization tasks for restarting a simulation
     *
     * Specify all tasks to be performed at fist timestep after a stop to
     * initialize not saved variables to the DataWarehouse
     *
     * @remark only subproblems need to be reinitialized all other variables
     * should be retrieved from saved checkpoints
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleRestartInitialize (
        const LevelP & level,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_compute_stable_timestep
     *
     * Specify all tasks to be performed before each time advance to compute a
     * timestep size which ensures numerical stability
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleComputeStableTimeStep (
        const LevelP & level,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule the time advance tasks
     *
     * Specify all tasks to be performed at each timestep to update the
     * simulation variables in the DataWarehouse
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleTimeAdvance (
        const LevelP & level,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_time_advance_solution
     *
     * Defines the dependencies and output of the task which updates u allowing
     * sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleTimeAdvance_solution (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_forward_euler
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which updates u
     * allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_forward_euler (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_forward_euler
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which updates u
     * allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_forward_euler (
        const LevelP & level,
        SchedulerP & sched
    );

#ifdef HAVE_HYPRE
    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble
     * (non AMR implementation)
     *
     * Switches between available implementations depending on the given time
     * solver
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble
     * (AMR implementation)
     *
     * Switches between available implementations depending on the given time
     * solver
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble_hypre
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble_hypre
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble_hypresstruct
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_semi_implicit_assemble_hypresstruct
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG, TS SI >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule semi-implicit solve task
     *
     * Defines the dependencies and output of the task which solve the implicit
     * system to update u
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleTimeAdvance_solve (
        const LevelP & level,
        SchedulerP & sched
    );
#endif

    /**
     * @brief Schedule task_time_advance_postprocess
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which updates energy and
     * solution value at the center allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_postprocess (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_postprocess
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which updates energy and
     * solution value at the center allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_postprocess (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the refinement tasks
     *
     * Specify all tasks to be performed after an AMR regrid in order to populate
     * variables in the DataWarehouse at newly created patches
     *
     * @remark If regridding happens at initial time step scheduleInitialize is
     * called instead
     *
     * @param new_patches patches to be populated
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleRefine (
        const PatchSet * new_patches,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_refine_solution
     *
     * Defines the dependencies and output of the task which interpolates the
     * solution from the coarser level to each one of the new_patches
     * allowing sched to control its execution order
     *
     * @param new_patches patches to be populated
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleRefine_solution (
        const PatchSet * new_patches,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the refinement tasks
     *
     * Do nothing
     *
     * @param level_fine unused
     * @param sched unused
     * @param need_old_coarse unused
     * @param need_new_coarse unused
     */
    virtual void
    scheduleRefineInterface (
        const LevelP & level_fine,
        SchedulerP & sched,
        bool need_old_coarse,
        bool need_new_coarse
    ) override;

    /**
     * @brief Schedule the time coarsen tasks
     *
     * Specify all tasks to be performed after each timestep to restrict the
     * computed variables from finer to coarser levels
     *
     * @param level_coarse level to be updated
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleCoarsen (
        const LevelP & level_coarse,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_coarsen_solution
     *
     * Defines the dependencies and output of the task which restrict the
     * solution to level_coarse from its finer level allowing sched to control
     * its execution order
     *
     * @param level_coarse level to be updated
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleCoarsen_solution (
        const LevelP & level_coarse,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the error estimate tasks
     *
     * Specify all tasks to be performed before each timestep to estimate the
     * spatial discretization error on the solution update in order to decide
     * where to refine the grid
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleErrorEstimate (
        const LevelP & level,
        SchedulerP & sched
    ) override;

    /**
     * @brief Schedule task_error_estimate_solution (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which estimates the
     * spatial discretization error allowing sched to controvaluel its execution order
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleErrorEstimate_solution (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_error_estimate_solution (refinement level implementation)
     *
     * Defines the dependencies and output of the task which estimates the
     * spatial discretization error allowing sched to control its execution order
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleErrorEstimate_solution (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the initial error estimate tasks
     *
     * Specify all tasks to be performed before the first timestep to estimate
     * the spatial discretization error on the solution update in order to decide
     * where to refine the grid
     *
     * @remark forward to scheduleErrorEstimate
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    virtual void
    scheduleInitialErrorEstimate (
        const LevelP & level,
        SchedulerP & sched
    ) override;

protected: // TASKS

    /**
     * @brief Initialize solution task
     *
     * Allocate and save variables u for each one of the patches
     * and save them to dw_new
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_initialize_solution (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Compute timestep task
     *
     * Puts into the new DataWarehouse the constant value specified in input (delt)
     * of the timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_compute_stable_timestep (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Advance solution task (Forward Euler implementation)
     *
     * Computes new value of u using the value of the solution and at
     * previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_forward_euler (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

#ifdef HAVE_HYPRE
    /**
     * @brief Assemble hypre system task (Semi Implicit implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypre (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Semi Implicit, full implementation)
     *
     * Assemble both semi-implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypre_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Semi Implicit, rhs implementation)
     *
     * Assemble only semi-implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypre_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Semi Implicit implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypresstruct (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Semi Implicit, full implementation)
     *
     * Assemble both semi-implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypresstruct_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Semi Implicit, rhs implementation)
     *
     * Assemble only semi-implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @tparam SI semi implicit scheme
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    template<TS SI>
    void
    task_time_advance_solution_semi_implicit_assemble_hypresstruct_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );
#endif

    /**
     * @brief Advance postprocess task
     *
     * Computes new value of the system energy and of the solution at the center
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_postprocess (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Refine solution task
     *
     * Computes interpolated value of u on new refined patched
     *
     * @param myworld data structure to manage mpi processes
     * @param patches_fine list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_refine_solution (
        const ProcessorGroup * myworld,
        const PatchSubset * patches_fine,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Coarsen solution task
     *
     * Restricted value of u from refined regions to coarse patches
     * underneath
     *
     * @param myworld data structure to manage mpi processes
     * @param patches_coarse list of patches to be updated
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_coarsen_solution (
        const ProcessorGroup * myworld,
        const PatchSubset * patches_coarse,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief ErrorEstimate solution task
     *
     * Computes the gradient of the solution using its value at the previous
     * timestep and set refinement flag where it is above the threshold given
     * in input
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_error_estimate_solution (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

protected: // IMPLEMENTATIONS

    /**
     * @brief Initialize solution implementation
     *
     * compute initial condition for u at a given grid position
     *
     * @param id grid index
     * @param patch grid patch
     * @param[out] u view of the solution field in the new dw
     */
    void
    initialize_solution (
        const IntVector & id,
        const Patch * patch,
        View < ScalarField<double> > & u
    );

    /**
     * @brief Advance solution implementation
     *
     * compute new value for u at a given grid position using the value of the
     * solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] u_new view of the solution field in the new dw
     */
    virtual void
    time_advance_solution_forward_euler (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & u_old,
        View< ScalarField<double> > & u_new
    );

#ifdef HAVE_HYPRE
    /**
     * @brief Assemble hypre system implementation
     * (SemiImplicit0, full implementation)
     *
     * compute both semi-implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * \f[
     * [1 - k\epsilon^2 \nabla^2] u_{k+1} = [1 + k - k u_k^2] u_k
     * \f]
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    template<TS SI>
    typename std::enable_if < SI == TS::SemiImplicit0, void >::type
    time_advance_solution_semi_implicit_assemble_hypre_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & u_old,
        View < ScalarField<Stencil7> > & A,
        View < ScalarField<double> > & b
    );

    /**
     * @brief Assemble hypre system implementation
     * (SemiImplicit0, rhs implementation)
     *
     * compute only semi-implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * \f[
     * [1 - k\epsilon^2 \nabla^2] u_{k+1} = [1 + k - k u_k^2] u_k
     * \f]
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    template<TS SI>
    typename std::enable_if < SI == TS::SemiImplicit0, void >::type
    time_advance_solution_semi_implicit_assemble_hypre_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & u_old,
        View < ScalarField<double> > & b
    );

    /**
     * @brief Assemble hypre sstruct system implementation
     * (SemiImplicit0, full implementation)
     *
     * compute both semi-implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * \f[
     * [1 - k\epsilon^2 \nabla^2] u_{k+1} = [1 + k - k u_k^2] u_k
     * \f]
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] A_additional view of the implicit matrix non-stencil per patch entries in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    template<TS SI>
    typename std::enable_if < SI == TS::SemiImplicit0, void >::type
    time_advance_solution_semi_implicit_assemble_hypresstruct_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & u_old,
        View < ScalarField<Stencil7> > & A,
        HypreSStruct::AdditionalEntries * A_additional,
        View < ScalarField<double> > & b
    );

    /**
     * @brief Assemble hypre sstruct system implementation
     * (SemiImplicit0, rhs implementation)
     *
     * compute only semi-implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * \f[
     * [1 - k\epsilon^2 \nabla^2] u_{k+1} = [1 + k - k u_k^2] u_k
     * \f]
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    template<TS SI>
    typename std::enable_if < SI == TS::SemiImplicit0, void >::type
    time_advance_solution_semi_implicit_assemble_hypresstruct_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & u_old,
        View < ScalarField<double> > & b
    );
#endif

    /**
     * @brief Advance system energy implementation
     *
     * compute new value for energy at a given grid position and accumulate it
     *
     * @param id grid index
     * @param patch grid patch
     * @param u_new view of the solution field in the new dw
     * @param[out] energy system energy
     */
    virtual void
    time_advance_postprocess_energy (
        const IntVector & id,
        const Patch * patch,
        const FDView < ScalarField<const double>, STN > & u_new,
        double & energy
    );

    /**
     * @brief Refine solution implementation
     *
     * Computes interpolated value of u at a given grid position

     * @param id_fine fine grid index
     * @param u_coarse_interp interpolator of the aolution field on the coarse level
     * @param[out] u_fine view of the aolution field on the fine level
     */
    void
    refine_solution (
        const IntVector id_fine,
        const View < ScalarField<const double> > & u_coarse_interp,
        View < ScalarField<double> > & u_fine
    );

    /**
     * @brief Coarsen solution implementation
     *
     * Computes restricted value of u at a given grid position

     * @param id_coarse coarse grid index
     * @param u_fine_restr restrictor of the aolution field on the fine level
     * @param[out] u_coarse view of the aolution field on the coarse level
     */
    void
    coarsen_solution (
        const IntVector id_coarse,
        const View < ScalarField<const double> > & u_fine_restr,
        View < ScalarField<double> > & u_coarse
    );

    /**
     * @brief ErrorEstimate solution implementation (cell centered implementation)
     *
     * Computes the gradient of the phase field using its value at the previous
     * timestep and set refinement flag where it is above the threshold given
     * in input
     *
     * @param id grid index
     * @param u view of the solution the old dw
     * @param[out] refine_flag view of refine flag (grid field) in the new dw
     * @param[out] refine_patch flag for patch refinement
     */
    template < VarType V >
    typename std::enable_if < V == CC, void >::type
    error_estimate_solution (
        const IntVector id,
        FDView < ScalarField<const double>, STN > & u,
        View < ScalarField<int> > & refine_flag,
        bool & refine_patch// SCHEDULINGS

    );

    /**
     * @brief ErrorEstimate solution implementation (vertex based implementation)
     *
     * Computes the gradient of the phase field using its value at the previous
     * timestep and set refinement flag where it is above the threshold given
     * in input
     *
     * @param id grid index
     * @param u view of the solution in the old dw
     * @param[out] refine_flag view of refine flag (grid field) in the new dw
     * @param[out] refine_patch flag for patch refinement
     */
    template < VarType V >
    typename std::enable_if < V == NC, void >::type
    error_estimate_solution (
        const IntVector id,
        FDView < ScalarField<const double>, STN > & u,
        View < ScalarField<int> > & refine_flag,
        bool & refine_patch
    );

}; // class Benchmark01

// CONSTRUCTORS/DESTRUCTOR

template < VarType VAR, StnType STN, bool AMR >
Benchmark01<VAR, STN, AMR>::Benchmark01 (
    const ProcessorGroup * myWorld,
    const MaterialManagerP materialManager,
    const std::string &,
    int verbosity
) : Application< ScalarProblem<VAR, STN>, AMR > ( myWorld, materialManager, verbosity )
{
    u_label = VarLabel::create ( "u", Variable<VAR, double>::getTypeDescription() );
    u0_label = VarLabel::create ( "u0", sum_vartype::getTypeDescription() );
    energy_label = VarLabel::create ( "energy", sum_vartype::getTypeDescription() );
#ifdef HAVE_HYPRE
    matrix_label = VarLabel::create ( "A", Variable<VAR, Stencil7>::getTypeDescription() );
    additional_entries_label = VarLabel::create ( "A" + HypreSStruct::Solver<DIM>::AdditionalEntriesSuffix, _AdditionalEntries::getTypeDescription() );
    rhs_label = VarLabel::create ( "b", Variable<VAR, double>::getTypeDescription() );
#endif
}

template < VarType VAR, StnType STN, bool AMR >
Benchmark01<VAR, STN, AMR>::~Benchmark01()
{
    VarLabel::destroy ( u_label );
    VarLabel::destroy ( u0_label );
    VarLabel::destroy ( energy_label );
#ifdef HAVE_HYPRE
    VarLabel::destroy ( matrix_label );
    VarLabel::destroy ( additional_entries_label );
    VarLabel::destroy ( rhs_label );
#endif
}

// SETUP

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::problemSetup (
    ProblemSpecP const & params,
    ProblemSpecP const &,
    GridP &
)
{
    this->m_materialManager->registerSimpleMaterial ( scinew SimpleMaterial() );

    ProblemSpecP benchmark = params->findBlock ( "PhaseField" );
    benchmark->require ( "delt", delt );
    benchmark->require ( "epsilon", epsilon );

    std::string scheme;
    benchmark->getWithDefault ( "scheme", scheme, "forward_euler" );
#ifdef HAVE_HYPRE
    time_scheme = str_to_ts ( scheme );

    if ( VAR == NC )
    {
        if ( time_scheme & TS::SemiImplicit )
            SCI_THROW ( InternalError ( "\n ERROR: semi-implicit solver not implemented for node centered variables", __FILE__, __LINE__ ) );
    }

    if ( time_scheme & TS::SemiImplicit )
    {
        ProblemSpecP solv = params->findBlock ( "Solver" );
        this->m_solver = dynamic_cast<SolverInterface *> ( this->getPort ( "solver" ) );
        if ( !this->m_solver )
        {
            SCI_THROW ( InternalError ( "Benchmark01:couldn't get solver port", __FILE__, __LINE__ ) );
        }
        this->m_solver->readParameters ( solv, "u" );
        this->m_solver->getParameters()->setSolveOnExtraCells ( false );
    }
#else
    if ( scheme != "forward_euler" )
        SCI_THROW ( InternalError ( "\n ERROR: SemiImplicit time schemes require HYPRE\n", __FILE__, __LINE__ ) );
#endif

    this->setBoundaryVariables ( u_label );

    if ( AMR )
    {
        this->setLockstepAMR ( true );

        // read amr parameters
        benchmark->require ( "refine_threshold", refine_threshold );

        std::map<std::string, FC> c2f;

        // default values
        c2f[u_label->getName()] = ( VAR == CC ) ? FC::FC0 : FC::FC1;

        ProblemSpecP amr, regridder, fci;
        if ( ! ( amr = params->findBlock ( "AMR" ) ) ) return;
        if ( ! ( fci = amr->findBlock ( "FineCoarseInterfaces" ) ) ) return;
        if ( ! ( fci = fci->findBlock ( "FCIType" ) ) ) return;
        do
        {
            std::string label, var;
            fci->getAttribute ( "label", label );
            fci->getAttribute ( "var", var );
            c2f[label] = str_to_fc ( var );
        }
        while ( ( fci = fci->findNextBlock ( "FCIType" ) ) );

        this->setC2F ( c2f );
    }
}

// SCHEDULINGS

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleInitialize (
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleInitialize_solution<AMR> ( level, sched );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleInitialize_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_initialize_solution", this, &Benchmark01::task_initialize_solution );
    task->computes ( u_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

/**
 * @remark we need to schedule all levels before task_error_estimate_solution to avoid
 * the error "Failure finding [u , coarseLevel, MI: none, NewDW
 * (mapped to dw index 1), ####] for Benchmark01::task_error_estimate_solution",
 * on patch #, Level #, on material #, resource (rank): #" while compiling the
 * TaskGraph
 */
template<VarType VAR, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleInitialize_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    // since the SimulationController is calling this scheduler starting from
    // the finest level we schedule only on the finest level
    if ( level->hasFinerLevel() ) return;

    GridP grid = level->getGrid();
    for ( int l = 0; l < grid->numLevels(); ++l )
        scheduleInitialize_solution < !MG > ( grid->getLevel ( l ), sched );
}

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::scheduleRestartInitialize (
    const LevelP & /*level*/,
    SchedulerP & /*sched*/
)
{
}

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::scheduleComputeStableTimeStep (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_compute_stable_timestep ", this, &Benchmark01::task_compute_stable_timestep );
    task->computes ( this->getDelTLabel(), level.get_rep() );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance (
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleTimeAdvance_solution ( level, sched );
    scheduleTimeAdvance_postprocess<AMR> ( level, sched );
};

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
#ifdef HAVE_HYPRE
    switch ( time_scheme )
    {
    case TS::ForwardEuler:
#endif
        scheduleTimeAdvance_solution_forward_euler<AMR> ( level, sched );
#ifdef HAVE_HYPRE
        break;
    case TS::SemiImplicit0:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit0> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
#if 0
    case TS::SemiImplicit1:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit1> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit2:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit2> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit3:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit3> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit4:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit4> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit5:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit5> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit6:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit6> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit7:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit7> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit8:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit8> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit9:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit9> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit10:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit10> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit11:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit11> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit12:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit12> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
    case TS::SemiImplicit13:
        scheduleTimeAdvance_solution_semi_implicit_assemble<AMR, TS::SemiImplicit13> ( level, sched );
        scheduleTimeAdvance_solve ( level, sched );
        break;
#endif
    default:
        SCI_THROW ( InternalError ( "\n ERROR: Unknown time scheme\n", __FILE__, __LINE__ ) );
    }
#endif
}

template < VarType VAR, StnType STN, bool AMR >
template < bool MG >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_forward_euler (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_solution_forward_euler", this, &Benchmark01<VAR, STN, AMR>::task_time_advance_solution_forward_euler );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->computes ( u_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template < VarType VAR, StnType STN, bool AMR >
template < bool MG >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_forward_euler (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( !level->hasCoarserLevel() )
        scheduleTimeAdvance_solution_forward_euler < !MG > ( level, sched );
    else
    {
        Task * task = scinew Task ( "Benchmark01:task_time_advance_solution_forward_euler", this, &Benchmark01::task_time_advance_solution_forward_euler );
        task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
        task->requires ( Task::OldDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::OldDW, u_label, FGT, FGN );
        task->requires ( Task::OldDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->computes ( u_label );
        sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
    }
}

#ifdef HAVE_HYPRE
template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble
(
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( this->m_solver->getName() == "hypre" )
        scheduleTimeAdvance_solution_semi_implicit_assemble_hypre<AMR, SI> ( level, sched );
    else
        SCI_THROW ( InternalError ( "\n ERROR: Unsupported implicit solver\n", __FILE__, __LINE__ ) );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble
(
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( this->m_solver->getName() == "hypre" )
    {
        if ( !level->hasCoarserLevel() )
            scheduleTimeAdvance_solution_semi_implicit_assemble_hypre < !MG, SI > ( level, sched );
        else
            scheduleTimeAdvance_solution_semi_implicit_assemble_hypre < MG, SI > ( level, sched );
    }
    else if ( this->m_solver->getName() == "hypre_sstruct" )
    {
        if ( level->hasCoarserLevel() ) return;

        GridP grid = level->getGrid();

        // all assemble task must be sent to the scheduler before the solve task
        scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct < !MG, SI > ( level, sched );
        for ( int l = 1; l < grid->numLevels(); ++l )
            scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct < MG, SI > ( grid->getLevel ( l ), sched );
    }
    else
        SCI_THROW ( InternalError ( "\n ERROR: Unsupported implicit solver\n", __FILE__, __LINE__ ) );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre", this, &Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre<SI> );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::OldDW, matrix_label, Ghost::None, 0 );
    task->computes ( matrix_label );
    task->computes ( rhs_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre", this, &Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre<SI> );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, matrix_label, Ghost::None, 0 );
    task->computes ( matrix_label );
    task->computes ( rhs_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct", this, &Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct<SI> );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::OldDW, matrix_label, Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label, Ghost::None, 0 );
    task->computes ( matrix_label );
    task->computes ( additional_entries_label );
    task->computes ( rhs_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG, TS SI >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solution_semi_implicit_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct", this, &Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct<SI> );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), 0, Task::CoarseLevel, 0, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::OldDW, matrix_label, Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label, Ghost::None, 0 );
    task->computes ( matrix_label );
    task->computes ( additional_entries_label );
    task->computes ( rhs_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_solve
(
    const LevelP & level,
    SchedulerP & sched
)
{
    this->m_solver->scheduleSolve ( level, sched, this->m_materialManager->allMaterials(),
                                    matrix_label, Task::NewDW, // A
                                    u_label, false,            // x
                                    rhs_label, Task::NewDW,    // b
                                    u_label, Task::OldDW       // guess
                                  );
}
#endif

template < VarType VAR, StnType STN, bool AMR >
template < bool MG >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_postprocess (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_time_advance_postprocess", this, &Benchmark01<VAR, STN, AMR>::task_time_advance_postprocess );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, u_label, FGT, FGN );
    task->computes ( u0_label );
    task->computes ( energy_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template < VarType VAR, StnType STN, bool AMR >
template < bool MG >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleTimeAdvance_postprocess (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( !level->hasCoarserLevel() )
        scheduleTimeAdvance_postprocess < !MG > ( level, sched );
    else
    {
        Task * task = scinew Task ( "Benchmark01::task_time_advance_postprocess", this, &Benchmark01<VAR, STN, AMR>::task_time_advance_postprocess );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::NewDW, u_label, FGT, FGN );
        task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->computes ( u0_label );
        task->computes ( energy_label );
        sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
    }
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleRefine
(
    const PatchSet * new_patches,
    SchedulerP & sched
)
{
    const Level * level = getLevel ( new_patches );

    // no need to refine on coarser level
    if ( level->hasCoarserLevel() )
        scheduleRefine_solution ( new_patches, sched );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleRefine_solution (
    const PatchSet * new_patches,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_refine_solution", this, &Benchmark01::task_refine_solution );
    task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->computes ( u_label );

#ifdef HAVE_HYPRE
    /*************************** WORKAROUND ***************************/
    /* on new patches of finer level need to create matrix variables  */
    task->computes ( matrix_label );
    /************************* END WORKAROUND *************************/
#endif

    sched->addTask ( task, new_patches, this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleRefineInterface (
    const LevelP & /*level_fine*/,
    SchedulerP & /*sched*/,
    bool /*need_old_coarse*/,
    bool /*need_new_coarse*/
)
{
};

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleCoarsen
(
    const LevelP & level_coarse,
    SchedulerP & sched
)
{
    scheduleCoarsen_solution ( level_coarse, sched );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleCoarsen_solution (
    const LevelP & level_coarse,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_coarsen_solution", this, &Benchmark01::task_coarsen_solution );
    task->requires ( Task::NewDW, u_label, nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );
    task->modifies ( u_label );
    sched->addTask ( task, level_coarse->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleErrorEstimate
(
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleErrorEstimate_solution<AMR> ( level, sched );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleErrorEstimate_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "Benchmark01::task_error_estimate_solution", this, &Benchmark01::task_error_estimate_solution );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::NewDW, u_label, FGT, FGN );
    task->modifies ( this->m_regridder->getRefineFlagLabel(), this->m_regridder->refineFlagMaterials() );
    task->modifies ( this->m_regridder->getRefinePatchFlagLabel(), this->m_regridder->refineFlagMaterials() );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
Benchmark01<VAR, STN, AMR>::scheduleErrorEstimate_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( !level->hasCoarserLevel() ) scheduleErrorEstimate_solution < !MG > ( level, sched );
    else
    {
        Task * task = scinew Task ( "Benchmark01::task_error_estimate_solution", this, &Benchmark01::task_error_estimate_solution );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), FGT, FGN );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::NewDW, u_label, FGT, FGN );
        task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->modifies ( this->m_regridder->getRefineFlagLabel(), this->m_regridder->refineFlagMaterials() );
        task->modifies ( this->m_regridder->getRefinePatchFlagLabel(), this->m_regridder->refineFlagMaterials() );
        sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
    }
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::scheduleInitialErrorEstimate
(
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleErrorEstimate ( level, sched );
}

// TASKS

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::task_initialize_solution (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_initialize_solution ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch );;

        BlockRange range ( this->get_range ( patch ) );
        DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over range " << range );;

        DWView < ScalarField<double>, VAR, DIM > u ( dw_new, u_label, material, patch );
        parallel_for ( range, [patch, &u, this] ( int i, int j, int k )->void { initialize_solution ( {i, j, k}, patch, u ); } );
    }

    DOUT ( this->m_dbg_lvl2, myrank );;
}

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::task_compute_stable_timestep (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_compute_stable_timestep ====" );;
    dw_new->put ( delt_vartype ( delt ), this->getDelTLabel(), getLevel ( patches ) );
    DOUT ( this->m_dbg_lvl2, myrank );;
}


template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_forward_euler (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_time_advance_solution_forward_euler ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch );;

        DWFDView < ScalarField<const double>, STN, VAR > u_old ( dw_old, u_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > u_new ( dw_new, u_label, material, patch );

        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );
        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            parallel_for ( p.get_range(), [&u_old, &u_new, this] ( int i, int j, int k )->void { time_advance_solution_forward_euler ( {i, j, k}, u_old, u_new ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );;
}

#ifdef HAVE_HYPRE
template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypre
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    timeStep_vartype timeStepVar;
    dw_old->get ( timeStepVar, VarLabel::find ( timeStep_name ) );
    double timeStep = timeStepVar;

    if ( timeStep == 1 || this->isRegridTimeStep() )
        task_time_advance_solution_semi_implicit_assemble_hypre_full<SI> ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_semi_implicit_assemble_hypre_rhs<SI> ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypre_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );

        DWView < ScalarField<Stencil7>, VAR, DIM > A ( dw_new, matrix_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > b ( dw_new, rhs_label, material, patch );
        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            parallel_for ( p.get_range(), [&u_old, &A, &b, this] ( int i, int j, int k )->void { time_advance_solution_semi_implicit_assemble_hypre_full<SI> ( {i, j, k}, u_old, A, b ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypre_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypre_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label, patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );

        DWView < ScalarField<double>, VAR, DIM > b ( dw_new, rhs_label, material, patch );
        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            parallel_for ( p.get_range(), [&u_old, &b, this] ( int i, int j, int k )->void { time_advance_solution_semi_implicit_assemble_hypre_rhs<SI> ( {i, j, k}, u_old, b ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypresstruct
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    timeStep_vartype timeStepVar;
    dw_old->get ( timeStepVar, VarLabel::find ( timeStep_name ) );
    double timeStep = timeStepVar;

    if ( timeStep == 1 || this->isRegridTimeStep() )
        task_time_advance_solution_semi_implicit_assemble_hypresstruct_full<SI> ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_semi_implicit_assemble_hypresstruct_rhs<SI> ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypresstruct_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * /*matls*/,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );

        DWView < ScalarField<Stencil7>, VAR, DIM > A_stencil ( dw_new, matrix_label, material, patch );
        HypreSStruct::AdditionalEntries * A_additional = scinew HypreSStruct::AdditionalEntries;
        DWView < ScalarField<double>, VAR, DIM > b ( dw_new, rhs_label, material, patch );
        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            parallel_for ( p.get_range(), [&u_old, &A_stencil, &A_additional, &b, this] ( int i, int j, int k )->void { time_advance_solution_semi_implicit_assemble_hypresstruct_full<SI> ( {i, j, k}, u_old, A_stencil, A_additional, b ); } );
        }

        _AdditionalEntries additional_entries;
        additional_entries.setData ( A_additional );
        dw_new->put ( additional_entries, additional_entries_label, material, patch );
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, StnType STN, bool AMR>
template<TS SI>
void
Benchmark01<VAR, STN, AMR>::task_time_advance_solution_semi_implicit_assemble_hypresstruct_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_time_advance_solution_semi_implicit_assemble_hypresstruct_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label, patches, matls );
    dw_new->transferFrom ( dw_old, additional_entries_label, patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );

        DWView < ScalarField<double>, VAR, DIM > b ( dw_new, rhs_label, material, patch );
        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            parallel_for ( p.get_range(), [&u_old, &b, this] ( int i, int j, int k )->void { time_advance_solution_semi_implicit_assemble_hypresstruct_rhs<SI> ( {i, j, k}, u_old, b ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}
#endif // HAVE_HYPRE

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::task_time_advance_postprocess (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01<VAR,STN>::task_time_advance_postprocess ====" );;

    double energy = 0.;

    IntVector i0;
    Point p0 {M_PI, M_PI, 0.};
    double u0 = 0;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        const Level * level ( patch->getLevel() );

        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch );;

        bool has0 = ( patch->findCell ( p0, i0 ) ) &&
                    ( !level->hasFinerLevel() ||
                      !level->getFinerLevel()->containsCell ( level->mapCellToFiner ( i0 ) )
                    );

        SubProblems < ScalarProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );
        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & u = p.template get_fd_view<U> ( dw_new );
            parallel_reduce_sum (
                p.get_range(),
                [patch, &u, this] ( int i, int j, int k, double & energy )->void { time_advance_postprocess_energy ( {i, j, k}, patch, u, energy ); },
                energy
            );

            if ( has0 && p.get_codim() == 0 ) // on interal problems look for center
            {
                Point p ( DWInterface<VAR, DIM>::get_position ( level, i0 ) );
                Vector dist = ( p.asVector() - p0.asVector() ) / level->dCell();
                double w[2][2] = {{ 1., 1. }, { 1., 1. }};
                IntVector n[2][2] = {{ i0, i0 }, { i0, i0 }};
                const double & dx = dist[X];
                const double & dy = dist[Y];
                if ( dx < 0. )
                {
                    n[0][0][X] = n[0][1][X] -= 1;
                    w[0][0] *= -dx;
                    w[0][1] *= -dx;
                    w[1][0] *= 1 + dx;
                    w[1][1] *= 1 + dx;
                }
                else if ( dx > 0. )
                {
                    n[1][0][X] = n[1][1][X] += 1;
                    w[0][0] *= 1 - dx;
                    w[0][1] *= 1 - dx;
                    w[1][0] *= dx;
                    w[1][1] *= dx;
                }
                else
                {
                    w[1][0] = 0.;
                    w[1][1] = 0.;
                }

                if ( dy < 0. )
                {
                    n[0][0][Y] = n[1][0][Y] -= 1;
                    w[0][0] *= -dy;
                    w[1][0] *= -dy;
                    w[0][1] *= 1 + dy;
                    w[1][1] *= 1 + dy;
                }
                else if ( dy > 0. )
                {
                    n[0][1][Y] = n[1][1][Y] += 1;
                    w[0][0] *= 1 - dy;
                    w[1][0] *= 1 - dy;
                    w[0][1] *= dy;
                    w[1][1] *= dy;
                }
                else
                {
                    w[0][1] = 0.;
                    w[1][1] = 0.;
                }

                u0 = w[0][0] * u[ n[0][0] ] +
                     w[0][1] * u[ n[0][1] ] +
                     w[1][0] * u[ n[1][0] ] +
                     w[1][1] * u[ n[1][1] ];

                if ( fabs ( u0 ) > 2 )
                    SCI_THROW ( AssertionFailed ( "\n ERROR: Unstable simulation\n", __FILE__, __LINE__ ) );
            }
        }
    }

    dw_new->put ( sum_vartype ( energy ), energy_label );
    dw_new->put ( sum_vartype ( u0 ), u0_label );

    DOUT ( this->m_dbg_lvl2, myrank );;
}


template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::task_refine_solution
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches_fine,
    const MaterialSubset * /*matls*/,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_refine_solution ====" );

    for ( int p = 0; p < patches_fine->size(); ++p )
    {
        const Patch * patch_fine = patches_fine->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Fine Patch: " << *patch_fine << " Level: " << patch_fine->getLevel()->getIndex() );

        DWView < ScalarField<double>, VAR, DIM > u_fine ( dw_new, u_label, material, patch_fine );

        AMRInterpolator < ScalarProblem<VAR, STN>, U, C2F > u_coarse_interp ( dw_new, u_label, this->getSubProblemsLabel(), material, patch_fine );

        BlockRange range_fine ( this->get_range ( patch_fine ) );
        DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over fine range" << range_fine );
        parallel_for ( range_fine, [&u_coarse_interp, &u_fine, this] ( int i, int j, int k )->void { refine_solution ( {i, j, k}, u_coarse_interp, u_fine ); } );
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::task_coarsen_solution (
    const ProcessorGroup * myworld,
    const PatchSubset * patches_coarse,
    const MaterialSubset * /*matls*/,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_coarsen_solution " );

    for ( int p = 0; p < patches_coarse->size(); ++p )
    {
        const Patch * patch_coarse = patches_coarse->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Coarse Patch: " << *patch_coarse << " Level: " << patch_coarse->getLevel()->getIndex() );

        DWView < ScalarField<double>, VAR, DIM > u_coarse ( dw_new, u_label, material, patch_coarse );

        AMRRestrictor < ScalarProblem<VAR, STN>, U, F2C > u_fine_restr ( dw_new, u_label, this->getSubProblemsLabel(), material, patch_coarse, false );

        for ( const auto & region : u_fine_restr.get_support() )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over coarse cells region " << region );
            BlockRange range_coarse (
                Max ( region.getLow(), this->get_low ( patch_coarse ) ),
                Min ( region.getHigh(), this->get_high ( patch_coarse ) )
            );

            parallel_for ( range_coarse, [&u_fine_restr, &u_coarse, this] ( int i, int j, int k )->void { coarsen_solution ( {i, j, k}, u_fine_restr, u_coarse ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::task_error_estimate_solution
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * /*matls*/,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== Benchmark01::task_error_estimate_solution " );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );

        Variable<PP, PatchFlag> refine_patch_flag;
        dw_new->get ( refine_patch_flag, this->m_regridder->getRefinePatchFlagLabel(), material, patch );

        PatchFlag * patch_flag_refine = refine_patch_flag.get().get_rep();

        bool refine_patch = false;

        DWView < ScalarField<int>, CC, DIM > refine_flag ( dw_new, this->m_regridder->getRefineFlagLabel(), material, patch );
        SubProblems< ScalarProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );
            FDView < ScalarField<const double>, STN > & u = p.template get_fd_view<U> ( dw_new );
            parallel_reduce_sum ( p.get_range(), [&u, &refine_flag, &refine_patch, this] ( int i, int j, int k, bool & refine_patch )->void { error_estimate_solution<VAR> ( {i, j, k}, u, refine_flag, refine_patch ); }, refine_patch );
        }

        if ( refine_patch )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Setting refine flag" );
            patch_flag_refine->set();
        }
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

// IMPLEMENTATIONS

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::initialize_solution (
    const IntVector & id,
    Patch const * patch,
    View< ScalarField<double> > & u
)
{
    Vector v ( this->get_position ( patch, id ).asVector() );

#if defined(__INTEL_COMPILER) && defined(BUG_WORKAROUND)
    // BUG workaround
    std::stringstream ss;
    ss << v << std::endl;
#endif

    v[0] -= M_PI;
    v[1] -= M_PI;
    v[2] = 0.;
    u[id] = tanh ( ( v.length() - 2. ) / ( epsilon * M_SQRT2 ) );
}

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::time_advance_solution_forward_euler (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & u_old,
    View< ScalarField<double> > & u_new
)
{
    const double & u = u_old[id];
    double lap = u_old.laplacian ( id );
    double src = u * ( u * u - 1. );
    double delta_u = delt * ( epsilon * epsilon * lap - src );
    u_new[id] = u + delta_u;
}

#ifdef HAVE_HYPRE
template < VarType VAR, StnType STN, bool AMR >
template<TS SI>
typename std::enable_if < SI == TS::SemiImplicit0, void >::type
Benchmark01<VAR, STN, AMR>::time_advance_solution_semi_implicit_assemble_hypre_full (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & u_old,
    View < ScalarField<Stencil7> > & A,
    View < ScalarField<double> > & b
)
{
    std::tuple<Stencil7, double> sys = u_old.laplacian_sys_hypre ( id );

    const Stencil7 & lap_stn = std::get<0> ( sys );
    const double & rhs = std::get<1> ( sys );
    const double & u = u_old[id];
    const double a = epsilon * epsilon * delt;
    const double s = ( 1 + delt - delt * u * u );

    for ( int i = 0; i < 7; ++i )
        A[id][i] = -a * lap_stn[i];
    A[id].p += 1;
    b[id] = s * u + a * rhs;
}

template < VarType VAR, StnType STN, bool AMR >
template<TS SI>
typename std::enable_if < SI == TS::SemiImplicit0, void >::type
Benchmark01<VAR, STN, AMR>::time_advance_solution_semi_implicit_assemble_hypre_rhs (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & u_old,
    View < ScalarField<double> > & b
)
{
    double rhs = u_old.laplacian_rhs_hypre ( id );
    const double & u = u_old[id];
    const double a = epsilon * epsilon * delt;
    const double s = ( 1 + delt - delt * u * u );

    b[id] = s * u + a * rhs;
}

template < VarType VAR, StnType STN, bool AMR >
template<TS SI>
typename std::enable_if < SI == TS::SemiImplicit0, void >::type
Benchmark01<VAR, STN, AMR>::time_advance_solution_semi_implicit_assemble_hypresstruct_full (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & u_old,
    View < ScalarField<Stencil7> > & A_stencil,
    HypreSStruct::AdditionalEntries * A_additional,
    View < ScalarField<double> > & b
)
{
    std::tuple<Stencil7, HypreSStruct::AdditionalEntries, double> sys = u_old.laplacian_sys_hypresstruct ( id );

    const Stencil7 & lap_stn = std::get<0> ( sys );
    HypreSStruct::AdditionalEntries & lap_extra = std::get<1> ( sys );
    const double & rhs = std::get<2> ( sys );
    const double & u = u_old[id];
    const double a = epsilon * epsilon * delt;
    const double s = ( 1 + delt - delt * u * u );

    for ( int i = 0; i < 7; ++i )
        A_stencil[id][i] = -a * lap_stn[i];
    A_stencil[id].p += 1;
    for ( auto & entry : lap_extra )
    {
        auto it = A_additional->find ( entry.first );
        if ( it != A_additional->end() ) it->second -= a * entry.second;
        else A_additional->emplace ( entry.first, -a * entry.second );
    }
    b[id] = s * u + a * rhs;
}

template < VarType VAR, StnType STN, bool AMR >
template<TS SI>
typename std::enable_if < SI == TS::SemiImplicit0, void >::type
Benchmark01<VAR, STN, AMR>::time_advance_solution_semi_implicit_assemble_hypresstruct_rhs (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & u_old,
    View < ScalarField<double> > & b
)
{
    double rhs = u_old.laplacian_rhs_hypresstruct ( id );
    const double & u = u_old[id];
    const double a = epsilon * epsilon * delt;
    const double s = ( 1 + delt - delt * u * u );

    b[id] = s * u + a * rhs;
}
#endif

template < VarType VAR, StnType STN, bool AMR >
void
Benchmark01<VAR, STN, AMR>::time_advance_postprocess_energy (
    const IntVector & id,
    const Patch * patch,
    const FDView < ScalarField<const double>, STN > & u_new,
    double & energy
)
{
    const Level * level = patch->getLevel();
    if ( !level->hasFinerLevel() || !level->getFinerLevel()->containsCell ( level->mapCellToFiner ( id ) ) )
    {
        const double & u = u_new[id];
        auto grad = u_new.gradient ( id );
        double A = level->dCell() [0] * level->dCell() [1];
        energy += A * ( epsilon * epsilon * ( grad[0] * grad[0] + grad[1] * grad[1] ) / 2. + ( u * u * u * u - 2 * u * u + 1. ) / 4. );
    }
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::refine_solution
(
    const IntVector id_fine,
    const View < ScalarField<const double> > & u_coarse_interp,
    View < ScalarField<double> > & u_fine
)
{
    u_fine[id_fine] = u_coarse_interp[id_fine];
}

template<VarType VAR, StnType STN, bool AMR>
void
Benchmark01<VAR, STN, AMR>::coarsen_solution
(
    const IntVector id_coarse,
    const View < ScalarField<const double> > & u_fine_restr,
    View < ScalarField<double> > & u_coarse
)
{
    u_coarse[id_coarse] = u_fine_restr[id_coarse];
}

template<VarType VAR, StnType STN, bool AMR>
template<VarType V>
typename std::enable_if < V == CC, void >::type
Benchmark01<VAR, STN, AMR>::error_estimate_solution
(
    const IntVector id,
    FDView < ScalarField<const double>, STN > & u,
    View < ScalarField<int> > & refine_flag,
    bool & refine_patch
)
{
    bool refine = false;
    auto grad = u.gradient ( id );
    double err2 = 0;
    for ( size_t d = 0; d < DIM; ++d )
        err2 += grad[d] * grad[d];
    refine = err2 > refine_threshold * refine_threshold;
    refine_flag[id] = refine;
    refine_patch |= refine;
}

template<VarType VAR, StnType STN, bool AMR>
template<VarType V>
typename std::enable_if < V == NC, void >::type
Benchmark01<VAR, STN, AMR>::error_estimate_solution
(
    const IntVector id,
    FDView < ScalarField<const double>, STN > & u,
    View < ScalarField<int> > & refine_flag,
    bool & refine_patch
)
{
    auto grad = u.gradient ( id );
    double err2 = 0;
    for ( size_t d = 0; d < DIM; ++d )
        err2 += grad[d] * grad[d];
    if ( err2 > refine_threshold * refine_threshold )
    {
        refine_patch = true;

        // loop over all cells sharing node id
        IntVector id0 = id - get_dim<DIM>::unit_vector();
        IntVector i;
        for ( i[Z] = id0[Z]; i[Z] <= id[Z]; ++i[Z] )
            for ( i[Y] = id0[Y]; i[Y] <= id[Y]; ++i[Y] )
                for ( i[X] = id0[X]; i[X] <= id[X]; ++i[X] )
                    if ( refine_flag.is_defined_at ( i ) )
                        refine_flag[i] = 1;
    }
}

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_Benchmark01_h
