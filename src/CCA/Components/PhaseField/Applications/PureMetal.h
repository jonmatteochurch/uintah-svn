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
 * @file CCA/Components/PhaseField/Applications/PureMetal.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Applications_PureMetal_h
#define Packages_Uintah_CCA_Components_PhaseField_Applications_PureMetal_h

#include <sci_defs/hypre_defs.h>

#ifdef HAVE_HYPRE
#   include <CCA/Components/Solvers/HypreSStruct/AdditionalEntriesP.h>
#   include <CCA/Components/Solvers/HypreSStruct/Solver.h>
#endif

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/Util/Expressions.h>
#include <CCA/Components/PhaseField/DataTypes/PureMetalProblem.h>
#include <CCA/Components/PhaseField/DataTypes/SubProblems.h>
#include <CCA/Components/PhaseField/DataTypes/ScalarField.h>
#include <CCA/Components/PhaseField/DataTypes/VectorField.h>
#include <CCA/Components/PhaseField/Applications/Application.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessModule.h>
#include <CCA/Components/PhaseField/Views/View.h>
#include <CCA/Components/PhaseField/Views/FDView.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWView.h>
#include <CCA/Components/PhaseField/AMR/AMRInterpolator.h>
#include <CCA/Components/PhaseField/AMR/AMRRestrictor.h>

#include <CCA/Ports/Regridder.h>

#include <Core/Util/Factory/Implementation.h>
#include <Core/Util/DebugStream.h>
#include <Core/Grid/SimpleMaterial.h>
#include <Core/Grid/Variables/PerPatchVars.h>
#include <Core/Grid/Variables/ReductionVariable.h>

namespace Uintah
{
namespace PhaseField
{

/// Debugging stream for component schedulings
static constexpr bool dbg_pure_metal_scheduling = false;

/**
 * @brief PureMetal PhaseField applications
 *
 * Implements a Finite Difference solver for the simulation of the anisotropic
 * solidification of an under-cooled metal around a solid seed using the model
 * from
 *
 * A, Karma and W.-J. Rappel,
 * "Phase-Field Method for Computationally Efficient Modeling of Solidification
 * with Arbitrary Interface Kinetics",
 * Phys.Rev.E, 1996.
 *
 * Phase-field equation for \f$\psi : \Omega \to [-1. 1]\f$
 * \f[
 * \tau_0 A^2 \dot \psi = \psi (1-\psi^2) - \lambda u (\psi^2 - 1)^2
 *                      + W_0^2 A^2 \nabla^2 \psi
 *                      + W_0^2 (A^2_x - \partial_y B_{xy} - \partial_y B_{xz}) \psi_x
 *                      + W_0^2 (A^2_y + \partial_y B_{xy} - \partial_y B_{yz}) \psi_y
 *                      + W_0^2 (A^2_z + \partial_y B_{xz} + \partial_y B_{yz}) \psi_z
 * \f]
 * (for 2D problems all quantities with \f$z\f$ in the subscript are null)
 *
 * non-dimensional temperature equation for \f$u:\Omega \to \mathbb R\f$
 * \f[
 * \dot u = \alpha \nabla^2 u + \frac 12 \dot \psi
 * \f]
 *
 * with the following anisotropy functions
 * \f[ A    = 1 - 3 \epsilon + 4\epsilon \sum n_i^4 \f]
 * \f[ B_{ij} = 16 \epsilon A \psi_i \psi_j (\psi_i^2-\psi_j^2) / |\nabla\psi|^4 \f]
 * where \f$ n=\frac{\nabla\psi}{|\nabla\psi|} \f$
 *
 * The model parameters are:
 * - \f$ \lambda \f$  coupling parameter
 * - \f$ W_0 \f$      width scaling (\f$ W = W_0 A \f$ is the interface width)
 * - \f$ \tau_0 \f$   time scaling (\f$ \tau = \tau_0 A^2 \f$ is the characteristic
 *                    time of attachment of atoms at the interface)
 * - \f$ \epsilon \f$ anisotropy strength
 * - \f$ \alpha \f$   thermal diffusivity
 *
 * The following non-dimensionalizations are performed
 * - \f$ x \mapsto x/W_0 \f$
 * - \f$ t \mapsto t/\tau_0 \f$
 * - \f$ \alpha \mapsto \alpha\tau_0/W_0^2 \f$
 *
 * @tparam VAR type of variable representation
 * @tparam DIM problem dimensions
 * @tparam STN finite-difference stencil
 * @tparam AMR whether to use adaptive mesh refinement
 */
template<VarType VAR, DimType DIM, StnType STN, bool AMR = false>
class PureMetal
    : public Application< PureMetalProblem<VAR, STN>, AMR >
    , public Implementation < PureMetal<VAR, DIM, STN, AMR>, UintahParallelComponent, const ProcessorGroup *, const MaterialManagerP, int>
{
public:

    /// Problem material index (only one SimpleMaterial)
    static constexpr int material = 0;

    /// Number of anisotropy functions B
    static constexpr size_t BSZ = combinations<DIM, 2>::value;

    /// Index for anisotropy functions B
    static constexpr size_t XY = 0;

    /// Index for anisotropy functions B
    static constexpr size_t XZ = 1;

    /// Index for anisotropy functions B
    static constexpr size_t YZ = 2;

    /// Index of variables within PureMetalProblem
    static constexpr size_t PSI = 0; ///< Index for phase-field
    static constexpr size_t U = 1;   ///< Index for non-dimensional temperature field
    static constexpr size_t A2 = 2;  ///< Index for the square of the anisotropy function
    static constexpr size_t B = 3;   ///< Index for the anisotropy terms \f$ B_ij \f$

    static constexpr int NVARS = 2;

private:

    /// Number of ghost elements required by STN (on the same level)
    using Application< PureMetalProblem<VAR, STN> >::FGN;

    /// Type of ghost elements required by VAR and STN (on coarser level)
    using Application< PureMetalProblem<VAR, STN> >::FGT;

    /// Number of ghost elements required by STN (on coarser level)
    using Application< PureMetalProblem<VAR, STN> >::CGN;

    /// Type of ghost elements required by VAR and STN (on the same level)
    using Application< PureMetalProblem<VAR, STN> >::CGT;

    /// Interpolation type for refinement
    using Application< PureMetalProblem<VAR, STN> >::C2F;

    /// Restriction type for coarsening
    using Application< PureMetalProblem<VAR, STN> >::F2C;

#ifdef HAVE_HYPRE
    /// Non-stencil entries type for implicit matrix
    using _AdditionalEntries = PerPatch<HypreSStruct::AdditionalEntriesP>;
#endif

    /// If grad_psi_norm2 is less than tol than psi is considered constant when computing anisotropy terms
    static constexpr double tol = 1.e-6;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    const static FactoryString Name;

protected: // MEMBERS

    /// Label for phase field in the DataWarehouse
    const VarLabel * psi_label;

    /// Label for non dimensional temperature field in the DataWarehouse
    const VarLabel * u_label;

    /// Label for the norm of the phase field gradient in the DataWarehouse
    const VarLabel * grad_psi_norm2_label;

    /// Label for anisotropy field A in the DataWarehouse
    const VarLabel * a_label;

    /// Label for anisotropy field A^2 in the DataWarehouse
    const VarLabel * a2_label;

    /// Label for the the phase field gradient in the DataWarehouse
    std::array<const VarLabel *, DIM> grad_psi_label;

    /// Label for anisotropy fields B in the DataWarehouse
    std::array<const VarLabel *, BSZ> b_label;

#ifdef HAVE_HYPRE
    /// Label for the implicit matrix stencil entries in the DataWarehouse
    const VarLabel *** matrix_label;

    /// Label for the implicit vector in the DataWarehouse
    const VarLabel ** rhs_label;

    /// Label for the implicit matrix non-stencil entries in the DataWarehouse
    const VarLabel ** additional_entries_label;
#endif // HAVE_HYPRE

    /// Time step size
    double delt;

    /// Coupling parameter
    double lambda;

    /// Non-dimensional thermal diffusivity
    double alpha;

    /// Anisotropy strength (\f$ \epsilon<0 \f$ favours growth along I bisector instead of along \f$ x \f$ axis
    double epsilon;

    /// Initial phase field interface width
    double gamma_psi;

    /// Initial temperature field interface width
    double gamma_u;

    /// Initial seed radius
    double r0;

    /// Initial undercooling
    double delta;

    /// Threshold for AMR
    double refine_threshold;

#ifdef HAVE_HYPRE
    /// Time advance scheme
    TS time_scheme;
#endif

    /// Module for post-processing tip info
    ArmPostProcessModule<VAR, DIM, STN, AMR> * post_process;

    bool is_time_advance_anisotropy_terms_scheduled;
    bool is_time_advance_solution_scheduled;

public: // CONSTRUCTORS/DESTRUCTOR

    /**
     * @brief Constructor
     *
     * Intantiate a PureMetal application
     *
     * @param myWorld data structure to manage mpi processes
     * @param materialManager data structure to manage materials
     * @param verbosity constrols amount of debugging output
     */
    PureMetal (
        const ProcessorGroup * myWorld,
        const MaterialManagerP materialManager,
        int verbosity = 0
    );

    /**
     * @brief Destructor
     */
    virtual ~PureMetal();

    /// Prevent copy (and move) constructor
    PureMetal ( const PureMetal & ) = delete;

    /// Prevent copy (and move) assignment
    /// @return deleted
    PureMetal & operator= ( const PureMetal & ) = delete;

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
     * solution derivatives allowing sched to control its execution order
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
     * @brief Schedule task_initialize_grad_psi (non AMR implementation)
     *
     * Defines the dependencies and output of the task which initializes psi
     * and grad_psi_norm2 allowing sched to control its execution order.
     * Initialize also the anisotropy terms to 0
     *
     * @param level grid level to be initialized
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleInitialize_grad_psi (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_initialize_grad_psi (AMR implementation)
     *
     * Does nothing since grad_psi and grad_psi_norm2 are computed by
     * error_estimate_grad_psi.
     *
     * @param level unused
     * @param sched unused
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleInitialize_grad_psi (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule the initialization tasks for restarting a simulation
     *
     * Specify all tasks to be performed at fist timestep after a stop to
     * initialize not saved variables to the DataWarehouse
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
     * @brief Schedule task_time_advance_grad_psi (non AMR implementation)
     *
     * Defines the dependencies and output of the task which updates grad psi
     * and grad_psi_norm2 allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_grad_psi (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_grad_psi (AMR implementation)
     *
     * Defines the dependencies and output of the task which updates grad psi
     * and grad_psi_norm2 allowing sched to control its execution order
     *
     * @remark does nothing since it grad_psi and grad_psi_norm2 are already
     * computed by error_estimate_grad_psi
     *
     * @param level unused
     * @param sched unused
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_grad_psi (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_anisotropy_terms
     *
     * Defines the dependencies and output of the task which updates the
     * anisotropy terms allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     *
     * @remark since no derivative is involved there is no need to have two
     * different implementation for non AMR and non AMR cases
     */
    void
    scheduleTimeAdvance_anisotropy_terms (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution
     *
     * Switches between available implementations depending on the given time
     * scheme
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
     * (non AMR implementation)
     *
     * Defines the dependencies and output of the task which updates psi
     * and u allowing sched to control its execution order
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
     * (AMR implementation)
     *
     * Defines the dependencies and output of the task which updates psi
     * and u allowing sched to control its execution order
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
     * @brief Schedule task_time_advance_solution_implicit
     * (non AMR implementation)
     *
     * Switches between available implementations depending on the given
     * time solver
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_implicit (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_implicit
     * (AMR implementation)
     *
     * Switches between available implementations depending on the given
     * time solver
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_implicit (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_backward_euler_assemble_hypre
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles
     * the implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_backward_euler_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_backward_euler_assemble_hypre
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles
     * the implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_backward_euler_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_backward_euler_assemble_hypresstruct
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_backward_euler_assemble_hypresstruct
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_crank_nicolson_assemble_hypre
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_crank_nicolson_assemble_hypre
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_crank_nicolson_assemble_hypresstruct
     * (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_time_advance_solution_crank_nicolson_assemble_hypresstruct
     * (refinement level implementation)
     *
     * Defines the dependencies and output of the task which assembles the
     * implicit system allowing sched to control its execution order
     *
     * @param level grid level to be updated
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );

    void
    scheduleTimeAdvannce_comunicate_before_solve_hypresstruct (
        const LevelP & level,
        SchedulerP & sched
    );
#endif // HAVE_HYPRE

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
     * @brief Schedule task_communicate_psi
     *
     * Defines the dependencies and output of the task which forces the
     * communication of ghost layers around refine coarse cells since they are
     * not triggered by the dependencies scheduleRefine_grad_psi
     *
     * @param level refined coarse level
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleRefine_communicate_psi (
        const LevelP & level,
        SchedulerP & sched
    )
    {
        Task * task = scinew Task ( "PureMetal::task_communicate_psi", this, &PureMetal::task_empty );
        task->requires ( Task::NewDW, psi_label, CGT, CGN );
        task->modifies ( psi_label );
        sched->addTask ( task, sched->getLoadBalancer()->getPerProcessorPatchSet ( level ), this->m_materialManager->allMaterials() );
    }

    /**
     * @brief Schedule task_refine_grad_psi
     *
     * Defines the dependencies and output of the task which computes the
     * derivatives of psi on each one of the new_patches allowing sched to
     * control its execution order
     *
     * @param new_patches patches to be populated
     * @param sched scheduler to manage the tasks
     */
    void
    scheduleRefine_grad_psi (
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
     * @brief Schedule task_error_estimate_grad_psi (coarsest level implementation)
     *
     * Defines the dependencies and output of the task which estimates the
     * spatial discretization error allowing sched to control its execution order
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleErrorEstimate_grad_psi (
        const LevelP & level,
        SchedulerP & sched
    );

    /**
     * @brief Schedule task_error_estimate_grad_psi (refinement level implementation)
     *
     * Defines the dependencies and output of the task which estimates the
     * spatial discretization error allowing sched to control its execution order
     *
     * @param level level to check
     * @param sched scheduler to manage the tasks
     */
    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleErrorEstimate_grad_psi (
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
     * Allocate and save variables for psi and u for each one of the patches
     * and save them to dw_new
     * @remark initialize also anisotropy terms to 0
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

#ifdef HAVE_HYPRE
    /**
     * @brief Initialize hypre sstruct structures on restart task
     *
     * Allocate and save variables for non-stencil variables for each
     * one of the patches and save them to dw_new
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_restart_initialize_hypresstruct (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );
#endif // HAVE_HYPRE

    /**
     * @brief Initialize grad psi task
     *
     * Allocate and save variable for psi derivatives for each one of the patches
     * and save it to dw_new
     * @remark initialize also anisotropy terms to 0
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_initialize_grad_psi (
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
     * @brief Advance grad_psi task
     *
     * Computes the gradient of the phase field using its value at the previous
     * timestep
     *
     * @remark this task is scheduled only if AMR is disabled since the the
     * gradient of the phase field for AMR simulations is already computed by
     * task_error_estimate_grad_psi
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_grad_psi (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Advance anisotropy terms task
     *
     * Computes A, A2 and B using the value of grad_psi at the previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_anisotropy_terms (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Advance solution task (Forward Euler implementation)
     *
     * Computes new value of psi and u using the newly computed
     * anisotropy terms together with the value of the solution and
     * grad_psi at the previous timestep
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
     * @brief Assemble hypre system task (Backward Euler implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypre (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Backward Euler, full implementation)
     *
     * Assemble both implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypre_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Backward Euler, rhs implementation)
     *
     * Assemble only implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypre_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Backward Euler implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypresstruct (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Backward Euler, full implementation)
     *
     * Assemble both implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypresstruct_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Backward Euler, rhs implementation)
     *
     * Assemble only implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_backward_euler_assemble_hypresstruct_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Crank Nicolson implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypre (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Crank Nicolson, full implementation)
     *
     * Assemble both implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypre_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre system task (Crank Nicolson, rhs implementation)
     *
     * Assemble only implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypre_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Crank Nicolson implementation)
     *
     * Switches between available implementations to avoid assembling the matrix
     * when not required
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypresstruct (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Crank Nicolson, full implementation)
     *
     * Assemble both implicit matrix and vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypresstruct_full (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Assemble hypre sstruct system task (Crank Nicolson, rhs implementation)
     *
     * Assemble only implicit vector for hypre using the value of the
     * solution and at previous timestep
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    void
    task_empty (
        const ProcessorGroup * /*myworld*/,
        const PatchSubset * /*patches*/,
        const MaterialSubset * /*matls*/,
        DataWarehouse * /*dw_old*/,
        DataWarehouse * /*dw_new*/
    ) {}
#endif // HAVE_HYPRE

    /**
     * @brief Refine solution task
     *
     * Computes interpolated value of u and psi on new refined patched
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
     * @brief Refine grad psi task
     *
     * Computes value of grad psi on new refined patched from the interpolated
     * values of psi
     *
     * @param myworld data structure to manage mpi processes
     * @param patches_fine list of patches to be initialized
     * @param matls unused
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_refine_grad_psi (
        const ProcessorGroup * myworld,
        const PatchSubset * patches_fine,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

    /**
     * @brief Coarsen solution task
     *
     * Restricted value of u and psi from refined regions to coarse patches
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
     * @brief ErrorEstimate grad_psi task
     *
     * Computes the gradient of the phase field using its value at the previous
     * timestep and set refinement flag where it is above the threshold given
     * in input
     *
     * @remark this task replaces task_time_advance_grad_psi for AMR simulations
     *
     * @param myworld data structure to manage mpi processes
     * @param patches list of patches to be initialized
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_error_estimate_grad_psi (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

protected: // IMPLEMENTATIONS

    template <int D = DIM>
    typename std::enable_if<D == D2, double>::type
    anisotropy_term
    (
        const IntVector & id,
        const View < VectorField<const double, DIM> > & grad_psi,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b
    )
    {
        return ( ( a2.dx ( id ) - b[XY].dy ( id ) ) * grad_psi[X][id]
                 + ( a2.dy ( id ) + b[XY].dx ( id ) ) * grad_psi[Y][id] );
    }

    template <int D = DIM>
    typename std::enable_if<D == D3, double>::type
    anisotropy_term
    (
        const IntVector & id,
        const View < VectorField<const double, DIM> > & grad_psi,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b
    )
    {
        return ( ( a2.dx ( id ) - b[XY].dy ( id ) - b[XZ].dz ( id ) ) * grad_psi[X][id]
                 + ( a2.dy ( id ) + b[XY].dx ( id ) - b[YZ].dz ( id ) ) * grad_psi[Y][id]
                 + ( a2.dz ( id ) + b[XZ].dx ( id ) + b[YZ].dy ( id ) ) * grad_psi[Z][id] );
    }

    /**
     * @brief Initialize solution implementation
     *
     * compute initial condition for psi and u at a given grid position
     *
     * @param id grid index
     * @param patch grid patch
     * @param[out] psi view of the phase field in the new dw
     * @param[out] u view of the temperature field in the new dw
     */
    void
    initialize_solution (
        const IntVector & id,
        const Patch * patch,
        View < ScalarField<double> > & psi,
        View < ScalarField<double> > & u
    );

    /**
     * @brief Advance grad_psi implementation
     *
     * compute new value for grad_psi at a given grid position using its value
     * at the previous timestep
     *
     * @param id grid index
     * @param psi view of the phase field in the old dw
     * @param[out] grad_psi view of the phase gradient field in the new dw
     * @param[out] grad_psi_norm2 view of the norm of the phase gradient field
     * in the new dw
     */
    void
    time_advance_grad_psi (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi,
        View < VectorField<double, DIM > > & grad_psi,
        View < ScalarField<double> > & grad_psi_norm2
    );

    /**
     * @brief Advance anisotropy terms implementation
     *
     * compute new value for grad_psi at a given grid position using its value
     * at the previous timestep
     *
     * @param id grid index
     * @param grad_psi view of the phase gradient field in the old dw
     * @param grad_psi_norm2 view of the norm of the phase gradient field
     * in the old dw
     * @param[out] a view of A in the new dw
     * @param[out] a2 view of A2 in the new dw
     * @param[out] b view of B in the new dw
     */
    void
    time_advance_anisotropy_terms_dflt (
        const IntVector & id,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & grad_psi_norm2,
        View < ScalarField<double> > & a,
        View < ScalarField<double> > & a2,
        View < VectorField<double, BSZ> > & b
    );

    void
    time_advance_anisotropy_terms_diag (
        const IntVector & id,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & grad_psi_norm2,
        View < ScalarField<double> > & a,
        View < ScalarField<double> > & a2,
        View < VectorField<double, BSZ> > & b
    );

    void ( PureMetal::*time_advance_anisotropy_terms ) (
        const IntVector &,
        const View < VectorField<const double, DIM> > &,
        const View < ScalarField<const double> > &,
        View < ScalarField<double> > &,
        View < ScalarField<double> > &,
        View < VectorField<double, BSZ> > &
    );

    /**
     * @brief Advance solution terms implementation
     *
     * compute new value for psi and u at a given grid position using the newly
     * computed anisotropy terms together with the value of the solution and
     * grad_psi at the previous timestep
     *
     * @param id grid index
     * @param psi_old view of the phase field in the old dw
     * @param u_old view of the temperature field in the old dw
     * @param grad_psi view of the phase gradient field in the old dw
     * @param a view of A in the new dw
     * @param a2 view of A2 in the new dw
     * @param b view of B in the new dw
     * @param[out] psi_new view of the phase field in the new dw
     * @param[out] u_new view of the temperature field in the new dw
     */
    void
    time_advance_solution_forward_euler (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<double> > & psi_new,
        View < ScalarField<double> > & u_new
    );

#ifdef HAVE_HYPRE
    /**
     * @brief Assemble hypre system implementation (Backward Euler, full implementation)
     *
     * compute both implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_backward_euler_assemble_hypre_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre system implementation (Backward Euler, rhs implementation)
     *
     * compute only implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_backward_euler_assemble_hypre_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre sstruct system implementation (Backward Euler, full implementation)
     *
     * compute both implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] A_additional view of the implicit matrix non-stencil field in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_backward_euler_assemble_hypresstruct_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
        HypreSStruct::AdditionalEntries * ( & A_additional ) [NVARS],
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre sstruct system implementation (Backward Euler, rhs implementation)
     *
     * compute only implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_backward_euler_assemble_hypresstruct_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre system implementation (Crank Nicolson, full implementation)
     *
     * compute both implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_crank_nicolson_assemble_hypre_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre system implementation (Crank Nicolson, rhs implementation)
     *
     * compute only implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_crank_nicolson_assemble_hypre_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre sstruct system implementation (Crank Nicolson, full implementation)
     *
     * compute both implicit matrix stencil and vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] A view of the implicit matrix stencil field in the new dw
     * @param[out] A_additional view of the implicit matrix non-stencil field in the new dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_crank_nicolson_assemble_hypresstruct_full (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
        HypreSStruct::AdditionalEntries * ( & A_additional ) [NVARS],
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );

    /**
     * @brief Assemble hypre sstruct system implementation (Crank Nicolson, rhs implementation)
     *
     * compute only implicit vector entries at a given grid
     * position using the value of the solution and at previous timestep
     *
     * @param id grid index
     * @param u_old view of the solution field in the old dw
     * @param[out] b view of the implicit vector field in the new dw
     */
    void
    time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs (
        const IntVector & id,
        const FDView < ScalarField<const double>, STN > & psi_old,
        const FDView < ScalarField<const double>, STN > & u_old,
        const View < VectorField<const double, DIM> > & grad_psi,
        const View < ScalarField<const double> > & a,
        const FDView < ScalarField<const double>, STN > & a2,
        const FDView < VectorField<const double, BSZ>, STN > & b,
        View < ScalarField<double> > * ( & rhs ) [NVARS]
    );
#endif // HAVE_HYPRE

    /**
     * @brief Refine solution implementation
     *
     * Computes interpolated value of u and psi at a given grid position

     * @param id_fine fine grid index
     * @param psi_coarse_interp interpolator of the phase field on the coarse level
     * @param u_coarse_interp interpolator of the temperature field on the coarse level
     * @param[out] psi_fine view of the phase field on the fine level
     * @param[out] u_fine view of the temperature field on the fine level
     */
    void
    refine_solution (
        const IntVector id_fine,
        const View < ScalarField<const double> > & psi_coarse_interp,
        const View < ScalarField<const double> > & u_coarse_interp,
        View < ScalarField<double> > & psi_fine,
        View < ScalarField<double> > & u_fine
    );

    /**
     * @brief Coarsen solution implementation
     *
     * Computes restricted value of u and psi at a given grid position

     * @param id_coarse coarse grid index
     * @param psi_fine_restr restrictor of the phase field on the fine level
     * @param u_fine_restr restrictor of the temperature field on the fine level
     * @param[out] psi_coarse view of the phase field on the coarse level
     * @param[out] u_coarse view of the temperature field on the coarse level
     */
    void
    coarsen_solution (
        const IntVector id_coarse,
        const View < ScalarField<const double> > & psi_fine_restr,
        const View < ScalarField<const double> > & u_fine_restr,
        View < ScalarField<double> > & psi_coarse,
        View < ScalarField<double> > & u_coarse
    );

    /**
     * @brief ErrorEstimate grad_psi implementation
     *
     * Computes the gradient of the phase field using its value at the previous
     * timestep and set refinement flag where it is above the threshold given
     * in input
     *
     * @param id grid index
     * @param psi view of the phase field in the old dw
     * @param[out] grad_psi view of the phase gradient field in the new dw
     * @param[out] grad_psi_norm2 view of the norm of the phase gradient field
     * in the new dw
     * @param[out] refine_flag view of refine flag (grid field) in the new dw
     */
    void
    error_estimate_grad_psi (
        const IntVector & id,
        FDView < ScalarField<const double>, STN > & psi,
        View < VectorField<double, DIM> > & grad_psi,
        View < ScalarField<double> > & grad_psi_norm2,
        View< ScalarField<int> > & refine_flag
    );

}; // class PureMetal

// CONSTRUCTORS/DESTRUCTOR

using namespace std::placeholders;

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
PureMetal<VAR, DIM, STN, AMR>::PureMetal (
    const ProcessorGroup * my_world,
    MaterialManagerP const material_manager,
    int verbosity
) : Application< PureMetalProblem<VAR, STN>, AMR > ( my_world, material_manager, verbosity ),
    post_process ( nullptr ),
    time_advance_anisotropy_terms ( &PureMetal::time_advance_anisotropy_terms_dflt )
{
    psi_label = VarLabel::create ( "psi", Variable<VAR, double>::getTypeDescription() );
    u_label = VarLabel::create ( "u", Variable<VAR, double>::getTypeDescription() );
    grad_psi_norm2_label = VarLabel::create ( "grad_psi_norm2", Variable<VAR, double>::getTypeDescription() );
    a_label = VarLabel::create ( "A", Variable<VAR, double>::getTypeDescription() );
    a2_label = VarLabel::create ( "A2", Variable<VAR, double>::getTypeDescription() );

    grad_psi_label[X] = VarLabel::create ( "psi_x", Variable<VAR, double>::getTypeDescription() );
    if ( DIM > D1 )
    {
        grad_psi_label[Y] = VarLabel::create ( "psi_y", Variable<VAR, double>::getTypeDescription() );
        b_label[XY] = VarLabel::create ( "Bxy", Variable<VAR, double>::getTypeDescription() );
    }
    if ( DIM > D2 )
    {
        grad_psi_label[Z] = VarLabel::create ( "psi_z", Variable<VAR, double>::getTypeDescription() );
        b_label[XZ] = VarLabel::create ( "Bxz", Variable<VAR, double>::getTypeDescription() );
        b_label[YZ] = VarLabel::create ( "Byz", Variable<VAR, double>::getTypeDescription() );
    }

#ifdef HAVE_HYPRE
    matrix_label = scinew const VarLabel ** [NVARS];
    for ( int i = 0; i < NVARS; ++i ) matrix_label[i] = scinew const VarLabel * [NVARS];
    additional_entries_label = scinew const VarLabel * [NVARS];
    rhs_label = scinew const VarLabel * [NVARS];
    matrix_label[PSI][PSI] = VarLabel::create ( "A_psipsi", Variable<VAR, Stencil7>::getTypeDescription() );
    matrix_label[PSI][U] = VarLabel::create ( "A_psiu", Variable<VAR, Stencil7>::getTypeDescription() );
    matrix_label[U][PSI] = VarLabel::create ( "A_upsi", Variable<VAR, Stencil7>::getTypeDescription() );
    matrix_label[U][U] = VarLabel::create ( "A_uu", Variable<VAR, Stencil7>::getTypeDescription() );
    additional_entries_label[PSI] = VarLabel::create ( "A_psipsi" + HypreSStruct::Solver<DIM>::AdditionalEntriesSuffix, _AdditionalEntries::getTypeDescription() );
    additional_entries_label[U] = VarLabel::create ( "A_uu" + HypreSStruct::Solver<DIM>::AdditionalEntriesSuffix, _AdditionalEntries::getTypeDescription() );
    rhs_label[PSI] = VarLabel::create ( "b_psi", Variable<VAR, double>::getTypeDescription() );
    rhs_label[U] = VarLabel::create ( "b_u", Variable<VAR, double>::getTypeDescription() );
#endif
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
PureMetal<VAR, DIM, STN, AMR>::~PureMetal()
{
    delete post_process;

    VarLabel::destroy ( psi_label );
    VarLabel::destroy ( u_label );
    VarLabel::destroy ( grad_psi_norm2_label );
    VarLabel::destroy ( a_label );
    VarLabel::destroy ( a2_label );
    for ( size_t d = 0; d < DIM; ++d )
        VarLabel::destroy ( grad_psi_label[d] );
    for ( size_t d = 0; d < BSZ; ++d )
        VarLabel::destroy ( b_label[d] );

#ifdef HAVE_HYPRE
    VarLabel::destroy ( matrix_label[PSI][PSI] );
    VarLabel::destroy ( matrix_label[PSI][U] );
    VarLabel::destroy ( matrix_label[U][PSI] );
    VarLabel::destroy ( matrix_label[U][U] );
    VarLabel::destroy ( additional_entries_label[PSI] );
    VarLabel::destroy ( additional_entries_label[U] );
    VarLabel::destroy ( rhs_label[PSI] );
    VarLabel::destroy ( rhs_label[U] );
    for ( int i = 0; i < NVARS; ++i ) delete[] matrix_label[i];
    delete[] matrix_label;
    delete[] additional_entries_label;
    delete[] rhs_label;
#endif
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::problemSetup (
    const ProblemSpecP & params,
    const ProblemSpecP & /*restart_prob_spec*/,
    GridP & /*grid*/
)
{
    // register default material
    this->m_materialManager->registerSimpleMaterial ( scinew SimpleMaterial() );

    // read model parameters
    ProblemSpecP pure_metal = params->findBlock ( "PhaseField" );
    pure_metal->require ( "delt", delt );
    pure_metal->require ( "alpha", alpha );
    pure_metal->require ( "R0", r0 );
    pure_metal->require ( "Delta", delta );
    pure_metal->require ( "epsilon", epsilon );
    pure_metal->getWithDefault ( "gamma_psi", gamma_psi, 1. );
    pure_metal->getWithDefault ( "gamma_u", gamma_u, 1. );

    if ( DIM == D3 && epsilon < 0 )
    {
        epsilon = -epsilon;
        time_advance_anisotropy_terms = &PureMetal::time_advance_anisotropy_terms_diag;
    }

    post_process = scinew ArmPostProcessModule<VAR, DIM, STN, AMR> ( this, this->m_regridder, params, psi_label );
    post_process->problemSetup();

    // coupling parameter
    lambda = alpha / 0.6267;

    std::string scheme;
    pure_metal->getWithDefault ( "scheme", scheme, "forward_euler" );
#ifdef HAVE_HYPRE
    time_scheme = str_to_ts ( scheme );

    if ( VAR == NC )
    {
        if ( time_scheme & TS::Implicit )
            SCI_THROW ( InternalError ( "\n ERROR: implicit solver not implemented for node centered variables", __FILE__, __LINE__ ) );
    }

    if ( time_scheme & TS::Implicit )
    {
        ProblemSpecP solv = params->findBlock ( "Solver" );
        this->m_solver = dynamic_cast<SolverInterface *> ( this->getPort ( "solver" ) );
        if ( !this->m_solver )
        {
            SCI_THROW ( InternalError ( "PureMetal:couldn't get solver port", __FILE__, __LINE__ ) );
        }
        this->m_solver->readParameters ( solv, "" );
        this->m_solver->getParameters()->setSymmetric ( false );
        this->m_solver->getParameters()->setSolveOnExtraCells ( false );
    }
#else
    if ( scheme != "forward_euler" )
        SCI_THROW ( InternalError ( "\n ERROR: Implicit time scheme requires HYPRE\n", __FILE__, __LINE__ ) );
#endif
    this->setBoundaryVariables ( psi_label, u_label, a2_label, b_label );

    if ( AMR )
    {
        this->setLockstepAMR ( true );

        // read amr parameters
        pure_metal->require ( "refine_threshold", refine_threshold );

        std::map<std::string, FC> c2f;

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

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleInitialize (
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleInitialize_solution<AMR> ( level, sched );
    scheduleInitialize_grad_psi<AMR> ( level, sched );
    if ( !AMR ) post_process->scheduleInitialize ( sched, level );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleInitialize_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_initialize_solution", this, &PureMetal::task_initialize_solution );
    task->computes ( psi_label );
    task->computes ( u_label );
    task->computes ( a_label );
    task->computes ( a2_label );
    for ( size_t d = 0; d < BSZ; ++d )
        task->computes ( b_label[d] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleInitialize_grad_psi (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_initialize_grad_psi", this, &PureMetal::task_initialize_grad_psi );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, psi_label, FGT, FGN );
    task->computes ( grad_psi_norm2_label );
    for ( size_t d = 0; d < DIM; ++d )
        task->computes ( grad_psi_label[d] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleInitialize_grad_psi (
    const LevelP & /*level*/,
    SchedulerP & /*sched*/
) {}

/**
 * @remark we need to schedule all levels before task_error_estimate_grad_psi to
 * avoid the error "Failure finding [u , coarseLevel, MI: none, NewDW
 * (mapped to dw index 1), ####] for PureMetal::task_error_estimate_grad_psi",
 * on patch #, Level #, on material #, resource (rank): #" while compiling the
 * TaskGraph
 */
template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleInitialize_solution (
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

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleRestartInitialize (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( !level->hasCoarserLevel() ) return;
    Task * task = scinew Task ( "PureMetal::task_communicate_subproblems", this, &PureMetal::task_empty );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->modifies ( this->getSubProblemsLabel() );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void PureMetal<VAR, DIM, STN, AMR>::scheduleComputeStableTimeStep ( LevelP const & level, SchedulerP & sched )
{
    Task * task = scinew Task ( "PureMetal::task_compute_stable_timestep", this, &PureMetal::task_compute_stable_timestep );
    task->computes ( this->getDelTLabel(), level.get_rep() );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );

    is_time_advance_anisotropy_terms_scheduled = false;
    is_time_advance_solution_scheduled = false;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance (
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleTimeAdvance_grad_psi<AMR> ( level, sched );
    scheduleTimeAdvance_anisotropy_terms ( level, sched );
    scheduleTimeAdvance_solution ( level, sched );
    if ( !AMR ) post_process->scheduleDoAnalysis ( sched, level );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_grad_psi (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_time_advance_grad_psi", this, &PureMetal::task_time_advance_grad_psi );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->computes ( grad_psi_norm2_label );
    for ( size_t d = 0; d < DIM; ++d )
        task->computes ( grad_psi_label[d] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_grad_psi (
    const LevelP &,
    SchedulerP &
)
{
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_anisotropy_terms (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( is_time_advance_anisotropy_terms_scheduled ) return;

    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_anisotropy_terms on all levels " );
    GridP grid = level->getGrid();

    for ( int l = 0; l < grid->numLevels(); ++l )
    {
        DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_dbg_derivatives on level " << l << " " );

        Task * task = scinew Task ( "PureMetal::task_time_advance_anisotropy_terms", this, &PureMetal::task_time_advance_anisotropy_terms );
        task->requires ( Task::OldDW, grad_psi_norm2_label, Ghost::None, 0 );
        for ( size_t d = 0; d < DIM; ++d )
            task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
        task->computes ( a_label );
        task->computes ( a2_label );
        for ( size_t d = 0; d < BSZ; ++d )
            task->computes ( b_label[d] );
        sched->addTask ( task, grid->getLevel ( l )->eachPatch(), this->m_materialManager->allMaterials() );
    }

    is_time_advance_anisotropy_terms_scheduled = true;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution (
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution on level " << level->getIndex() << " " );

#ifdef HAVE_HYPRE
    if ( time_scheme == TS::ForwardEuler )
#endif
        scheduleTimeAdvance_solution_forward_euler<AMR> ( level, sched );
#ifdef HAVE_HYPRE
    else
        scheduleTimeAdvance_solution_implicit<AMR> ( level, sched );
#endif
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_forward_euler (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_time_advance_solution", this, &PureMetal::task_time_advance_solution_forward_euler );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    for ( size_t d = 0; d < BSZ; ++d )
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
    task->computes ( psi_label );
    task->computes ( u_label );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_forward_euler (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( is_time_advance_solution_scheduled ) return;

    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_forward_euler on all levels " );
    GridP grid = level->getGrid();

    scheduleTimeAdvance_solution_forward_euler < !MG > ( grid->getLevel ( 0 ), sched );
    for ( int l = 1; l < grid->numLevels(); ++l )
    {
        Task * task = scinew Task ( "PureMetal::task_time_advance_solution_forward_euler", this, &PureMetal::task_time_advance_solution_forward_euler );
        task->requires ( Task::OldDW, this->getSubProblemsLabel(), Ghost::None, 0 );
        task->requires ( Task::OldDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::OldDW, psi_label, FGT, FGN );
        task->requires ( Task::OldDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        for ( size_t d = 0; d < DIM; ++d )
            task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
        task->requires ( Task::OldDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::OldDW, u_label, FGT, FGN );
        task->requires ( Task::OldDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
        task->requires ( Task::NewDW, a2_label, FGT, FGN );
        task->requires ( Task::NewDW, a2_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        for ( size_t d = 0; d < BSZ; ++d )
        {
            task->requires ( Task::NewDW, b_label[d], FGT, FGN );
            task->requires ( Task::NewDW, b_label[d], nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        }
        task->computes ( psi_label );
        task->computes ( u_label );
        sched->addTask ( task, grid->getLevel ( l )->eachPatch(), this->m_materialManager->allMaterials() );
    }
}

#ifdef HAVE_HYPRE
template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_implicit
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_implicit on level " << level->getIndex() << " " );

    void ( PureMetal::*scheduleTimeAdvance_assemble ) ( const LevelP &, SchedulerP & );

    switch ( time_scheme )
    {
    case TS::BackwardEuler:
        scheduleTimeAdvance_assemble = &PureMetal::scheduleTimeAdvance_solution_backward_euler_assemble_hypre<AMR>;
        break;
    case TS::CrankNicolson:
        scheduleTimeAdvance_assemble = &PureMetal::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre<AMR>;
        break;
    default:
        SCI_THROW ( InternalError ( "\n ERROR: Unknown time scheme\n", __FILE__, __LINE__ ) );
    }

    ( this->*scheduleTimeAdvance_assemble ) ( level, sched );

    const VarLabel * labels[NVARS] = { psi_label, u_label };
    this->m_solver->scheduleSolve ( NVARS, level, sched,
                                    this->m_materialManager->allMaterials(),
                                    matrix_label, Task::NewDW, // A
                                    labels, false,             // x
                                    rhs_label, Task::NewDW,    // b
                                    labels, Task::OldDW        // guess
                                  );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_implicit
(
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( is_time_advance_solution_scheduled ) return;

    if ( this->m_solver->getName() == "hypre_sstruct" )
    {
        void ( PureMetal::*scheduleTimeAdvance_assemble_sg ) ( const LevelP &, SchedulerP & );
        void ( PureMetal::*scheduleTimeAdvance_assemble_mg ) ( const LevelP &, SchedulerP & );

        switch ( time_scheme )
        {
        case TS::BackwardEuler:
            scheduleTimeAdvance_assemble_sg = &PureMetal::scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct < !MG >;
            scheduleTimeAdvance_assemble_mg = &PureMetal::scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct<MG>;
            break;
        case TS::CrankNicolson:
            scheduleTimeAdvance_assemble_sg = &PureMetal::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct < !MG >;
            scheduleTimeAdvance_assemble_mg = &PureMetal::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct<MG>;
            break;
        default:
            SCI_THROW ( InternalError ( "\n ERROR: Unknown time scheme\n", __FILE__, __LINE__ ) );
        }

        DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_implicit on all levels " );
        GridP grid = level->getGrid();

        ( this->*scheduleTimeAdvance_assemble_sg ) ( grid->getLevel ( 0 ), sched );

        for ( int l = 1; l < grid->numLevels(); ++l )
            ( this->*scheduleTimeAdvance_assemble_mg ) ( grid->getLevel ( l ), sched );

        for ( int l = 0; l < grid->numLevels(); ++l )
            scheduleTimeAdvannce_comunicate_before_solve_hypresstruct ( grid->getLevel ( l ), sched );

        const VarLabel * labels[NVARS] = { psi_label, u_label };
        this->m_solver->scheduleSolve ( NVARS, level, sched,
                                        this->m_materialManager->allMaterials(),
                                        matrix_label, Task::NewDW, // A
                                        labels, false,             // x
                                        rhs_label, Task::NewDW,    // b
                                        labels, Task::OldDW        // guess
                                      );
    }
    else
    {
        void ( PureMetal::*scheduleTimeAdvance_assemble_sg ) ( const LevelP &, SchedulerP & );
        void ( PureMetal::*scheduleTimeAdvance_assemble_mg ) ( const LevelP &, SchedulerP & );

        switch ( time_scheme )
        {
        case TS::BackwardEuler:
            scheduleTimeAdvance_assemble_sg = &PureMetal::scheduleTimeAdvance_solution_backward_euler_assemble_hypre < !MG >;
            scheduleTimeAdvance_assemble_mg = &PureMetal::scheduleTimeAdvance_solution_backward_euler_assemble_hypre < MG >;
            break;
        case TS::CrankNicolson:
            scheduleTimeAdvance_assemble_sg = &PureMetal::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre < !MG >;
            scheduleTimeAdvance_assemble_mg = &PureMetal::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre < MG >;
            break;
        default:
            SCI_THROW ( InternalError ( "\n ERROR: Unknown time scheme\n", __FILE__, __LINE__ ) );
        }

        GridP grid = level->getGrid();

        ( this->*scheduleTimeAdvance_assemble_sg ) ( grid->getLevel ( 0 ), sched );

        const VarLabel * labels[NVARS] = { psi_label, u_label };
        this->m_solver->scheduleSolve ( NVARS, grid->getLevel ( 0 ), sched,
                                        this->m_materialManager->allMaterials(),
                                        matrix_label, Task::NewDW, // A
                                        labels, false,             // x
                                        rhs_label, Task::NewDW,    // b
                                        labels, Task::OldDW        // guess
                                      );

        for ( int l = 1; l < grid->numLevels(); ++l )
        {
            ( this->*scheduleTimeAdvance_assemble_mg ) ( grid->getLevel ( l ), sched );

            const VarLabel * labels[NVARS] = { psi_label, u_label };
            this->m_solver->scheduleSolve ( NVARS, grid->getLevel ( l ), sched,
                                            this->m_materialManager->allMaterials(),
                                            matrix_label, Task::NewDW, // A
                                            labels, false,             // x
                                            rhs_label, Task::NewDW,    // b
                                            labels, Task::OldDW        // guess
                                          );
        }
    }

    is_time_advance_solution_scheduled = true;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_backward_euler_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_backward_euler_assemble_hypre on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_backward_euler_assemble_hypre", this, &PureMetal::task_time_advance_solution_backward_euler_assemble_hypre );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    for ( size_t d = 0; d < BSZ; ++d )
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_backward_euler_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_backward_euler_assemble_hypre on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_backward_euler_assemble_hypre", this, &PureMetal::task_time_advance_solution_backward_euler_assemble_hypre );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::NewDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    task->requires ( Task::NewDW, a2_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < BSZ; ++d )
    {
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
        task->requires ( Task::NewDW, b_label[d], nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    }
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct", this, &PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    for ( size_t d = 0; d < BSZ; ++d )
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( additional_entries_label[PSI] );
    task->computes ( additional_entries_label[U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_backward_euler_assemble_hypresstruct on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct", this, &PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    task->requires ( Task::NewDW, a2_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < BSZ; ++d )
    {
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
        task->requires ( Task::NewDW, b_label[d], nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    }
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( additional_entries_label[PSI] );
    task->computes ( additional_entries_label[U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre", this, &PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    for ( size_t d = 0; d < BSZ; ++d )
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_crank_nicolson_assemble_hypre on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre", this, &PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::NewDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    task->requires ( Task::NewDW, a2_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < BSZ; ++d )
    {
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
        task->requires ( Task::NewDW, b_label[d], nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    }
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct", this, &PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    for ( size_t d = 0; d < BSZ; ++d )
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( additional_entries_label[PSI] );
    task->computes ( additional_entries_label[U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleTimeAdvance_solution_crank_nicolson_assemble_hypresstruct on level " << level->getIndex() << " " );

    Task * task = scinew Task ( "PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct", this, &PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), FGT, FGN );
    task->requires ( Task::OldDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, psi_label, FGT, FGN );
    task->requires ( Task::OldDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::OldDW, u_label, FGT, FGN );
    task->requires ( Task::OldDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < DIM; ++d )
        task->requires ( Task::OldDW, grad_psi_label[d], Ghost::None, 0 );
    task->requires ( Task::NewDW, a_label, Ghost::None, 0 );
    task->requires ( Task::NewDW, a2_label, FGT, FGN );
    task->requires ( Task::NewDW, a2_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    for ( size_t d = 0; d < BSZ; ++d )
    {
        task->requires ( Task::NewDW, b_label[d], FGT, FGN );
        task->requires ( Task::NewDW, b_label[d], nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    }
    task->requires ( Task::OldDW, matrix_label[PSI][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[PSI][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, matrix_label[U][U], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[PSI], Ghost::None, 0 );
    task->requires ( Task::OldDW, additional_entries_label[U], Ghost::None, 0 );
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    task->computes ( additional_entries_label[PSI] );
    task->computes ( additional_entries_label[U] );
    task->computes ( rhs_label[PSI] );
    task->computes ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleTimeAdvannce_comunicate_before_solve_hypresstruct
(
    const LevelP & level,
    SchedulerP & sched
)
{
    // add empty task to force receiving of all pending mpi messages before hypre solve
    Task * task = scinew Task ( "task_communicate_before_solve", this, &PureMetal::task_empty );

    task->requires ( Task::NewDW, this->getSubProblemsLabel(), FGT, FGN );
    if ( level->hasCoarserLevel() )
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    if ( level->hasFinerLevel() )
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );

    task->modifies ( matrix_label[PSI][PSI] );
    task->modifies ( matrix_label[PSI][U] );
    task->modifies ( matrix_label[U][PSI] );
    task->modifies ( matrix_label[U][U] );
    task->modifies ( additional_entries_label[PSI] );
    task->modifies ( additional_entries_label[U] );
    task->modifies ( rhs_label[PSI] );
    task->modifies ( rhs_label[U] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}
#endif // HAVE_HYPRE

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleRefine
(
    const PatchSet * new_patches,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleRefine on: " << *new_patches );

    const Level * level = getLevel ( new_patches );

    // no need to refine on coarser level
    if ( level->hasCoarserLevel() )
    {
        scheduleRefine_solution ( new_patches, sched );
        scheduleRefine_communicate_psi ( level->getCoarserLevel(), sched );
        scheduleRefine_grad_psi ( new_patches, sched );
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleRefine_solution (
    const PatchSet * patches,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleRefine_solution on: " << *patches );

    Task * task = scinew Task ( "PureMetal::task_refine_solution", this, &PureMetal::task_refine_solution );
    task->requires ( Task::NewDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::NewDW, u_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->computes ( psi_label );
    task->computes ( u_label );
    sched->addTask ( task, patches, this->m_materialManager->allMaterials() );

#ifdef HAVE_HYPRE
    /*************************** WORKAROUND ***************************/
    /* on new patches of finer level need to create matrix variables  */
    task->computes ( matrix_label[PSI][PSI] );
    task->computes ( matrix_label[PSI][U] );
    task->computes ( matrix_label[U][PSI] );
    task->computes ( matrix_label[U][U] );
    /************************* END WORKAROUND *************************/
#endif
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleRefine_grad_psi (
    const PatchSet * patches,
    SchedulerP & sched
)
{
    DOUTR ( dbg_pure_metal_scheduling, "scheduleRefine_grad_psi on: " << *patches );

    Task * task = scinew Task ( "PureMetal::task_refine_grad_psi", this, &PureMetal::task_refine_grad_psi );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->requires ( Task::NewDW, psi_label, FGT, FGN );
    task->requires ( Task::NewDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
    task->computes ( grad_psi_norm2_label );
    for ( size_t d = 0; d < DIM; ++d )
        task->computes ( grad_psi_label[d] );
    sched->addTask ( task, patches, this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleRefineInterface (
    const LevelP & /*level_fine*/,
    SchedulerP & /*sched*/,
    bool /*need_old_coarse*/,
    bool /*need_new_coarse*/
) {}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleCoarsen
(
    const LevelP & level_coarse,
    SchedulerP & sched
)
{
    scheduleCoarsen_solution ( level_coarse, sched );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleCoarsen_solution (
    const LevelP & level_coarse,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_coarsen_solution", this, &PureMetal::task_coarsen_solution );
    task->requires ( Task::NewDW, psi_label, nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );
    task->requires ( Task::NewDW, u_label, nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::FineLevel, nullptr, Task::NormalDomain, Ghost::None, 0 );
    task->modifies ( psi_label );
    task->modifies ( u_label );
    sched->addTask ( task, level_coarse->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::scheduleErrorEstimate
(
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleErrorEstimate_grad_psi<AMR> ( level, sched );
    if ( !level->hasFinerLevel() ) post_process->scheduleDoAnalysis ( sched, level );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < !MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleErrorEstimate_grad_psi (
    const LevelP & level,
    SchedulerP & sched
)
{
    Task * task = scinew Task ( "PureMetal::task_error_estimate_grad_psi", this, &PureMetal::task_error_estimate_grad_psi );
    task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
    task->requires ( Task::NewDW, psi_label, FGT, FGN );
    task->modifies ( this->m_regridder->getRefineFlagLabel(), this->m_regridder->refineFlagMaterials() );
    task->computes ( grad_psi_norm2_label );
    for ( size_t d = 0; d < DIM; ++d )
        task->computes ( grad_psi_label[d] );
    sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
PureMetal<VAR, DIM, STN, AMR>::scheduleErrorEstimate_grad_psi (
    const LevelP & level,
    SchedulerP & sched
)
{
    if ( !level->hasCoarserLevel() ) scheduleErrorEstimate_grad_psi < !MG > ( level, sched );
    else
    {
        Task * task = scinew Task ( "PureMetal::task_error_estimate_grad_psi", this, &PureMetal::task_error_estimate_grad_psi );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), Ghost::None, 0 );
        task->requires ( Task::NewDW, this->getSubProblemsLabel(), nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->requires ( Task::NewDW, psi_label, FGT, FGN );
        task->requires ( Task::NewDW, psi_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
        task->modifies ( this->m_regridder->getRefineFlagLabel(), this->m_regridder->refineFlagMaterials() );
        task->computes ( grad_psi_norm2_label );
        for ( size_t d = 0; d < DIM; ++d )
            task->computes ( grad_psi_label[d] );
        sched->addTask ( task, level->eachPatch(), this->m_materialManager->allMaterials() );
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void PureMetal<VAR, DIM, STN, AMR>::scheduleInitialErrorEstimate
(
    const LevelP & level,
    SchedulerP & sched
)
{
    scheduleErrorEstimate ( level, sched );
}

// TASKS

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_initialize_solution (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * /*matls*/,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_initialize_solution ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch );;

        // Allocate solution variables into the new DataWarehouse
        DWView < ScalarField<double>, VAR, DIM > psi ( dw_new, psi_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > u ( dw_new, u_label, material, patch );

        // Get patch range
        BlockRange range ( this->get_range ( patch ) );
        DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over range " << range );;

        // Initialize solution variables in range
        parallel_for ( range, [patch, &psi, &u, this] ( int i, int j, int k )->void { initialize_solution ( {i, j, k}, patch, psi, u ); } );

        // Allocate anisotropy terms variables into the new DataWarehouse
        DWView < ScalarField<double>, VAR, DIM > a ( dw_new, a_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > a2 ( dw_new, a2_label, material, patch );
        DWView < VectorField<double, BSZ>, VAR, DIM > b ( dw_new, b_label, material, patch );

        // Initialize anisotropy terms variables
        a.initialize ( 0. );
        a2.initialize ( 0. );
        b.initialize ( 0. );
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_initialize_grad_psi (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * /*matls*/,
    DataWarehouse * /*dw_old*/,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_initialize_grad_psi ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch );;

        // Allocate grad_psi and grad_psi_norm2 variables into the new DataWarehouse
        DWView < ScalarField<double>, VAR, DIM > grad_psi_norm2 ( dw_new, grad_psi_norm2_label, material, patch );
        DWView < VectorField<double, DIM>, VAR, DIM > grad_psi ( dw_new, grad_psi_label, material, patch );

        // Retrieve subproblems from the DataWarehouse
        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );

        // Iterate over each subproblem
        for ( const auto & p : subproblems )
        {
            // Get a view of psi that implements finite differences approximations
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over " << p );;
            FDView < ScalarField<const double>, STN > & psi = p.template get_fd_view<PSI> ( dw_new );

            // Compute psi derivatives on subproblem range
            parallel_for ( p.get_range(), [&psi, &grad_psi, &grad_psi_norm2, this] ( int i, int j, int k )->void { time_advance_grad_psi ( {i, j, k}, psi, grad_psi, grad_psi_norm2 ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_compute_stable_timestep (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_compute_stable_timestep ====" );;

    dw_new->put ( delt_vartype ( delt ), this->getDelTLabel(), getLevel ( patches ) );

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_grad_psi (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_time_advance_grad_psi ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch );;

        DWView < ScalarField<double>, VAR, DIM > grad_psi_norm2 ( dw_new, grad_psi_norm2_label, material, patch );
        DWView < VectorField<double, DIM>, VAR, DIM > grad_psi ( dw_new, grad_psi_label, material, patch );

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over " << p );;
            FDView < ScalarField<const double>, STN > & psi = p.template get_fd_view<PSI> ( dw_old );
            parallel_for ( p.get_range(), [&psi, &grad_psi, &grad_psi_norm2, this] ( int i, int j, int k )->void { time_advance_grad_psi ( {i, j, k}, psi, grad_psi, grad_psi_norm2 ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_anisotropy_terms (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_time_advance_anisotropy_terms ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch );;

        DWView < ScalarField<const double>, VAR, DIM > grad_psi_norm2 ( dw_old, grad_psi_norm2_label, material, patch );
        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );

        DWView < ScalarField<double>, VAR, DIM > a ( dw_new, a_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > a2 ( dw_new, a2_label, material, patch );
        DWView < VectorField<double, BSZ>, VAR, DIM > b ( dw_new, b_label, material, patch );

        BlockRange range ( this->get_range ( patch ) );
        DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over range " << range );;
        parallel_for ( range, [&grad_psi, &grad_psi_norm2, &a, &a2, &b, this] ( int i, int j, int k )->void { ( this->*time_advance_anisotropy_terms ) ( {i, j, k}, grad_psi, grad_psi_norm2, a, a2, b ); } );
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_forward_euler (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_time_advance_solution_forward_euler ====" );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch );;

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        DWView < ScalarField<double>, VAR, DIM > psi_new ( dw_new, psi_label, material, patch );
        DWView < ScalarField<double>, VAR, DIM > u_new ( dw_new, u_label, material, patch );

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over " << p );;
            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &psi_new, &u_new, this] ( int i, int j, int k )->void { time_advance_solution_forward_euler ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, psi_new, u_new ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

#ifdef HAVE_HYPRE
template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypre
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
        task_time_advance_solution_backward_euler_assemble_hypre_full ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_backward_euler_assemble_hypre_rhs ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypre_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_backward_euler_assemble_hypre_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View< ScalarField<Stencil7> > * A [NVARS][NVARS] =
        {
            {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][U], material, patch )
            }, {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][U], material, patch )
            }
        };

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &A, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_backward_euler_assemble_hypre_full ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, A, rhs ); } );
        }

        delete A[PSI][PSI];
        delete A[PSI][U];
        delete A[U][PSI];
        delete A[U][U];

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypre_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_backward_euler_assemble_hypre_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label[PSI][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[PSI][U], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][U], patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_backward_euler_assemble_hypre_rhs ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, rhs ); } );
        }

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypresstruct
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
        task_time_advance_solution_backward_euler_assemble_hypresstruct_full ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_backward_euler_assemble_hypresstruct_rhs ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypresstruct_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View< ScalarField<Stencil7> > * A_stencil [NVARS][NVARS] =
        {
            {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][U], material, patch )
            }, {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][U], material, patch )
            }
        };

        HypreSStruct::AdditionalEntries * A_additional[NVARS]
        {
            scinew HypreSStruct::AdditionalEntries,
            scinew HypreSStruct::AdditionalEntries
        };

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &A_stencil, &A_additional, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_backward_euler_assemble_hypresstruct_full ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, A_stencil, A_additional, rhs ); } );
        }

        delete A_stencil[PSI][PSI];
        delete A_stencil[PSI][U];
        delete A_stencil[U][PSI];
        delete A_stencil[U][U];

        delete rhs[PSI];
        delete rhs[U];

        _AdditionalEntries additional_entries[NVARS];
        additional_entries[PSI].setData ( A_additional[PSI] );
        additional_entries[U].setData ( A_additional[U] );
        dw_new->put ( additional_entries[PSI], additional_entries_label[PSI], material, patch );
        dw_new->put ( additional_entries[U], additional_entries_label[U], material, patch );
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_backward_euler_assemble_hypresstruct_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_backward_euler_assemble_hypresstruct_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label[PSI][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[PSI][U], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][U], patches, matls );
    dw_new->transferFrom ( dw_old, additional_entries_label[PSI], patches, matls );
    dw_new->transferFrom ( dw_old, additional_entries_label[U], patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_backward_euler_assemble_hypresstruct_rhs ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, rhs ); } );
        }

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypre
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
        task_time_advance_solution_crank_nicolson_assemble_hypre_full ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_crank_nicolson_assemble_hypre_rhs ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypre_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View< ScalarField<Stencil7> > * A [NVARS][NVARS] =
        {
            {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][U], material, patch )
            }, {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][U], material, patch )
            }
        };

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &A, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_crank_nicolson_assemble_hypre_full ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, A, rhs ); } );
        }

        delete A[PSI][PSI];
        delete A[PSI][U];
        delete A[U][PSI];
        delete A[U][U];

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypre_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypre_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label[PSI][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[PSI][U], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][U], patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_crank_nicolson_assemble_hypre_rhs ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, rhs ); } );
        }

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypresstruct
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
        task_time_advance_solution_crank_nicolson_assemble_hypresstruct_full ( myworld, patches, matls, dw_old, dw_new );
    else
        task_time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs ( myworld, patches, matls, dw_old, dw_new );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypresstruct_full
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * /*matls*/,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct_full ====" );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View< ScalarField<Stencil7> > * A_stencil [NVARS][NVARS] =
        {
            {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[PSI][U], material, patch )
            }, {
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][PSI], material, patch ),
                scinew DWView< ScalarField<Stencil7>, VAR, DIM > ( dw_new, matrix_label[U][U], material, patch )
            }
        };

        HypreSStruct::AdditionalEntries * A_additional[NVARS]
        {
            scinew HypreSStruct::AdditionalEntries,
            scinew HypreSStruct::AdditionalEntries
        };

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &A_stencil, &A_additional, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_backward_euler_assemble_hypresstruct_full ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, A_stencil, A_additional, rhs ); } );
        }

        delete A_stencil[PSI][PSI];
        delete A_stencil[PSI][U];
        delete A_stencil[U][PSI];
        delete A_stencil[U][U];

        delete rhs[PSI];
        delete rhs[U];

        _AdditionalEntries additional_entries[NVARS];
        additional_entries[PSI].setData ( A_additional[PSI] );
        additional_entries[U].setData ( A_additional[U] );
        dw_new->put ( additional_entries[PSI], additional_entries_label[PSI], material, patch );
        dw_new->put ( additional_entries[U], additional_entries_label[U], material, patch );
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset * matls,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1, myrank << "==== PureMetal::task_time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs ====" );

    dw_new->transferFrom ( dw_old, matrix_label[PSI][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[PSI][U], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][PSI], patches, matls );
    dw_new->transferFrom ( dw_old, matrix_label[U][U], patches, matls );
    dw_new->transferFrom ( dw_old, additional_entries_label[PSI], patches, matls );
    dw_new->transferFrom ( dw_old, additional_entries_label[U], patches, matls );

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " " );

        DWView < VectorField<const double, DIM>, VAR, DIM > grad_psi ( dw_old, grad_psi_label, material, patch );
        DWView < ScalarField<const double>, VAR, DIM > a ( dw_new, a_label, material, patch );

        View < ScalarField<double> > * rhs[NVARS] =
        {
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[PSI], material, patch ),
            scinew DWView < ScalarField<double>, VAR, DIM > ( dw_new, rhs_label[U], material, patch )
        };

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_old, this->getSubProblemsLabel(), material, patch );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3, myrank << "= Iterating over " << p );

            FDView < ScalarField<const double>, STN > & psi_old = p.template get_fd_view<PSI> ( dw_old );
            FDView < ScalarField<const double>, STN > & u_old = p.template get_fd_view<U> ( dw_old );
            FDView < ScalarField<const double>, STN > & a2 = p.template get_fd_view<A2> ( dw_new );
            FDView < VectorField<const double, BSZ>, STN > b = p.template get_fd_view<B> ( dw_new );
            parallel_for ( p.get_range(), [&psi_old, &u_old, &grad_psi, &a, &a2, &b, &rhs, this] ( int i, int j, int k )->void { time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs ( {i, j, k}, psi_old, u_old, grad_psi, a, a2, b, rhs ); } );
        }

        delete rhs[PSI];
        delete rhs[U];
    }

    DOUT ( this->m_dbg_lvl2, myrank );
}
#endif // HAVE_HYPRE

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_refine_solution
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches_fine,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_refine_solution ====" );;

    for ( int p = 0; p < patches_fine->size(); ++p )
    {
        const Patch * patch_fine = patches_fine->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Fine Patch: " << *patch_fine << " Level: " << patch_fine->getLevel()->getIndex() );;

        DWView < ScalarField<double>, VAR, DIM > psi_fine ( dw_new, psi_label, material, patch_fine );
        DWView < ScalarField<double>, VAR, DIM > u_fine ( dw_new, u_label, material, patch_fine );

        AMRInterpolator < PureMetalProblem<VAR, STN>, PSI, C2F > psi_coarse_interp ( dw_new, psi_label, this->getSubProblemsLabel(), material, patch_fine );
        AMRInterpolator < PureMetalProblem<VAR, STN>, U, C2F > u_coarse_interp ( dw_new, u_label, this->getSubProblemsLabel(), material, patch_fine );

        BlockRange range_fine ( this->get_range ( patch_fine ) );
        DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over fine range" << range_fine );;
        parallel_for ( range_fine, [&psi_coarse_interp, &u_coarse_interp, &psi_fine, &u_fine, this] ( int i, int j, int k )->void { refine_solution ( {i, j, k}, psi_coarse_interp, u_coarse_interp, psi_fine, u_fine ); } );
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_refine_grad_psi
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches_fine,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_refine_grad_psi ====" );;

    for ( int p = 0; p < patches_fine->size(); ++p )
    {
        const Patch * patch_fine = patches_fine->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Fine Patch: " << *patch_fine << " Level: " << patch_fine->getLevel()->getIndex() );;

        DWView < ScalarField<double>, VAR, DIM > grad_psi_norm2 ( dw_new, grad_psi_norm2_label, material, patch_fine );
        DWView < VectorField<double, DIM>, VAR, DIM > grad_psi ( dw_new, grad_psi_label, material, patch_fine );

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch_fine );

        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over " << p );;
            FDView < ScalarField<const double>, STN > & psi = p.template get_fd_view<PSI> ( dw_new );
            parallel_for ( p.get_range(), [&psi, &grad_psi, &grad_psi_norm2, this] ( int i, int j, int k )->void { time_advance_grad_psi ( {i, j, k}, psi, grad_psi, grad_psi_norm2 ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_coarsen_solution (
    const ProcessorGroup * myworld,
    const PatchSubset * patches_coarse,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_coarsen_solution " );;

    for ( int p = 0; p < patches_coarse->size(); ++p )
    {
        const Patch * patch_coarse = patches_coarse->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Coarse Patch: " << *patch_coarse << " Level: " << patch_coarse->getLevel()->getIndex() );;

        DWView < ScalarField<double>, VAR, DIM > psi_coarse ( dw_new, psi_label, material, patch_coarse );
        DWView < ScalarField<double>, VAR, DIM > u_coarse ( dw_new, u_label, material, patch_coarse );

        AMRRestrictor < PureMetalProblem<VAR, STN>, PSI, F2C > psi_fine_restr ( dw_new, psi_label, this->getSubProblemsLabel(), material, patch_coarse, false );
        AMRRestrictor < PureMetalProblem<VAR, STN>, U, F2C > u_fine_restr ( dw_new, u_label, this->getSubProblemsLabel(), material, patch_coarse, false );

        for ( const auto & region : u_fine_restr.get_support() )
        {
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over coarse cells region " << region );;
            BlockRange range_coarse (
                Max ( region.getLow(), this->get_low ( patch_coarse ) ),
                Min ( region.getHigh(), this->get_high ( patch_coarse ) )
            );

            parallel_for ( range_coarse, [&psi_fine_restr, &u_fine_restr, &psi_coarse, &u_coarse, this] ( int i, int j, int k )->void { coarsen_solution ( {i, j, k}, psi_fine_restr, u_fine_restr, psi_coarse, u_coarse ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::task_error_estimate_grad_psi
(
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse *,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();

    DOUT ( this->m_dbg_lvl1,  myrank << "==== PureMetal::task_error_estimate_grad_psi " );;

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( this->m_dbg_lvl2,  myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );;

        DWView < ScalarField<double>, VAR, DIM > grad_psi_norm2 ( dw_new, grad_psi_norm2_label, material, patch );
        DWView < VectorField<double, DIM>, VAR, DIM > grad_psi ( dw_new, grad_psi_label, material, patch );

        DWView < ScalarField<int>, CC, DIM > refine_flag ( dw_new, this->m_regridder->getRefineFlagLabel(), material, patch );
        refine_flag.initialize ( 0 );

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, this->getSubProblemsLabel(), material, patch );
        for ( const auto & p : subproblems )
        {
            DOUT ( this->m_dbg_lvl3,  myrank << "= Iterating over " << p );;
            FDView < ScalarField<const double>, STN > & psi = p.template get_fd_view<PSI> ( dw_new );
            parallel_for ( p.get_range(), [&psi, &grad_psi, &grad_psi_norm2, &refine_flag, this] ( int i, int j, int k )->void { error_estimate_grad_psi ( {i, j, k}, psi, grad_psi, grad_psi_norm2, refine_flag ); } );
        }
    }

    DOUT ( this->m_dbg_lvl2,  myrank );;
}

// IMPLEMENTATIONS

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::initialize_solution (
    const IntVector & id,
    const Patch * patch,
    View < ScalarField<double> > & psi,
    View < ScalarField<double> > & u
)
{
    Vector v ( this->get_position ( patch, id ).asVector() );
    double r2 = 0;

    for ( size_t d = 0; d < DIM; ++d )
        r2 += v[d] * v[d];

    double tmp = r2 - r0 * r0;

    psi[id] = - tanh ( gamma_psi * tmp );
    u[id] = -delta * ( 1. + tanh ( gamma_u * tmp ) ) / 2.;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_grad_psi (
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi,
    View < VectorField<double, DIM> > & grad_psi,
    View < ScalarField<double> > & grad_psi_norm2
)
{
    auto grad = psi.gradient ( id );

    grad_psi_norm2[id] = 0;
    for ( size_t d = 0; d < DIM; ++d )
    {
        grad_psi[d][id] = grad[d];
        grad_psi_norm2[id] += grad[d] * grad[d];
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void PureMetal<VAR, DIM, STN, AMR>::time_advance_anisotropy_terms_dflt (
    const IntVector & id,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & grad_psi_norm2,
    View < ScalarField<double> > & a,
    View < ScalarField<double> > & a2,
    View < VectorField<double, BSZ> > & b
)
{
    double tmp = 1. + epsilon;
    double n2 = grad_psi_norm2[id];
    if ( n2 < tol )
    {
        a[id] = tmp;
        a2[id] = tmp * tmp;
        for ( size_t d = 0; d < BSZ; ++d )
            b[d][id] = 0.;
    }
    else
    {
        double n4 = n2 * n2;

        double tmp4 = 0.;
        double grad[DIM], grad2[DIM];
        for ( size_t d = 0; d < DIM; ++d )
        {
            grad[d] = grad_psi[d][id];
            grad2[d] = grad[d] * grad[d];
            tmp4 += grad2[d] * grad2[d];
        }
        tmp4 *= 4. / n4;

        double tmp = 1. + epsilon * ( tmp4 - 3. );
        a[id] = tmp;
        a2[id] = tmp * tmp;
        if ( DIM > 1 )
        {
            b[XY][id] = 16. * epsilon * tmp * ( grad[X] * grad[Y] )  * ( grad2[X] - grad2[Y] ) / n4;
        }
        if ( DIM > 2 ) // Compile time if
        {
            b[XZ][id] = 16. * epsilon * tmp * ( grad[X] * grad[Z] ) * ( grad2[X] - grad2[Z] ) / n4;
            b[YZ][id] = 16. * epsilon * tmp * ( grad[Y] * grad[Z] ) * ( grad2[Y] - grad2[Z] ) / n4;
        }
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void PureMetal<VAR, DIM, STN, AMR>::time_advance_anisotropy_terms_diag (
    const IntVector & id,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & grad_psi_norm2,
    View < ScalarField<double> > & a,
    View < ScalarField<double> > & a2,
    View < VectorField<double, BSZ> > & b
)
{
    double tmp = 1. + epsilon;
    double n2 = grad_psi_norm2[id];
    if ( n2 < tol )
    {
        a[id] = tmp;
        a2[id] = tmp * tmp;
        for ( size_t d = 0; d < BSZ; ++d )
            b[d][id] = 0.;
    }
    else
    {
        double n4 = n2 * n2 ;

        double grad[DIM], grad2[DIM];
        for ( size_t d = 0; d < DIM; ++d )
        {
            grad[d] = grad_psi[d][id];
            grad2[d] = grad[d] * grad[d];
        }

        double sum4 = grad[X] + grad[Y];
        sum4 *= sum4;
        sum4 *= sum4;
        double dif4 = grad[X] - grad[Y];
        dif4 *= dif4;
        dif4 *= dif4;
        double tmp4 = grad2[Z];
        tmp4 *= tmp4;
        tmp4 *= 4.;
        tmp4 += sum4 + dif4;
        tmp4 /= n4;

        double tmp = 1. + epsilon * ( tmp4 - 3. );
        a[id] = tmp;
        a2[id] = tmp * tmp;
        b[XY][id] = ( 16. * epsilon * grad[X] * grad[Y] * ( grad2[Y] - grad2[X] ) ) / n4;
        b[XZ][id] = ( 8. * epsilon * grad[X] * grad[Z] * ( grad2[X] + 3. * grad2[Y] - 2. * grad2[Z] ) ) / n4;
        b[YZ][id] = ( 8. * epsilon * grad[Y] * grad[Z] * ( 3. * grad2[X] + grad2[Y] - 2. * grad2[Z] ) ) / n4;
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_forward_euler
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<double> > & psi_new,
    View < ScalarField<double> > & u_new
)
{
    double source = 1. - psi_old[id] * psi_old[id];
    source *= ( psi_old[id] - lambda * u_old[id] * source );

    double delta_psi = delt * ( psi_old.laplacian ( id ) * a2[id]
                                + anisotropy_term ( id, grad_psi, a2, b )
                                + source ) / a[id];

    double delta_u = delt * u_old.laplacian ( id ) * alpha + delta_psi / 2.;

    psi_new[id] = psi_old[id] + delta_psi;
    u_new[id] = u_old[id] + delta_u;
}

#ifdef HAVE_HYPRE
template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_backward_euler_assemble_hypre_full
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = delt * a[id];
    const double c1_psi = delt / a[id];
    const double c_u = delt * alpha;

    auto psi_sys = psi_old.laplacian_sys_hypre ( id );
    Stencil7 psi_stn = std::get<0> ( psi_sys );
    double psi_rhs = std::get<1> ( psi_sys );

    auto u_sys = u_old.laplacian_sys_hypre ( id );
    Stencil7 u_stn = std::get<0> ( u_sys );
    double u_rhs = std::get<1> ( u_sys );

    for ( int i = 0; i < 7; ++i )
    {
        ( *A[PSI][PSI] ) [id][i] = -c0_psi * psi_stn[i];
        ( *A[U][U] ) [id][i] = -c_u * u_stn[i];
    }

    ( *A[PSI][PSI] ) [id].p += 1. - c1_psi * source_psi;
    ( *A[PSI][U] ) [id].p += - c1_psi * source_u;
    ( *A[U][PSI] ) [id].p += -1.;
    ( *A[U][U] ) [id].p += 1.;
    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs;
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_backward_euler_assemble_hypre_rhs
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    const double c0_psi = delt * a[id];
    const double c1_psi = delt / a[id];
    const double c_u = delt * alpha;

    const double & psi_rhs = psi_old.laplacian_rhs_hypre ( id );
    const double & u_rhs = u_old.laplacian_rhs_hypre ( id );

    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs;
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_backward_euler_assemble_hypresstruct_full
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
    HypreSStruct::AdditionalEntries * ( & A_additional ) [NVARS],
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = delt * a[id];
    const double c1_psi = delt / a[id];
    const double c_u = delt * alpha;

    auto psi_sys = psi_old.laplacian_sys_hypresstruct ( id );
    Stencil7 psi_stn = std::get<0> ( psi_sys );
    HypreSStruct::AdditionalEntries psi_add = std::get<1> ( psi_sys );
    double psi_rhs = std::get<2> ( psi_sys );

    auto u_sys = u_old.laplacian_sys_hypresstruct ( id );
    Stencil7 u_stn = std::get<0> ( u_sys );
    HypreSStruct::AdditionalEntries u_add = std::get<1> ( u_sys );
    double u_rhs = std::get<2> ( u_sys );

    for ( int i = 0; i < 7; ++i )
    {
        ( *A[PSI][PSI] ) [id][i] = -c0_psi * psi_stn[i];
        ( *A[U][U] ) [id][i] = -c_u * u_stn[i];
    }

    for ( auto & entry : psi_add )
        *A_additional[PSI] += -c0_psi * entry;

    for ( auto & entry : u_add )
        *A_additional[U] += -c_u * entry;

    ( *A[PSI][PSI] ) [id].p += 1. - c1_psi * source_psi;
    ( *A[PSI][U] ) [id].p += - c1_psi * source_u;
    ( *A[U][PSI] ) [id].p += -1.;
    ( *A[U][U] ) [id].p += 1.;

    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs;
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_backward_euler_assemble_hypresstruct_rhs
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    const double c0_psi = delt * a[id];
    const double c1_psi = delt / a[id];
    const double c_u = delt * alpha;

    double psi_rhs = psi_old.laplacian_rhs_hypresstruct ( id );
    double u_rhs = u_old.laplacian_rhs_hypresstruct ( id );

    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs;
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs;
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_crank_nicolson_assemble_hypre_full
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = 0.5 * delt * a[id];
    const double c1_psi = 0.5 * delt / a[id];
    const double c_u = 0.5 * delt * alpha;

    auto psi_sys = psi_old.laplacian_sys_hypre ( id );
    Stencil7 psi_stn = std::get<0> ( psi_sys );
    double psi_rhs = std::get<1> ( psi_sys );

    auto u_sys = u_old.laplacian_sys_hypre ( id );
    Stencil7 u_stn = std::get<0> ( u_sys );
    double u_rhs = std::get<1> ( u_sys );

    for ( int i = 0; i < 7; ++i )
    {
        ( *A[PSI][PSI] ) [id][i] = -c0_psi * psi_stn[i];
        ( *A[U][U] ) [id][i] = -c_u * u_stn[i];
    }

    ( *A[PSI][PSI] ) [id].p += 1. - c1_psi * source_psi;
    ( *A[PSI][U] ) [id].p += - c1_psi * source_u;
    ( *A[U][PSI] ) [id].p += -1.;
    ( *A[U][U] ) [id].p += 1.;
    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs
                         + c0_psi * psi_old.laplacian ( id ) + c1_psi * ( source_psi * psi_old[id] + source_u * u_old[id] );
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs + c_u * u_old.laplacian ( id );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_crank_nicolson_assemble_hypre_rhs
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = 0.5 * delt * a[id];
    const double c1_psi = 0.5 * delt / a[id];
    const double c_u = 0.5 * delt * alpha;

    double psi_rhs = psi_old.laplacian_rhs_hypre ( id );
    double u_rhs = u_old.laplacian_rhs_hypre ( id );

    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs
                         + c0_psi * psi_old.laplacian ( id ) + c1_psi * ( source_psi * psi_old[id] + source_u * u_old[id] );
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs + c_u * u_old.laplacian ( id );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_crank_nicolson_assemble_hypresstruct_full
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<Stencil7> > * ( & A ) [NVARS][NVARS],
    HypreSStruct::AdditionalEntries * ( & A_additional ) [NVARS],
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = 0.5 * delt * a[id];
    const double c1_psi = 0.5 * delt / a[id];
    const double c_u = 0.5 * delt * alpha;

    auto psi_sys = psi_old.laplacian_sys_hypresstruct ( id );
    Stencil7 psi_stn = std::get<0> ( psi_sys );
    HypreSStruct::AdditionalEntries psi_add = std::get<1> ( psi_sys );
    double psi_rhs = std::get<2> ( psi_sys );

    auto u_sys = u_old.laplacian_sys_hypresstruct ( id );
    Stencil7 u_stn = std::get<0> ( u_sys );
    HypreSStruct::AdditionalEntries u_add = std::get<1> ( u_sys );
    double u_rhs = std::get<2> ( u_sys );

    for ( int i = 0; i < 7; ++i )
    {
        ( *A[PSI][PSI] ) [id][i] = -c0_psi * psi_stn[i];
        ( *A[U][U] ) [id][i] = -c_u * u_stn[i];
    }

    for ( auto & entry : psi_add )
        *A_additional[PSI] += -c0_psi * entry;

    for ( auto & entry : u_add )
        *A_additional[U] += -c_u * entry;

    ( *A[PSI][PSI] ) [id].p += 1. - c1_psi * source_psi;
    ( *A[PSI][U] ) [id].p += - c1_psi * source_u;
    ( *A[U][PSI] ) [id].p += -1.;
    ( *A[U][U] ) [id].p += 1.;
    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs
                         + c0_psi * psi_old.laplacian ( id ) + c1_psi * ( source_psi * psi_old[id] + source_u * u_old[id] );
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs + c_u * u_old.laplacian ( id );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::time_advance_solution_crank_nicolson_assemble_hypresstruct_rhs
(
    const IntVector & id,
    const FDView < ScalarField<const double>, STN > & psi_old,
    const FDView < ScalarField<const double>, STN > & u_old,
    const View < VectorField<const double, DIM> > & grad_psi,
    const View < ScalarField<const double> > & a,
    const FDView < ScalarField<const double>, STN > & a2,
    const FDView < VectorField<const double, BSZ>, STN > & b,
    View < ScalarField<double> > * ( & rhs ) [NVARS]
)
{
    double source_psi = 1. - psi_old[id] * psi_old[id];
    double source_u = - lambda * source_psi * source_psi;

    const double c0_psi = 0.5 * delt * a[id];
    const double c1_psi = 0.5 * delt / a[id];
    const double c_u = 0.5 * delt * alpha;

    double psi_rhs = psi_old.laplacian_rhs_hypresstruct ( id );
    double u_rhs = u_old.laplacian_rhs_hypresstruct ( id );

    ( *rhs[PSI] ) [id] = psi_old[id] + c1_psi * anisotropy_term ( id, grad_psi, a2, b ) + c0_psi * psi_rhs
                         + c0_psi * psi_old.laplacian ( id ) + c1_psi * ( source_psi * psi_old[id] + source_u * u_old[id] );
    ( *rhs[U] ) [id] = u_old[id] - psi_old[id] + c_u * u_rhs + c_u * u_old.laplacian ( id );
}
#endif

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::refine_solution
(
    const IntVector id_fine,
    const View < ScalarField<const double> > & psi_coarse_interp,
    const View < ScalarField<const double> > & u_coarse_interp,
    View < ScalarField<double> > & psi_fine,
    View < ScalarField<double> > & u_fine
)
{
    psi_fine[id_fine] = psi_coarse_interp[id_fine];
    u_fine[id_fine] = u_coarse_interp[id_fine];
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::coarsen_solution
(
    const IntVector id_coarse,
    const View < ScalarField<const double> > & psi_fine_restr,
    const View < ScalarField<const double> > & u_fine_restr,
    View < ScalarField<double> > & psi_coarse,
    View < ScalarField<double> > & u_coarse
)
{
    psi_coarse[id_coarse] = psi_fine_restr[id_coarse];
    u_coarse[id_coarse] = u_fine_restr[id_coarse];
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
PureMetal<VAR, DIM, STN, AMR>::error_estimate_grad_psi (
    const IntVector & id,
    FDView < ScalarField<const double>, STN > & psi,
    View < VectorField<double, DIM> > & grad_psi,
    View < ScalarField<double> > & grad_psi_norm2,
    View< ScalarField<int> > & refine_flag
)
{
    bool refine = false;
    auto grad = psi.gradient ( id );

    grad_psi_norm2[id] = 0;
    for ( size_t d = 0; d < DIM; ++d )
    {
        grad_psi[d][id] = grad[d];
        grad_psi_norm2[id] += grad[d] * grad[d];
    }

    refine = grad_psi_norm2[id] > refine_threshold * refine_threshold;
    if ( VAR == CC ) // static if
    {
        if ( refine )
            refine_flag[id] = 1;
    }
    else
    {
        if ( refine )
        {
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
}

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_Applications_PureMetal_h

