/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
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
 * @file CCA/Components/PhaseField/Applications/Application.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Applications_Application_h
#define Packages_Uintah_CCA_Components_PhaseField_Applications_Application_h

#include <CCA/Components/Application/ApplicationCommon.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWInterface.h>
#include <CCA/Components/PhaseField/BoundaryConditions/BCInterface.h>
#include <CCA/Components/PhaseField/AMR/AMRInterface.h>
#include <CCA/Components/PhaseField/DataTypes/SubProblems.h>

namespace Uintah
{
namespace PhaseField
{

// Circular inclusion SubProblems/Application
template < typename Problem > class SubProblems;

/// Debugging stream for component schedulings
static constexpr bool dbg_scheduling = false;

/**
 * @brief Virtual base for PhaseField applications
 *
 * Wrapper of ApplicationCommon and interfaces
 *
 * @remark Application < Problem, false > is used as base for Application < Problem, true >
 *
 * @tparam Problem type of PhaseField problem
 * @tparam AMR whether to use adaptive mesh refinement
*/
template < typename Problem, bool AMR = false > class Application;

/**
 * @brief Virtual base for PhaseField applications (non-AMR implementation)
 *
 * Wrapper of ApplicationCommon and interfaces
 *
 * @implements Application < Problem, AMR >
 *
 * @tparam VAR type of variable representation
 * @tparam STN finite-difference stencil
 */
template < VarType VAR, StnType STN, typename... Field >
class Application < Problem<VAR, STN, Field... >, false >
    : public ApplicationCommon
    , public DWInterface < VAR, get_stn<STN>::dim >
    , public BCInterface < VAR, STN >
{
protected: // STATIC MEMBERS

    /// If SubProblems are required (i.e. if Problem has any boundary Field)
    static constexpr bool use_subprblems = sizeof... ( Field );

    /// Number of ghost elements required by STN (on the same level)
    static constexpr int FGN = get_stn<STN>::ghosts;

    /// Type of ghost elements required by VAR and STN (on the same level)
    static constexpr Ghost::GhostType FGT = FGN ? get_var<VAR>::ghost_type : Ghost::None;

    /// Number of ghost elements required by STN (on coarser level)
    /// @remark this should depend on FCI bc type but if fixed for simplicity
    static constexpr int CGN = 1;

    /// Type of ghost elements required by VAR and STN (on coarser level)
    static constexpr Ghost::GhostType CGT = CGN ? get_var<VAR>::ghost_type : Ghost::None;

    /// Interpolation type for refinement
    static constexpr FCIType C2F = ( VAR == CC ) ? I0 : I1; // TODO make template parameter

    /// Restriction type for coarsening
    static constexpr FCIType F2C = ( VAR == CC ) ? I1 : I0; // TODO make template parameter

protected: // MEMBERS

    /// Verbosity level 1
    const bool m_dbg_lvl1;

    /// Verbosity level 2
    const bool m_dbg_lvl2;

    /// Verbosity level 3
    const bool m_dbg_lvl3;

    /// Verbosity level 3
    const bool m_dbg_lvl4;

    /// List of labels for boundary variables
    std::tuple< typename Field::label_type... > * m_boundary_labels;

    /// Label for subproblems in the DataWarehouse
    const VarLabel * m_subproblems_label;

    /// Store which fine/coarse interface conditions to use on each variable
    std::map<std::string, FC> * m_c2f;

#ifdef HAVE_HYPRE
    /// Time advance scheme
    TS m_time_scheme;

    /// Implicit solver
    SolverInterface * m_solver;
#endif

public: // CONSTRUCTORS/DESTRUCTOR

    /**
     * @brief Constructor
     *
     * Intantiate a PhaseField application
     *
     * @remark m_subproblems_label is initialized only if problem has any
     *         boundary field
     * @remark m_boundary_labels and m_c2f are left uninitialized
     *         for being set later by Application implementation problemSetup
     *
     * @param myWorld data structure to manage mpi processes
     * @param materialManager data structure to manage materials
     * @param verbosity constrols amount of debugging output
     */
    Application (
        const ProcessorGroup * myWorld,
        const MaterialManagerP materialManager,
        int verbosity = 0
    ) : ApplicationCommon ( myWorld, materialManager ),
        m_dbg_lvl1 ( verbosity > 0 ),
        m_dbg_lvl2 ( verbosity > 1 ),
        m_dbg_lvl3 ( verbosity > 2 ),
        m_dbg_lvl4 ( verbosity > 3 ),
        m_boundary_labels ( nullptr ),
        m_subproblems_label ( nullptr ),
        m_c2f ( nullptr )
#ifdef HAVE_HYPRE
        , m_solver ( nullptr )
#endif
    {
        if ( use_subprblems )
            m_subproblems_label = VarLabel::create ( "subproblems", SubProblems< Problem<VAR, STN, Field... > >::getTypeDescription() );
    }

    /**
     * @brief Destructor
     */
    ~Application ()
    {
        VarLabel::destroy ( m_subproblems_label );
    }

public: // GETTERS

    /**
     * @brief Get subproblems variable label
     *
     * @return label for subproblems variable in the DataWarehouse
     */
    virtual const VarLabel *
    getSubProblemsLabel() const
    {
        return m_subproblems_label;
    }

    /**
     * @brief Get boundary variables labels
     *
     * @return list of labels for bounrary variables in the DataWarehouse
     */
    std::tuple< typename Field::label_type... > *
    getBoundaryLabels() const
    {
        return m_boundary_labels;
    }

protected: // SETTERS/GETTERS

    /**
     * @brief Set the list of labels for boundary variables
     *
     * Initialize m_subproblems_label
     *
     * @remark must be called by Application implementation problemSetup if problem has any
     *         boundary field
     *
     * @param labels list of labels for bounrary variables
     */
    void
    setBoundaryVariables (
        const typename Field::label_type & ... labels
    )
    {
        ASSERT ( !m_boundary_labels );
        m_boundary_labels = scinew std::tuple< typename Field::label_type...  > { labels... };
    }

    /**
     * @brief Set the mapping between variable names and C2F condition for amr grids
     *
     * Initialize m_c2f
     *
     * @remark must be called by Application implementation problemSetup if problem has any
     *         boundary field
     *
     * @param c2f mapping between variable names and C2F condition for amr grids
     */
    void
    setC2F (
        const std::map<std::string, FC> & c2f
    )
    {
        ASSERT ( !m_c2f );
        m_c2f = scinew std::map<std::string, FC> ( c2f );
    }

protected: // SCHEDULINGS

    /**
     * @brief Schedule task_initialize_subproblems
     *
     * Defines the dependencies and output of the task which initializes the
     * subproblems allowing sched to control its execution order
     *
     * @param grid grid to be initialized
     * @param perProcPatchSet set of all patches assigned to process
     * @param scheduler scheduler to manage the tasks
     */
    virtual void
    scheduleInitializeSystemVars (
        const GridP & grid,
        const PatchSet * perProcPatchSet,
        SchedulerP & scheduler
    ) override
    {
        ApplicationCommon::scheduleInitializeSystemVars ( grid, perProcPatchSet, scheduler );

        if ( use_subprblems )
        {
            ASSERTMSG ( m_boundary_labels, "Application uses subproblems. Missing call to setBoundaryVariables()" );

            // set behaviour noCheckpoint
            DOUTR ( dbg_scheduling, "scheduleInitializeSystemVars" );

            scheduler->overrideVariableBehavior ( m_subproblems_label->getName(), false, false, false, false, true );

            for ( int idx = 0; idx < grid->numLevels(); ++idx )
            {
                Task * task = scinew Task ( "Application::task_initialize_subproblems", this, &Application::task_initialize_subproblems );
                task->computes ( m_subproblems_label );
                scheduler->addTask ( task, grid->getLevel ( idx )->eachPatch(), this->m_materialManager->allMaterials() );
            }

#ifdef HAVE_HYPRE
            if ( m_solver )
            {
                for ( int idx = 0; idx < grid->numLevels(); ++idx )
                {
                    m_solver->scheduleInitialize ( grid->getLevel ( idx ), scheduler, this->m_materialManager->allMaterials() );
                }
            }
#endif
        }
    }

    /**
     * @brief Schedule task_time_advance_subproblems and task_time_advance_subproblems_foreign
     *
     * Defines the dependencies and output of the task which initializes the
     * subproblems allowing sched to control its execution order
     *
     * @param grid grid to be initialized
     * @param perProcPatchSet set of all patches assigned to process
     * @param scheduler scheduler to manage the tasks
     */
    virtual void
    scheduleAdvanceSystemVars (
        const GridP & grid,
        const PatchSet * perProcPatchSet,
        SchedulerP & scheduler
    ) override
    {
        ApplicationCommon::scheduleAdvanceSystemVars ( grid, perProcPatchSet, scheduler );

        if ( use_subprblems )
        {
            DOUTR ( dbg_scheduling, "scheduleAdvanceSystemVars" );

            for ( int idx = 0; idx < grid->numLevels(); ++idx )
            {
                Task * task = scinew Task ( "Application::task_time_advance_subproblems", this, &Application::task_time_advance_subproblems );
                task->requires ( Task::OldDW, m_subproblems_label, Ghost::None, 0 );
                task->computes ( m_subproblems_label );
                scheduler->addTask ( task, grid->getLevel ( idx )->eachPatch(), this->m_materialManager->allMaterials() );
            }

            Task * task = scinew Task ( "Application::task_time_advance_subproblems_foreign", this, &Application::task_time_advance_subproblems_foreign );
            task->setType ( Task::OncePerProc );
            scheduler->addTask ( task, perProcPatchSet, this->m_materialManager->allMaterials() );
        }
    }

protected: // TASKS

    /**
     * @brief Initialize subproblems task (indexed implementation)
     *
     * Create the SubProblems for each patch/material and save t to dw_new
     * The index sequence is used to flatten the boundary labels tuple
     *
     * @param unused to allow template argument deduction
     * @param myworld data structure to manage mpi processes
     * @param patches set of patches to be initialized
     * @param matls set of materials to be initialized
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    template < size_t... I >
    void
    task_initialize_subproblems (
        index_sequence<I...> _DOXYARG ( unused ),
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * _DOXYARG ( dw_old ),
        DataWarehouse * dw_new
    )
    {
        int myrank = myworld->myRank();
        DOUT ( m_dbg_lvl1, myrank << "==== Application::task_initialize_subproblems ====" );
        for ( int m = 0; m < matls->size(); m++ )
        {
            const int material = matls->get ( m );
            for ( int p = 0; p < patches->size(); ++p )
            {
                const Patch * patch = patches->get ( p );
                DOUT ( m_dbg_lvl2, myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() << " Material: " << material );

                SubProblems < Problem<VAR, STN, Field...> > subproblems ( std::get<I> ( *m_boundary_labels )..., m_subproblems_label, material, patch, m_c2f );
                dw_new->put ( subproblems, m_subproblems_label, material, patch );
            }
        }

        DOUT ( m_dbg_lvl2, myrank );
    }

    /**
     * @brief Initialize subproblems task
     *
     * Create the SubProblems for each patch/material and save it to dw_new
     * Calls the indexed implementation
     *
     * @param myworld data structure to manage mpi processes
     * @param patches set of patches to be initialized
     * @param matls set of materials to be initialized
     * @param dw_old unused
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_initialize_subproblems (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    )
    {
        task_initialize_subproblems ( index_sequence_for<Field...> {}, myworld, patches, matls, dw_old, dw_new );
    }

    /**
     * @brief Advance subproblems task
     *
     * Move SubProblems for each one of the patches and from dw_old to dw_new
     *
     * @param myworld data structure to manage mpi processes
     * @param patches set of patches to be initialized
     * @param matls set of materials to be initialized
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_subproblems (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    )
    {
        int myrank = myworld->myRank();
        const Level * level = getLevel ( patches );
        DOUT ( m_dbg_lvl1, myrank << "==== Application::task_time_advance_subproblems ====" );
        DOUT ( m_dbg_lvl2, myrank << "== Transfer From OldDW Patches: " << *patches << " Level: " << level->getIndex() << " Materials: " << *matls );
        dw_new->transferFrom ( dw_old, m_subproblems_label, patches, matls );

        DOUT ( m_dbg_lvl2, myrank );
    }

    /**
     * @brief Advance subproblems marked as foreign variables task
     *
     * Move SubProblems that are markes as foreign variables and from dw_old to dw_new
     *
     * @remark Foreign SubProblems' exist only when more than one mpi rank is used
     *
     * @param myworld data structure to manage mpi processes
     * @param patches unused
     * @param matls unused
     * @param dw_old DataWarehouse for previous timestep
     * @param dw_new DataWarehouse to be initialized
     */
    void
    task_time_advance_subproblems_foreign (
        const ProcessorGroup * myworld,
        const PatchSubset * _DOXYARG ( patches ),
        const MaterialSubset * _DOXYARG ( matls ),
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    )
    {
        int myrank = myworld->myRank();
        if ( myworld->nRanks() > 1 ) // look for foreign subproblems
        {
            DOUT ( m_dbg_lvl2, myrank << "== Transfer Foreign Vars From OldDW: " );
            dw_new->transferForeignFrom ( dw_old, m_subproblems_label );
        }

        DOUT ( m_dbg_lvl2, myrank );
    }

}; // class Application

/**
 * @brief Virtual base for PhaseField applications (AMR implementation)
 *
 * Wrapper of ApplicationCommon and interfaces
 *
 * @implements Application < Problem, AMR >
 *
 * @tparam VAR type of variable representation
 * @tparam STN finite-difference stencil
 */
template < VarType VAR, StnType STN, typename... Field >
class Application < Problem<VAR, STN, Field...>, true >
    : public Application < Problem<VAR, STN, Field...>, false >
    , public AMRInterface < VAR, get_stn<STN>::dim >
{
protected: // STATIC MEMBERS

    /// If SubProblems are required (i.e. if Problem has any boundary Field)
    using Application< Problem<VAR, STN, Field...>, false >::use_subprblems;

    /// Number of ghost elements required by STN (on the same level)
    using Application< Problem<VAR, STN, Field...> >::FGN;

    /// Type of ghost elements required by VAR and STN (on the same level)
    using Application< Problem<VAR, STN, Field...> >::FGT;

    /// Number of ghost elements required by STN (on the coarser level)
    using Application< Problem<VAR, STN, Field...> >::CGN;

    /// Type of ghost elements required by VAR and STN (on the coarser level)
    using Application< Problem<VAR, STN, Field...> >::CGT;

public: // CONSTRUCTORS/DESTRUCTOR

    /// Constructor
    using Application< Problem<VAR, STN, Field...>, false >::Application;

protected: // SCHEDULINGS

    /**
     * @brief Schedule task_initialize_subproblems and task_communicate_subproblems after regridding
     *
     * Defines the dependencies and output of the task which initializes the
     * subproblems allowing sched to control its execution order
     *
     * @remark subproblems need to be reinitialized on all patches because
     * even preexisting patches may have different neighbors
     * @remark task_communicate_subproblems is scheduled just for inducing MPI
     * send of reinitialized subproblems to other processors that may need them
     *
     * @param grid grid to be initialized
     * @param perProcPatchSet set of all patches assigned to process
     * @param scheduler scheduler to manage the tasks
     */
    virtual void
    scheduleRefineSystemVars (
        const GridP & grid,
        const PatchSet * perProcPatchSet,
        SchedulerP & scheduler
    ) override
    {
        Application< Problem<VAR, STN, Field...>, false >::scheduleRefineSystemVars ( grid, perProcPatchSet, scheduler );

        if ( use_subprblems )
        {
            DOUTR ( dbg_scheduling, "scheduleRefineSystemVars" );

            for ( int idx = 0; idx < grid->numLevels(); ++idx )
            {
                Task * task = scinew Task ( "Application::task_initialize_subproblems", this, &Application::task_initialize_subproblems );
                task->computes ( this->m_subproblems_label );
                scheduler->addTask ( task, grid->getLevel ( idx )->eachPatch(), this->m_materialManager->allMaterials() );
            }
            for ( int idx = 0; idx < grid->numLevels(); ++idx )
            {
                Task * task = scinew Task ( "Application::task_communicate_subproblems", this, &Application::task_communicate_subproblems );
                task->requires ( Task::NewDW, this->m_subproblems_label, FGT, FGN );
                if (idx) task->requires ( Task::NewDW, this->m_subproblems_label, nullptr, Task::CoarseLevel, nullptr, Task::NormalDomain, CGT, CGN );
                task->modifies ( this->m_subproblems_label );
                scheduler->addTask ( task, grid->getLevel ( idx )->eachPatch(), this->m_materialManager->allMaterials() );
            }
        }

#ifdef HAVE_HYPRE
        if ( this->m_solver )
        {
            if ( this->m_solver->getName() == "hypre" )
            {
                scheduler->setRestartInitTimestep(true);
            }
            else if ( this->m_solver->getName() == "hyprefac" )
            {
                for ( int idx = 0; idx < grid->numLevels(); ++idx )
                {
                    this->m_solver->scheduleInitialize ( grid->getLevel ( idx ), scheduler, this->m_materialManager->allMaterials() );
                }
            }
        }
#endif
    }

protected: // TASKS

    /**
     * @brief MPI communicate subproblems task
     *
     * Does nothing
     *
     * @param myworld unused
     * @param patches unused
     * @param matls unused
     * @param dw_old unused
     * @param dw_new unused
     */
    void
    task_communicate_subproblems (
        const ProcessorGroup * _DOXYARG ( myworld ),
        const PatchSubset * _DOXYARG ( patches ),
        const MaterialSubset * _DOXYARG ( matls ),
        DataWarehouse * _DOXYARG ( dw_old ),
        DataWarehouse * _DOXYARG ( dw_new )
    )
    {}

}; // class Application

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_Applications_Application_h
