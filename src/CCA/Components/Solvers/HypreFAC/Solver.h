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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreFAC_Solver_h
#define Packages_Uintah_CCA_Components_Solvers_HypreFAC_Solver_h

// #define PRINTSYSTEM

#include <CCA/Components/Solvers/HypreFAC/AdditionalEntriesP.h>
#include <CCA/Components/Solvers/HypreFAC/SolverStructP.h>
#include <CCA/Components/Solvers/HypreFAC/GlobalDataP.h>
#include <CCA/Components/Solvers/HypreFAC/AdditionalEntries.h>
#include <CCA/Components/Solvers/HypreFAC/SolverStruct.h>
#include <CCA/Components/Solvers/HypreFAC/GlobalData.h>
#include <CCA/Components/Solvers/HypreFAC/SolverParams.h>
#include <CCA/Components/Solvers/HypreFAC/PartData.h>
#include <CCA/Components/Solvers/SolverCommon.h>
#include <Core/Grid/Variables/PerPatch.h> // must be included after GlobalDataP/SolverStructP/AdditionalEntriesP where swapbytes override is defined
#include <Core/Grid/Variables/SoleVariable.h> // must be included after GlobalDataP/SolverStructP/AdditionalEntriesP where swapbytes override is defined
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Util/Timers/Timers.hpp>
#include <Core/Exceptions/ConvergenceFailure.h>

#include <HYPRE_sstruct_ls.h>
#include <HYPRE_krylov.h>

#if 1
#define HYPRE(fn) HYPRE_##fn
#else
#define HYPREDBG_NDIM 3
#include "hypre_dbg.hpp"
#define HYPRE(fn) HYPREDBG_##fn
#endif

#include <string>
#include <cstddef>
#include <numeric>

/**
 *  @class  HypreFAC::Solver
 *  @author Jon Matteo Church
 *  @brief  Uintah hypre solver interface.
 *  Allows the solution of a linear system of the form \[ \mathbf{A} \mathbf{x} = \mathbf{b}\] where \[\mathbf{A}\] is
 *  stencil7 matrix.
 */

namespace Uintah
{
namespace HypreFAC
{

static constexpr bool dbg_doing = false;
static constexpr bool dbg_assembling = false;
extern DebugStream cout_doing;

enum SolverType { FAC, PFMG_PCG, PFMG };

template<SolverType> class SolverImpl;

template<>
class SolverImpl<FAC> : public SolverBase
{
    HYPRE_SStructSolver    solver;

public:
    SolverImpl (
        MPI_Comm comm,
        SolverStruct * solver_struct,
        SolverParams * solver_params,
        int * plevels,
        int ( * prefinements ) [HYPRE_MAXDIM],
        bool guess
    )
    {
        const auto & nparts = solver_struct->data->nparts;

        HYPRE(SStructFACCreate) ( comm, &solver );
        HYPRE(SStructFACSetPLevels) ( solver, nparts, plevels );
        HYPRE(SStructFACSetPRefinements) ( solver, nparts, prefinements );

        if ( guess )
            HYPRE(SStructFACSetNonZeroGuess) ( solver );
        else
            HYPRE(SStructFACSetZeroGuess) ( solver );

        if ( solver_params->max_levels > 0 )
            HYPRE(SStructFACSetMaxLevels) ( solver, solver_params->max_levels );
        if ( solver_params->max_iter > 0 )
            HYPRE(SStructFACSetMaxIter) ( solver, solver_params->max_iter );
        if ( solver_params->tol > 0 )
            HYPRE(SStructFACSetTol) ( solver, solver_params->tol );
        if ( solver_params->rel_change > -1 )
            HYPRE(SStructFACSetRelChange) ( solver, solver_params->rel_change );
        if ( solver_params->relax_type >  -1 )
            HYPRE(SStructFACSetRelaxType) ( solver, solver_params->relax_type );
        if ( solver_params->num_pre_relax >  -1 )
            HYPRE(SStructFACSetNumPreRelax) ( solver, solver_params->num_pre_relax );
        if ( solver_params->num_post_relax >  -1 )
            HYPRE(SStructFACSetNumPostRelax) ( solver, solver_params->num_post_relax );
        if ( solver_params->csolver_type >  -1 )
            HYPRE(SStructFACSetCoarseSolverType) ( solver, solver_params->csolver_type );
        if ( solver_params->logging >  -1 )
            HYPRE(SStructFACSetLogging) ( solver, solver_params->logging );

        HYPRE(SStructFACSetup2) ( solver,
                                     solver_struct->A,
                                     solver_struct->b,
                                     solver_struct->x
                                   );
    }

    virtual ~SolverImpl()
    {
        HYPRE(SStructFACDestroy2) ( solver );
    };

    virtual void
    solve (
        SolverStruct * solver_struct,
        int & num_iterations, double & final_res_norm
    ) override
    {
        HYPRE(SStructFACSolve3) (
            solver,
            solver_struct->A,
            solver_struct->b,
            solver_struct->x
        );
        HYPRE(SStructFACGetNumIterations) ( solver, &num_iterations );
        HYPRE(SStructFACGetFinalRelativeResidualNorm) ( solver, &final_res_norm );
    }

    virtual void destroy () override
    {
        ASSERTFAIL ( "hehe" );
    }
};

template<>
class SolverImpl<PFMG_PCG> : public SolverBase
{
    HYPRE_SStructSolver    solver;
    HYPRE_SStructSolver    precond;

public:
    SolverImpl (
        MPI_Comm comm,
        SolverStruct * solver_struct,
        SolverParams * solver_params//,
//         int * plevels,
//         int ( * prefinements ) [HYPRE_MAXDIM],
//         bool guess
    )
    {
        HYPRE_SStructPCGCreate ( comm, &solver );
        HYPRE_PCGSetTwoNorm ( ( HYPRE_Solver ) solver, 1 );

        if ( solver_params->max_iter > 0 )
            HYPRE_PCGSetMaxIter ( ( HYPRE_Solver ) solver, solver_params->max_iter );
        if ( solver_params->tol > 0 )
            HYPRE_PCGSetTol ( ( HYPRE_Solver ) solver, solver_params->tol );

        /* use SysPFMG solver as preconditioner */
        HYPRE_SStructSysPFMGCreate ( comm, &precond );
        HYPRE_SStructSysPFMGSetMaxIter ( precond, 1 );
        HYPRE_SStructSysPFMGSetTol ( precond, 0.0 );
        HYPRE_SStructSysPFMGSetZeroGuess ( precond );
        /* weighted Jacobi = 1; red-black GS = 2 */
        HYPRE_SStructSysPFMGSetRelaxType ( precond, 3 );
        if ( solver_params->relax_type == WeightedJacobi )
            HYPRE_SStructFACSetJacobiWeight ( precond, solver_params->weight );
        if ( solver_params->num_pre_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPreRelax ( precond, solver_params->num_pre_relax );
        if ( solver_params->num_post_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPostRelax ( precond, solver_params->num_post_relax );

        HYPRE_PCGSetPrecond ( ( HYPRE_Solver ) solver,
                              ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSolve,
                              ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSetup,
                              ( HYPRE_Solver ) precond );

        HYPRE_PCGSetup ( ( HYPRE_Solver ) solver,
                         ( HYPRE_Matrix ) solver_struct->A,
                         ( HYPRE_Vector ) solver_struct->b,
                         ( HYPRE_Vector ) solver_struct->x );
    }

    virtual ~SolverImpl()
    {
        HYPRE_SStructPCGDestroy ( solver );
        HYPRE_SStructSysPFMGDestroy ( precond );
    };

    virtual void
    solve (
        SolverStruct * solver_struct,
        int & num_iterations,
        double & final_res_norm
    ) override
    {
        HYPRE_PCGSolve ( ( HYPRE_Solver ) solver,
                         ( HYPRE_Matrix ) solver_struct->A,
                         ( HYPRE_Vector ) solver_struct->b,
                         ( HYPRE_Vector ) solver_struct->x );
        HYPRE_PCGGetNumIterations ( ( HYPRE_Solver ) solver, &num_iterations );
        HYPRE_PCGGetFinalRelativeResidualNorm ( ( HYPRE_Solver ) solver, &final_res_norm );
    }

    virtual void destroy () override
    {
        ASSERTFAIL ( "hehe" );
    }

};

template<>
class SolverImpl<PFMG> : public SolverBase
{
    HYPRE_SStructSolver    solver;
    int max_iter;

public:
    SolverImpl (
        MPI_Comm comm,
        SolverStruct * solver_struct,
        SolverParams * solver_params,
//         int * plevels,
//         int ( * prefinements ) [HYPRE_MAXDIM],
        bool guess
    ) : max_iter (solver_params->max_iter)
    {
        HYPRE_SStructSysPFMGCreate ( comm, &solver );
        if ( solver_params->max_iter > 0 )
            HYPRE_SStructSysPFMGSetMaxIter ( solver, solver_params->max_iter + 1 );
        if ( solver_params->tol > 0 )
            HYPRE_SStructSysPFMGSetTol ( solver, solver_params->tol );
        if ( guess )
            HYPRE_SStructSysPFMGSetNonZeroGuess ( solver );
        else
            HYPRE_SStructSysPFMGSetZeroGuess ( solver );
        /* weighted Jacobi = 1; red-black GS = 2 */
        HYPRE_SStructSysPFMGSetRelaxType ( solver, solver_params->relax_type );
        if ( solver_params->relax_type == WeightedJacobi )
            HYPRE_SStructFACSetJacobiWeight ( solver, solver_params->weight );
        if ( solver_params->num_pre_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPreRelax ( solver, solver_params->num_pre_relax );
        if ( solver_params->num_post_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPostRelax ( solver, solver_params->num_post_relax );

        HYPRE_SStructSysPFMGSetPrintLevel ( solver, 1 );

        HYPRE_SStructSysPFMGSetup (
            solver,
            solver_struct->A,
            solver_struct->b,
            solver_struct->x
        );
    }

    virtual ~SolverImpl()
    {
        HYPRE_SStructSysPFMGDestroy ( solver );
    };

    virtual void 
    solve (
        SolverStruct * solver_struct,
        int & num_iterations,
        double & final_res_norm
    ) override
    {
        HYPRE_SStructSysPFMGSolve (
            solver,
            solver_struct->A,
            solver_struct->b,
            solver_struct->x
        );
        HYPRE_SStructSysPFMGGetNumIterations ( solver, &num_iterations );
        if ( max_iter > 0 && num_iterations <= max_iter )
            final_res_norm = 0.;
        HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm ( solver, &final_res_norm );  // BROKEN !
    }

    virtual void destroy () override
    {
        ASSERTFAIL ( "hehe" );
    }

};


template < size_t DIM >
class Solver : public SolverCommon
{
public:
    static std::string AdditionalEntriesSuffix;

private:
    static constexpr size_t stencil_size = 2 * DIM + 1;
    static constexpr Patch::FaceType face_start = ( Patch::FaceType ) ( 0 );
    static constexpr Patch::FaceType face_end = ( Patch::FaceType ) ( 2 * DIM - 1 );

    static int offsets[stencil_size][DIM];
    static int stn0, stnE, stnW, stnN, stnS, stnT, stnB;

    const VarLabel * m_timeStepLabel;
    const VarLabel * m_hypre_global_data_label;
    std::map <int, const VarLabel *> m_hypre_solver_struct_label;

    SolverParams * m_params;

    // setup in scheduleSolve
    bool m_modifies_solution;
    const VarLabel * m_stencil_entries_label;
    const VarLabel * m_additional_entries_label;
    Task::WhichDW m_matrix_dw;
    const VarLabel * m_solution_label;
    const VarLabel * m_rhs_label;
    Task::WhichDW  m_rhs_dw;
    const VarLabel * m_guess_label;
    Task::WhichDW m_guess_dw;

    // hypre timers - note that these variables do NOT store timings - rather, each corresponds to
    // a different timer index that is managed by Hypre. To enable the use and reporting of these
    // hypre timings, #define HYPRE_TIMING in HypreSolver.h
    int tHypreAll_;    // Tracks overall time spent in Hypre = matrix/vector setup & assembly + solve time.
    int tSolveOnly_;   // Tracks time taken by hypre to solve the system of equations
    int tMatVecSetup_; // Tracks the time taken by uintah/hypre to allocate and set matrix and vector box vaules

    double m_moving_average;

    void
    add_fine_to_fine_entry
    (
        const int & part,
        const IntVector & index,
        const IntVector & to_index,
        const double & value,
        CCVariable<Stencil7> & stencil_entries,
        AdditionalEntries & additional_entries
    )
    {
        // check if diagonal entry
        if ( index == to_index )
        {
            stencil_entries[index].p += value;
            return;
        }

        // check if stencil entry
        for ( size_t entry = 1; entry < stencil_size; ++entry )
        {
            bool is_stencil = true;
            for ( size_t d = 0; d < DIM; ++d )
                if ( to_index[d] - index[d] != offsets[entry][d] )
                {
                    is_stencil = false;
                    break;
                }

            if ( is_stencil )
            {
                stencil_entries[index][entry - 1] += value;
                return;
            }
        }

        // otherwise is additional entry
        MatrixEntry::first_type key = {index, part, to_index};
        auto res = additional_entries.emplace ( key, value );
        if ( !res.second ) res.first->second += value;
    }

    template <size_t NDIM>
    typename std::enable_if< NDIM == 2ul, std::vector<double> >::type
    get_stencil_values (
        const PartDataP & pdata,
        const int & box,
        const CCVariable<Stencil7> & stencil_entries
    )
    {
        std::vector<double> values ( pdata->boxsizes[box] * stencil_size );

        std::vector<double>::iterator it = values.begin();
        int k = pdata->ilowers[box][2];
        for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
            for ( int i = pdata->ilowers[box][0]; i <= pdata->iuppers[box][0]; ++i )
            {
                const auto & vals = stencil_entries[IntVector ( i, j, k )];
                ( *it++ ) = vals.p;
                ( *it++ ) = vals.w;
                ( *it++ ) = vals.e;
                ( *it++ ) = vals.s;
                ( *it++ ) = vals.n;
            }

        return values;
    }

    template <size_t NDIM>
    typename std::enable_if< NDIM == 3ul, std::vector<double> >::type
    get_stencil_values (
        const PartDataP & pdata,
        const int & box,
        const CCVariable<Stencil7> & stencil_entries
    )
    {
        std::vector<double> values ( pdata->boxsizes[box] * stencil_size );

        std::vector<double>::iterator it = values.begin();
        for ( int k = pdata->ilowers[box][2]; k <= pdata->iuppers[box][2]; ++k )
            for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
                for ( int i = pdata->ilowers[box][0]; i <= pdata->iuppers[box][0]; ++i )
                {
                    const auto & vals = stencil_entries[IntVector ( i, j, k )];
                    ( *it++ ) = vals.p;
                    ( *it++ ) = vals.w;
                    ( *it++ ) = vals.e;
                    ( *it++ ) = vals.s;
                    ( *it++ ) = vals.n;
                    ( *it++ ) = vals.b;
                    ( *it++ ) = vals.t;
                }

        return values;
    }

public:
    Solver (
        const ProcessorGroup * myworld
    )
        : SolverCommon ( myworld )
        , m_timeStepLabel ( VarLabel::create ( timeStep_name, timeStep_vartype::getTypeDescription() ) )
        , m_hypre_global_data_label ( VarLabel::create ( "hypre_global_data", SoleVariable<GlobalDataP>::getTypeDescription() ) )
        , m_hypre_solver_struct_label ()
        , m_params ( scinew SolverParams() )
        , m_moving_average ( 0. )
    {
    }

    virtual ~Solver()
    {
        VarLabel::destroy ( m_timeStepLabel );
        VarLabel::destroy ( m_hypre_global_data_label );
        for ( auto && l : m_hypre_solver_struct_label ) VarLabel::destroy ( l.second );
    }

    virtual
    void
    readParameters
    (
        ProblemSpecP & params_ps,
        const std::string  & varname
    ) override
    {
        bool found = false;
        if ( params_ps )
        {
            for ( ProblemSpecP param_ps = params_ps->findBlock ( "Parameters" ); param_ps != nullptr; param_ps = param_ps->findNextBlock ( "Parameters" ) )
            {
                std::string variable;
                if ( param_ps->getAttribute ( "variable", variable ) && variable != varname )
                {
                    continue;
                }
                int sFreq;
                int coefFreq;
                param_ps->getWithDefault ( "solveFrequency",      m_params->solveFrequency, 1 );
                param_ps->getWithDefault ( "setupFrequency",      sFreq,                    1 );
                param_ps->getWithDefault ( "updateCoefFrequency", coefFreq,                 1 );
                m_params->setSetupFrequency ( sFreq );
                m_params->setUpdateCoefFrequency ( coefFreq );

                ASSERTMSG ( param_ps->get ( "max_levels",             m_params->max_levels ),   "HypreFAC::Solver ERROR. Missing parameter: max_level" );
                param_ps->getWithDefault ( "maxiterations",          m_params->max_iter,       -1 );
                param_ps->getWithDefault ( "tolerance",              m_params->tol,            -1. );
                param_ps->getWithDefault ( "rel_change",             m_params->rel_change,     -1 );
                param_ps->getWithDefault ( "relax_type", ( int & ) m_params->relax_type,     -1 );
                param_ps->getWithDefault ( "weight",                 m_params->weight,         -1. );
                param_ps->getWithDefault ( "npre",                   m_params->num_pre_relax,  -1 );
                param_ps->getWithDefault ( "npost",                  m_params->num_post_relax, -1 );
                param_ps->getWithDefault ( "csolver_type", ( int & ) m_params->csolver_type,   -1 );
                param_ps->getWithDefault ( "logging",                m_params->logging,        -1 );
                param_ps->getWithDefault ( "C2F",                    m_params->C2F,            -1 );

                found = true;
            }
        }
        if ( !found )
        {
            m_params->solveFrequency = 1;
            m_params->setSetupFrequency ( 1 );
            m_params->setUpdateCoefFrequency ( 1 );

            m_params->max_levels = -1;
            m_params->tol = -1.;
            m_params->max_iter = -1;
            m_params->rel_change = -1;
            m_params->relax_type = DefaultRelaxType;
            m_params->weight = -1.;
            m_params->num_pre_relax = -1;
            m_params->num_post_relax = -1;
            m_params->logging = -1;
            m_params->C2F = -1;
        }

        if ( m_params->C2F == -1 )
            m_params->C2F = 1; // FIXME 0 (should be NC default), 1 (should be CC default)
    }

    virtual inline
    SolverParameters *
    getParameters()
    override
    {
        return m_params;
    }

    /**
     *  @brief Schedules the solution of the linear system \[ \mathbf{A} \mathbf{x} = \mathbf{b}\].
     *
     *  @param level A reference to the level on which the system is solved.
     *  @param scheduler A reference to the Uintah scheduler.
     *  @param materials A pointer to the MaterialSet.
     *  @param matrix_label Varlabel of the coefficient matrix \[\mathbf{A}\]
     *  @param matrix_dw The datawarehouse in which the coefficient matrix lives.
     *  @param solution_label The varlabel of the solutio1n vector.
     *  @param modifies_solution A boolean that specifies the behaviour of the task
                          associated with the ScheduleSolve. If set to true,
                          then the task will only modify x. Otherwise, it will
                          compute x. This is a key option when you are computing
                          x in another place.
     * @param rhs_label The VarLabel of the right hand side vector.
     * @param rhs_dw Specifies the datawarehouse in which b lives.
     * @param guess_label VarLabel of the initial guess.
     * @param guess_dw Specifies the datawarehouse of the initial guess.
     *
     */
    virtual
    void
    scheduleSolve
    (
        const LevelP & level,
        SchedulerP & scheduler,
        const MaterialSet * materials,
        const VarLabel * matrix_label,
        Task::WhichDW matrix_dw,
        const VarLabel * solution_label,
        bool modifies_solution,
        const VarLabel * rhs_label,
        Task::WhichDW rhs_dw,
        const VarLabel * guess_label,
        Task::WhichDW guess_dw,
        bool /*is_first_solve*/ = true
    ) override
    {
        if ( level->hasCoarserLevel() ) return;

        m_modifies_solution = modifies_solution;
        m_stencil_entries_label = matrix_label;
        m_additional_entries_label = VarLabel::find ( matrix_label->getName() + Solver::AdditionalEntriesSuffix );
        m_rhs_label = rhs_label;
        m_rhs_dw = rhs_dw;
        m_matrix_dw = matrix_dw;
        m_solution_label = solution_label;
        m_guess_label = guess_label;
        m_guess_dw = guess_dw;

        printSchedule ( level, cout_doing, "HypreFAC::Solver:scheduleSolve" );

        auto global_data_label = hypre_global_data_label();
        auto solver_struct_label = hypre_solver_struct_label ( d_myworld->myRank() );
        auto grid = level->getGrid();

        Task * task = scinew Task ( "HypreFAC::Solver::solve", this, &Solver::solve, global_data_label, solver_struct_label );
        task->requires ( Task::OldDW, solver_struct_label );
        task->requires ( Task::OldDW, global_data_label );

        task->computes ( global_data_label );
        task->computes ( solver_struct_label );

        task->modifies ( m_stencil_entries_label );
        task->modifies ( m_additional_entries_label );

        task->modifies ( m_rhs_label );
        task->requires ( m_rhs_dw, m_rhs_label, Ghost::None, 0 );
        task->requires ( Task::NewDW, m_timeStepLabel );

        if ( m_guess_label )
        {
            task->requires ( m_guess_dw, m_guess_label, Ghost::None, 0 );
        }

        task->computes ( m_solution_label );

        task->setType ( Task::OncePerProc );
        task->usesMPI ( true );

        auto * patches = scheduler->getLoadBalancer()->getPerProcessorPatchSet ( grid );
        scheduler->addTask ( task, patches, materials );

        scheduler->overrideVariableBehavior ( m_additional_entries_label->getName(), false, false, false, true, true );

    };

    virtual
    void
    scheduleInitialize
    (
        const LevelP & level,
        SchedulerP & scheduler,
        const MaterialSet * matls
    ) override
    {
        auto global_data_label = hypre_global_data_label();
        auto solver_struct_label = hypre_solver_struct_label ( d_myworld->myRank() );

        printSchedule ( level, cout_doing, "HypreFAC::Solver:scheduleInitialize" );

        Task * task = scinew Task ( "initialize_hypre", this, &Solver::initialize, level );

        task->computes ( global_data_label );
        task->computes ( solver_struct_label );

        task->setType ( Task::OncePerProc );

        ASSERT ( scheduler );
        ASSERT ( scheduler->getLoadBalancer() );
        auto * patches = scheduler->getLoadBalancer()->getPerProcessorPatchSet ( level );
        scheduler->addTask ( task, patches, matls );

        scheduler->overrideVariableBehavior ( global_data_label->getName(), false, false, false, true, true );
        scheduler->overrideVariableBehavior ( solver_struct_label->getName(), false, false, false, true, true );
    };


    virtual inline
    void
    scheduleRestartInitialize
    (
        const LevelP & level,
        SchedulerP & scheduler,
        const MaterialSet * matls
    )
    override
    {
        scheduleInitialize ( level, scheduler, matls );
    };

    virtual
    std::string
    getName()
    override
    {
        return "hyprefac";
    }

    void
    allocateHypreMatrices
    (
        DataWarehouse * new_dw,
        const Level * level
    );

private:
    void initialize
    (
        const ProcessorGroup * pg,
        const PatchSubset * patches,
        const MaterialSubset * /*matls*/,
        DataWarehouse * old_dw,
        DataWarehouse * new_dw,
        LevelP level
    )
    {
        const int nranks = pg->nRanks();
        const int rank = pg->myRank();
        const int nparts = level->getGrid()->numLevels();
        const int part = level->getIndex();
        const int nboxes = level->numPatches();
        const int npatches = patches->size();

        DOUT ( dbg_doing, "\nHypreFAC::Solver::initialize ( rank " << rank << " level " << part << " ) " );

        IntVector low, high;
        level->getGrid()->getLevel ( 0 )->findCellIndexRange ( low, high );
        IntVector periodic = level->getGrid()->getLevel ( 0 )->getPeriodicBoundaries();
        IntVector refinement = level->getGrid()->getLevel ( 0 )->getRefinementRatio();

        for ( int l = 1; l <= part; ++l )
        {
            auto r = level->getGrid()->getLevel ( l )->getRefinementRatio();
            for ( int i = 0; i < HYPRE_MAXDIM; ++i )
                refinement[i] *= r[i];
        }

        IntVector prefinement = level->getRefinementRatio();

        // time step is not there after regrid
        timeStep_vartype timeStep ( 0 );
        if ( !new_dw->exists ( m_timeStepLabel ) )
        {
            old_dw->get ( timeStep, m_timeStepLabel );
            new_dw->put ( timeStep, m_timeStepLabel );
        }
        else
            new_dw->get ( timeStep, m_timeStepLabel );

        SolverStruct * solver_struct;
        auto solver_struct_label = hypre_solver_struct_label ( rank );
        if ( new_dw->exists ( solver_struct_label ) )
        {
            SoleVariable<SolverStructP> solver_struct_var;
            new_dw->get ( solver_struct_var, solver_struct_label );
            solver_struct = solver_struct_var.get().get_rep();
        }
        else
        {
            GlobalData * global_data;
            auto global_data_label = hypre_global_data_label();
            if ( new_dw->exists ( global_data_label ) )
            {
                SoleVariable<GlobalDataP> global_data_var;
                new_dw->get ( global_data_var, global_data_label );
                global_data = global_data_var.get().get_rep();
            }
            else
            {
                global_data = scinew GlobalData;
                global_data->nparts = nparts;
                global_data->nboxes.resize ( nranks );
                for ( auto & nboxes : global_data->nboxes )
                    nboxes.resize ( nparts );
                global_data->nvars = 1;

                SoleVariable<GlobalDataP> global_data_var;
                global_data_var.setData ( global_data );
                new_dw->put ( global_data_var, global_data_label );
            }

            solver_struct = scinew SolverStruct;
            solver_struct->data = global_data;
            solver_struct->pdatas = scinew PartDataP[nparts];
            solver_struct->created = false;
            solver_struct->restart = true;

            SoleVariable<SolverStructP> solver_struct_var;
            solver_struct_var.setData ( solver_struct );
            new_dw->put ( solver_struct_var, solver_struct_label );
        }

        if ( !solver_struct->pdatas[part] )
        {
            solver_struct->data->nboxes[rank][part] = npatches;
            solver_struct->data->vartypes.emplace_back ( HYPRE_SSTRUCT_VARIABLE_CELL );

            PartDataP pdata = solver_struct->pdatas[part] = scinew PartData;
            pdata->nboxes = nboxes;

            pdata->plevel = part;

            pdata->low[0] = low[0] * refinement[0];
            pdata->low[1] = low[1] * refinement[1];
            pdata->low[2] = low[2] * refinement[2];

            pdata->high[0] = high[0] * refinement[0] - 1;
            pdata->high[1] = high[1] * refinement[1] - 1;
            pdata->high[2] = high[2] * refinement[2] - 1;

            pdata->prefinement[0] = prefinement[0];
            pdata->prefinement[1] = prefinement[1];
            pdata->prefinement[2] = prefinement[2];

            pdata->periodic[0] = periodic[0] * ( pdata->high[0] - pdata->low[0] + 1 );
            pdata->periodic[1] = periodic[1] * ( pdata->high[1] - pdata->low[1] + 1 );
            pdata->periodic[2] = periodic[2] * ( pdata->high[2] - pdata->low[2] + 1 );

            for ( int p = 0; p < patches->size(); p++ )
            {
                const Patch * patch = patches->get ( p );
                IntVector ilower = patch->getCellLowIndex();
                IntVector iupper = patch->getCellHighIndex() - IntVector { 1, 1, 1 };
                std::array<bool, 6> interface = {{ false, false, false, false, false, false }};

                int boxsize = 1;
                for ( int i = 0; i < HYPRE_MAXDIM; i++ )
                    boxsize *= ( iupper[i] - ilower[i] + 1 );

                for ( Patch::FaceType f = face_start; f <= face_end; f = Patch::nextFace ( f ) )
                    if ( patch->getBCType ( f ) == Patch::Coarse )
                        interface[f] = true;

                pdata->ilowers.emplace_back ( ilower );
                pdata->iuppers.emplace_back ( iupper );
                pdata->boxsizes.emplace_back ( boxsize );
                pdata->interfaces.emplace_back ( interface );

                int id = patch->getID();
                // std::map::operator[] leaks memory!
                if ( pdata->patch2box.find ( id ) == pdata->patch2box.end() )
                    pdata->patch2box.insert ( std::make_pair ( id, pdata->box2patch.size() ) );
                pdata->box2patch.emplace_back ( id );
            }
        }
    }

    void
    solve
    (
        const ProcessorGroup * pg,
        const PatchSubset * patches,
        const MaterialSubset * matls,
        DataWarehouse * old_dw,
        DataWarehouse * new_dw,
        const VarLabel * global_data_label,
        const VarLabel * solver_struct_label
    )
    {
        auto comm = pg->getComm();
        auto rank = pg->myRank();

        DOUT ( dbg_doing, "\nHypreFAC::Solver::solve (rank " << rank << ")" );

        struct SolverStruct * solver_struct = nullptr;

        ASSERTMSG ( old_dw->exists ( m_timeStepLabel ), "old time step should exist" );
        ASSERTMSG ( !new_dw->exists ( m_timeStepLabel ), "new time step should not exist" );
        timeStep_vartype timeStep ( 0 );
        old_dw->get ( timeStep, m_timeStepLabel );
        new_dw->put ( timeStep, m_timeStepLabel );

        ASSERTMSG ( old_dw->exists ( global_data_label ), "old hypre fac global data should exist" );
        ASSERTMSG ( !new_dw->exists ( global_data_label ), "new hypre fac global data should not exist" );
        SoleVariable<GlobalDataP> global_data_var;
        old_dw->get ( global_data_var, global_data_label );
        new_dw->put ( global_data_var, global_data_label );

        ASSERTMSG ( old_dw->exists ( solver_struct_label ), "old hypre fac solver struct should exist" );
        ASSERTMSG ( !new_dw->exists ( solver_struct_label ), "new hypre fac solver struct should not exist" );
        SoleVariable<SolverStructP> solver_struct_var;
        old_dw->get ( solver_struct_var, solver_struct_label );
        new_dw->put ( solver_struct_var, solver_struct_label );
        solver_struct = solver_struct_var.get().get_rep();

        bool do_setup, do_update, do_solve;
        //________________________________________________________
        // Matrix setup frequency - this will destroy and recreate a new Hypre matrix at the specified setupFrequency
        {
            int freq = m_params->getSetupFrequency();
            do_setup = ! ( freq == 0 || timeStep % freq ) || ( timeStep == 1 ) || ( m_application->isRegridTimeStep() );
        }

        //________________________________________________________
        // Update coefficient frequency - This will ONLY UPDATE the matrix coefficients without destroying/recreating the Hypre Matrix
        {
            int freq = m_params->getUpdateCoefFrequency();
            do_update = do_setup || ! ( freq == 0 || timeStep % freq ) || ( timeStep == 1 );
        }

        //________________________________________________________
        // Solve frequency
        {
            int freq = m_params->solveFrequency;
            do_solve = ! ( freq == 0 || timeStep % freq );
            // note - the first timeStep in hypre is timeStep 1
        }

        if ( !do_solve )
        {
            new_dw->transferFrom ( old_dw, m_solution_label, patches, matls, true );
            return;
        }

#ifdef HYPRE_TIMING
        tHypreAll_ = hypre_InitializeTiming ( "Total Hypre time" );
        hypre_BeginTiming ( tHypreAll_ );

        tMatVecSetup_ = hypre_InitializeTiming ( "Matrix + Vector setup" );
        tSolveOnly_   = hypre_InitializeTiming ( "Solve time" );
#endif

        DataWarehouse * matrix_dw = new_dw->getOtherDataWarehouse ( m_matrix_dw );
        DataWarehouse * rhs_dw = new_dw->getOtherDataWarehouse ( m_rhs_dw );
        DataWarehouse * guess_dw = new_dw->getOtherDataWarehouse ( m_guess_dw );

        ASSERTEQ ( sizeof ( Stencil7 ), 7 * sizeof ( double ) );

        Timers::Simple timer;
        timer.start();

        bool restart = solver_struct->restart;
        solver_struct->restart = false;

        for ( int m = 0; m < matls->size(); ++m )
        {
            int material = matls->get ( m );

#ifdef HYPRE_TIMING
            hypre_BeginTiming ( tMatVecSetup_ );
#endif

            const auto & nparts = solver_struct->data->nparts;
            const auto & nvars = solver_struct->data->nvars;
            auto & vartypes = solver_struct->data->vartypes;
            auto & nboxes = solver_struct->data->nboxes;

            int * plevels = nullptr;
            int ( *prefinements ) [HYPRE_MAXDIM] = nullptr;

            std::vector<std::vector<std::map<IntVector, std::vector<double>>>> extra_values ( nvars );
            for ( auto & v : extra_values )
                v.resize ( nparts );

            if ( nparts > 1 )
            {
                plevels = new int[nparts];
                prefinements = new int[nparts][HYPRE_MAXDIM];
                for ( int part = 0; part < nparts; ++part )
                {
                    auto pdata = solver_struct->pdatas[part];
                    plevels[part] = pdata->plevel;
                    prefinements[part][0] = pdata->prefinement[0];
                    prefinements[part][1] = pdata->prefinement[1];
                    prefinements[part][2] = pdata->prefinement[2];
                }
            }

            if ( do_setup || restart )
            {
//                 if ( solver_struct->created )
//                 {
//                     HYPRE(SStructVectorDestroy) ( solver_struct->x );
//                     HYPRE(SStructVectorDestroy) ( solver_struct->b );
//                     HYPRE(SStructMatrixDestroy) ( solver_struct->A );
//                     HYPRE(SStructGraphDestroy) ( solver_struct->graph );
//                     HYPRE(SStructStencilDestroy) ( solver_struct->stencil );
//                     HYPRE(SStructGridDestroy) ( solver_struct->grid );
//                     solver_struct->destroy_solver();
//                     solver_struct->created = false;
//                 }

                /*-----------------------------------------------------------
                 * Set up the grid
                 *-----------------------------------------------------------*/

                DOUT ( dbg_assembling, rank << "   hypre grid setup" );

                HYPRE(SStructGridCreate) ( comm, DIM, nparts, &solver_struct->grid );

                for ( int part = 0; part < nparts; ++part )
                {
                    auto & pdata = solver_struct->pdatas[part];
                    for ( int box = 0; box < nboxes[rank][part]; ++box )
                        HYPRE(SStructGridSetExtents) ( solver_struct->grid, part, pdata->ilowers[box].get_pointer(), pdata->iuppers[box].get_pointer() );
                    HYPRE(SStructGridSetVariables) ( solver_struct->grid, part, nvars, vartypes.data() );
                    HYPRE(SStructGridSetPeriodic) ( solver_struct->grid, part, pdata->periodic );
                }
                HYPRE(SStructGridAssemble) ( solver_struct->grid );

                /*-----------------------------------------------------------
                 * Set up the stencils
                 *-----------------------------------------------------------*/

                DOUT ( dbg_assembling, rank << "   hypre stencil setup" );

                HYPRE(SStructStencilCreate) ( DIM, stencil_size, &solver_struct->stencil );
                for ( int var = 0; var < nvars; ++var )
                    for ( size_t entry = 0; entry < stencil_size; ++entry )
                        HYPRE(SStructStencilSetEntry) ( solver_struct->stencil, entry, offsets[entry], var );

                /*-----------------------------------------------------------
                 * Set up the graph
                 *-----------------------------------------------------------*/

                DOUT ( dbg_assembling, rank << "   hypre graph setup" );

                HYPRE(SStructGraphCreate) ( comm, solver_struct->grid, &solver_struct->graph );
                HYPRE(SStructGraphSetObjectType) ( solver_struct->graph, HYPRE_SSTRUCT );

                /* set stencils */
                for ( int part = 0; part < nparts; ++part )
                    for ( int var = 0; var < nvars; ++var )
                        HYPRE(SStructGraphSetStencil) ( solver_struct->graph, part, var, solver_struct->stencil );

                if (patches->size())
                {
                GridP grd = patches->get ( 0 )->getLevel()->getGrid();

                AdditionalEntries *** additional_var = scinew AdditionalEntries ** [nparts];
                additional_var[0] = nullptr;
                for ( int fine_part = 1; fine_part < nparts; ++fine_part )
                {
                    additional_var[fine_part] = scinew AdditionalEntries * [nboxes[rank][fine_part]];
                    auto & pdata = solver_struct->pdatas[fine_part];
                    for ( int fine_box = 0; fine_box < nboxes[rank][fine_part]; ++fine_box )
                    {
                        PerPatch<AdditionalEntriesP> additional_entries;
                        const Patch * patch = grd->getPatchByID ( pdata->box2patch[fine_box], fine_part );
                        matrix_dw->get ( additional_entries, m_additional_entries_label, material, patch );
                        additional_var[fine_part][fine_box] = additional_entries.get().get_rep();
                    }
                }


#if 1
                /* add coarse to fine entries */
                for ( int fine_part = 1; fine_part < nparts; ++fine_part )
                {
                    auto & pdata = solver_struct->pdatas[fine_part];
                    int coarse_part = fine_part - 1;

                    std::function<void ( IntVector &, IntVector &, const int &, const int &, const int &, double & ) > refine_interface;

                    switch ( m_params->C2F )
                    {
                    case 0:  // interpolate coarse node with fine node
                    {
                        refine_interface = [&solver_struct, coarse_part, fine_part, &extra_values] (
                                               IntVector & coarse_index,
                                               IntVector & fine_index,
                                               const int & var,
                                               const int & d,
                                               const int & s,
                                               double & stencil_entry
                                           )
                        {
                            if ( stencil_entry )
                            {
                                IntVector to_index ( fine_index[0], fine_index[1], fine_index[2] );
                                HYPRE(SStructGraphAddEntries) ( solver_struct->graph, coarse_part, coarse_index.get_pointer(), var, fine_part, to_index.get_pointer(), var );
                                extra_values[var][coarse_part][coarse_index].emplace_back ( stencil_entry );
                                stencil_entry = 0.;
                            }
                        };
                    }
                    break;
                    case 1:
                    {
                        double n_fine = pdata->prefinement[0] * pdata->prefinement[1] * pdata->prefinement[2];
                        refine_interface = [&solver_struct, coarse_part, fine_part, pdata, n_fine, &extra_values] (
                                               IntVector & coarse_index,
                                               IntVector & fine_index,
                                               const int & var,
                                               const int & d,
                                               const int & s,
                                               double & stencil_entry
                                           )
                        {
                            if ( stencil_entry )
                            {
                                for ( int k = 0; k < pdata->prefinement[2]; ++k )
                                    for ( int j = 0; j < pdata->prefinement[1]; ++j )
                                        for ( int i = 0; i < pdata->prefinement[0]; ++i )
                                        {
                                            IntVector to_index ( fine_index[0] + i, fine_index[1] + j, fine_index[2] + k );
                                            if ( s > 0 ) to_index[d] -= pdata->prefinement[d] - 1;
                                            HYPRE(SStructGraphAddEntries) ( solver_struct->graph, coarse_part, coarse_index.get_pointer(), var, fine_part, to_index.get_pointer(), var );
                                            extra_values[var][coarse_part][coarse_index].emplace_back ( stencil_entry / n_fine );
                                        }
                                stencil_entry = 0.;
                            }
                        };
                    }
                    break;
                    default:
                        ASSERTFAIL ( "C2F value not handeld" );
                    }

                    // add zero connection between fine patch and refined coarse patch
                    // fac setup will otherwise fail when forming the coarse operators
                    // since it reassign the coarse box to the processors owning ts refined boxes
                    for ( int coarse_box = 0; coarse_box < nboxes[rank][coarse_part]; ++coarse_box )
                    {
                        auto & pdata = solver_struct->pdatas[coarse_part];
                        const Patch * coarse_patch = grd->getPatchByID ( pdata->box2patch[coarse_box], coarse_part );
                        const Level * coarse_level = coarse_patch->getLevel();
                        const Level * fine_level = coarse_level->getFinerLevel().get_rep();
                        int * prefinement = solver_struct->pdatas[fine_part]->prefinement;

                        IntVector coarse_lower = pdata->ilowers[coarse_box];
                        IntVector coarse_upper = pdata->iuppers[coarse_box] + IntVector ( 1, 1, 1 );
                        IntVector fine_lower, fine_upper;
                        std::vector<const Patch *> fine_patches;

                        for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                        {
                            fine_lower[d] = coarse_lower[d] * prefinement[d];
                            fine_upper[d] = coarse_upper[d] * prefinement[d];
                        }

                        fine_level->selectPatches ( fine_lower, fine_upper, fine_patches );

                        for ( const Patch * fine_patch : fine_patches )
                        {
                            IntVector fine_index = fine_patch->getCellLowIndex();
                            IntVector coarse_index;
                            for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                                coarse_index[d] = fine_index[d] / prefinement[d] - ( fine_index[d] < 0 );

                            for ( int var = 0; var < nvars; ++var )
                            {
                                HYPRE(SStructGraphAddEntries) ( solver_struct->graph, coarse_part, coarse_index.get_pointer(), var, fine_part, fine_index.get_pointer(), var );
                                extra_values[var][coarse_part][coarse_index].emplace_back ( 0 );
                            }
                        }

                        for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                        {
                            if ( d < DIM )
                            {
                                coarse_lower[d] -= 1;
                                coarse_upper[d] += 1;
                            }
                            fine_lower[d] = coarse_lower[d] * prefinement[d];
                            fine_upper[d] = coarse_upper[d] * prefinement[d];
                        }

                        fine_level->selectPatches ( fine_lower, fine_upper, fine_patches );

                        if ( fine_patches.empty() )
                            continue;

                        CCVariable<Stencil7> stencil_var;
                        matrix_dw->getModifiable ( stencil_var, m_stencil_entries_label, material, coarse_patch );

                        for ( const Patch * fine_patch : fine_patches )
                        {
                            IntVector fine_lower = fine_patch->getCellLowIndex();
                            IntVector fine_upper = fine_patch->getCellHighIndex();

                            for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                            {
                                fine_upper[d] -= 1;
                                coarse_lower[d] = fine_lower[d] / prefinement[d] - ( fine_lower[d] < 0 );
                                coarse_upper[d] = fine_upper[d] / prefinement[d] - ( fine_upper[d] < 0 );
                            }

                            IntVector fine_index, coarse_index, ilower, iupper;
                            for ( size_t d = 0; d < DIM; ++d )
                                for ( int s :
                                        {
                                            -1, 1
                                        } )
                                {
                                    size_t f = 2 * d + ( s + 1 ) / 2;
                                    if ( fine_patch->getBCType ( ( Patch::FaceType ) f ) == Patch::Coarse )
                                    {
                                        ilower = coarse_lower;
                                        iupper = coarse_upper;

                                        if ( fine_patch->isVirtual() )
                                        {
                                            fine_lower = fine_patch->getRealPatch()->getCellLowIndex();
                                            fine_upper = fine_patch->getRealPatch()->getCellHighIndex();
                                        }
                                        else
                                        {
                                            fine_lower = fine_patch->getCellLowIndex();
                                            fine_upper = fine_patch->getCellHighIndex();
                                        }
                                        if ( s < 0 )
                                        {
                                            iupper[d] = ilower[d] += s;
                                        }
                                        else
                                        {
                                            ilower[d] = iupper[d] += s;
                                            fine_lower[d] = fine_upper[d] - prefinement[d] + 1;
                                        }

                                        ilower = Max ( ilower, pdata->ilowers[coarse_box] );
                                        iupper = Min ( iupper, pdata->iuppers[coarse_box] );

                                        for ( coarse_index[0] = ilower[0], fine_index[0] = fine_lower[0]; coarse_index[0] <= iupper[0]; ++coarse_index[0], fine_index[0] += prefinement[0] )
                                            for ( coarse_index[1] = ilower[1], fine_index[1] = fine_lower[1]; coarse_index[1] <= iupper[1]; ++coarse_index[1], fine_index[1] += prefinement[1] )
                                                for ( coarse_index[2] = ilower[2], fine_index[2] = fine_lower[2]; coarse_index[2] <= iupper[2]; ++coarse_index[2], fine_index[2] += prefinement[2] )
                                                    for ( int var = 0; var < nvars; ++var )
                                                        refine_interface ( coarse_index, fine_index, var, d, s, stencil_var[coarse_index][f - s] );


                                    }
                                }

                        }
                    }
                }
#endif

#if 1
                /* add fine to coarse entries */
                for ( int fine_part = 1; fine_part < nparts; ++fine_part )
                {
                    auto & pdata = solver_struct->pdatas[fine_part];

                    std::function<void ( const IntVector &, const IntVector &, const double &, CCVariable<Stencil7> &, AdditionalEntries & ) > refine_interface;

                    switch ( m_params->C2F )
                    {
                    case 0:  // interpolate coarse node with fine node
                    {
                        refine_interface = [this, fine_part] (
                                               const IntVector & index,
                                               const IntVector & to_index,
                                               const double & value,
                                               CCVariable<Stencil7> & stencil_entries,
                                               AdditionalEntries & additional_entries
                                           )
                        {
                            if ( value )
                            {
                                add_fine_to_fine_entry ( fine_part, index, to_index, value, stencil_entries, additional_entries );
                            }
                        };
                    }
                    break;
                    case 1:
                    {
                        double n_fine = pdata->prefinement[0] * pdata->prefinement[1] * pdata->prefinement[2];
                        refine_interface = [this, fine_part, pdata, n_fine, &extra_values] (
                                               const IntVector & index,
                                               const IntVector & to_index,
                                               const double & value,
                                               CCVariable<Stencil7> & stencil_entries,
                                               AdditionalEntries & additional_entries
                                           )
                        {
                            if ( value )
                            {
                                double entry = value / n_fine;
                                for ( int k = 0; k < pdata->prefinement[2]; ++k )
                                    for ( int j = 0; j < pdata->prefinement[1]; ++j )
                                        for ( int i = 0; i < pdata->prefinement[0]; ++i )
                                        {
                                            add_fine_to_fine_entry ( fine_part, index, {to_index[0] + i, to_index[1] + j, to_index[2] + k}, entry, stencil_entries, additional_entries );
                                        }
                            }
                        };
                    }
                    break;
                    default:
                        ASSERTFAIL ( "C2F value not handeld" );
                    }

                    for ( int fine_box = 0; fine_box < nboxes[rank][fine_part]; ++fine_box )
                    {
                        CCVariable<Stencil7> stencil_var;
                        const Patch * patch = grd->getPatchByID ( pdata->box2patch[fine_box], fine_part );
                        AdditionalEntries f2f_entries;

                        auto iter = additional_var[fine_part][fine_box]->begin();
                        while ( iter != additional_var[fine_part][fine_box]->end() )
                        {
                            // look if entries the coarse_index is refined
                            const IntVector & fine_index = std::get<0> ( iter->first );
                            const int & coarse_part = std::get<1> ( iter->first );
                            const IntVector & coarse_index = std::get<2> ( iter->first );
                            const double & value = iter->second;

                            bool is_refined = false;
                            if ( fine_part - coarse_part == 1 )
                            {
                                IntVector to_fine_index;
                                for ( size_t i = 0; i < HYPRE_MAXDIM; ++i ) // till last dimension otherwise getPatchFromIndex fails!
                                    to_fine_index[i] = coarse_index[i] * pdata->prefinement[i];

                                for ( int to_fine_box = 0; to_fine_box < nboxes[rank][fine_part]; ++to_fine_box )
                                    if ( pdata->ilowers[to_fine_box] <= to_fine_index && to_fine_index <= pdata->iuppers[to_fine_box] )
                                    {
                                        is_refined = true;

                                        matrix_dw->getModifiable ( stencil_var, m_stencil_entries_label, material, patch, Ghost::None, 0 );
                                        refine_interface ( fine_index, to_fine_index, value, stencil_var, f2f_entries );

                                        break;
                                    }
                            }
                            if ( is_refined )
                                iter = additional_var[fine_part][fine_box]->erase ( iter );
                            else
                                ++iter;
                        }

                        auto hint = additional_var[fine_part][fine_box]->end();
                        for ( const auto & entry : f2f_entries )
                            hint = additional_var[fine_part][fine_box]->emplace_hint ( hint, entry );

                        for ( const auto & entry : *additional_var[fine_part][fine_box] )
                        {
                            IntVector fine_index = std::get<0> ( entry.first );
                            const int & coarse_part = std::get<1> ( entry.first );
                            IntVector coarse_index = std::get<2> ( entry.first );
                            const double & value = entry.second;

                            for ( int var = 0; var < nvars; ++var )
                            {
                                HYPRE(SStructGraphAddEntries) ( solver_struct->graph, fine_part, fine_index.get_pointer(), var, coarse_part, coarse_index.get_pointer(), var );
                                extra_values[var][fine_part][fine_index].emplace_back ( value );
                            }
                        }
                    }
                }
#endif

                for ( int fine_part = 0; fine_part < nparts; ++fine_part )
                    delete[] additional_var[fine_part];
                delete[] additional_var;
                }
                HYPRE(SStructGraphAssemble) ( solver_struct->graph );

                /*-----------------------------------------------------------
                 * Set up the matrix
                 *-----------------------------------------------------------*/

                DOUT ( dbg_assembling, rank << "   hypre matrix setup" );

                HYPRE(SStructMatrixCreate) ( comm, solver_struct->graph, &solver_struct->A );
                HYPRE(SStructMatrixSetObjectType) ( solver_struct->A, HYPRE_SSTRUCT );
                HYPRE(SStructMatrixInitialize) ( solver_struct->A );
            }

            if ( do_update || restart )
            {
                /* set stencil values */
                DOUT ( dbg_assembling, rank << "   hypre update matrix stencil entries" );
                for ( int p = 0; p < patches->size(); ++p )
                    for ( int var = 0; var < nvars; ++var )
                    {
                        const Patch * patch = patches->get ( p );

                        int part = patch->getLevel()->getIndex();
                        auto & pdata = solver_struct->pdatas[part];
                        int box = pdata->patch2box[patch->getID()];

                        constCCVariable<Stencil7> stencil_var;
                        matrix_dw->get ( stencil_var, m_stencil_entries_label, material, patch, Ghost::None, 0 );

                        std::vector<double> values = get_stencil_values<DIM> ( pdata, box, stencil_var );

                        std::vector<int> stencil_indices ( stencil_size );
                        std::iota ( std::begin ( stencil_indices ), std::end ( stencil_indices ), 0 );
                        HYPRE(SStructMatrixSetBoxValues) ( solver_struct->A, part, pdata->ilowers[box].get_pointer(), pdata->iuppers[box].get_pointer(), var, stencil_size, stencil_indices.data(), values.data() );
                    }

                if ( nparts > 1 )
                {
                    for ( int var = 0; var < nvars; ++var )
                        for ( int part = 0; part < nparts; ++part )
                            for ( std::pair<IntVector, std::vector<double>> && entry : extra_values[var][part] )
                            {
                                auto & index = entry.first;
                                auto & values = entry.second;
                                auto nvalues = values.size();
                                std::vector<int> entries ( nvalues );
                                std::iota ( entries.begin(), entries.end(), stencil_size );
                                HYPRE(SStructMatrixSetValues) ( solver_struct->A, part, index.get_pointer(), var, nvalues, entries.data(), values.data() );
                            }

                    /* reset matrix values so that stencil connections between two parts are zeroed */
                    DOUT ( dbg_assembling, rank << "   hypre reset matix dangling entries" );
                    for ( int part = nparts - 1; part > 0; part-- )
                    {
                        auto & pdata = solver_struct->pdatas[part];
                        HYPRE(SStructFACZeroCFSten) ( solver_struct->A, solver_struct->grid, part, pdata->prefinement );
                        HYPRE(SStructFACZeroFCSten) ( solver_struct->A, solver_struct->grid, part );
                        HYPRE(SStructFACZeroAMRMatrixData) ( solver_struct->A, part - 1, pdata->prefinement );
                    }
                }
            }

            if ( do_setup || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre matrix assemble" );
                HYPRE(SStructMatrixAssemble) ( solver_struct->A );

                /*-----------------------------------------------------------
                 * Set up the linear system
                 *-----------------------------------------------------------*/
                DOUT ( dbg_assembling, rank << "   hypre rhs vector setup" );
                HYPRE(SStructVectorCreate) ( comm, solver_struct->grid, &solver_struct->b );
//              HYPRE(SStructVectorSetObjectType) ( solver_struct->b, object_type );
                HYPRE(SStructVectorInitialize) ( solver_struct->b );
            }

            DOUT ( dbg_assembling, rank << "   hypre rhs vector update coefficients" );
            for ( int p = 0; p < patches->size(); ++p )
                for ( int var = 0; var < nvars; ++var )
                {
                    const Patch * patch = patches->get ( p );
                    int part = patch->getLevel()->getIndex();

                    auto & pdata = solver_struct->pdatas[part];
                    int box = pdata->patch2box[patch->getID()];

                    constCCVariable<double> rhs_var;
                    rhs_dw->get ( rhs_var, m_rhs_label, material, patch, Ghost::None, 0 );

                    for ( int k = pdata->ilowers[box][2]; k <= pdata->iuppers[box][2]; ++k )
                        for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
                        {
                            const double * vals = &rhs_var[IntVector ( pdata->ilowers[box][0], j, k )];
                            IntVector ll ( pdata->ilowers[box][0], j, k );
                            IntVector hh ( pdata->iuppers[box][0], j, k );
                            HYPRE(SStructVectorSetBoxValues) ( solver_struct->b, part, ll.get_pointer(), hh.get_pointer(), var, const_cast<double *> ( vals ) );
                        }
                }

            if ( nparts > 1 )
                HYPRE(SStructFACZeroAMRVectorData) ( solver_struct->b, plevels, prefinements );

            if ( do_setup || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre rhs vector assemble" );
                HYPRE(SStructVectorAssemble) ( solver_struct->b );

                DOUT ( dbg_assembling, rank << "   hypre solution vector setup" );
                HYPRE(SStructVectorCreate) ( comm, solver_struct->grid, &solver_struct->x );
//              SStructVectorSetObjectType ( solver_struct->x, object_type );
                HYPRE(SStructVectorInitialize) ( solver_struct->x );
            }

            if ( m_guess_label )
            {
                DOUT ( dbg_assembling, rank << "   hypre solution update coefficients" );
                for ( int p = 0; p < patches->size(); ++p )
                    for ( int var = 0; var < nvars; ++var )
                    {
                        const Patch * patch = patches->get ( p );

                        int part = patch->getLevel()->getIndex();
                        auto & pdata = solver_struct->pdatas[part];
                        int box = pdata->patch2box[patch->getID()];

                        constCCVariable<double> guess_var;
                        guess_dw->get ( guess_var, m_guess_label, material, patch, Ghost::None, 0 );

                        for ( int k = pdata->ilowers[box][2]; k <= pdata->iuppers[box][2]; k++ )
                            for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
                            {
                                const double * vals = &guess_var[IntVector ( pdata->ilowers[box][0], j, k )];
                                IntVector ll ( pdata->ilowers[box][0], j, k );
                                IntVector hh ( pdata->iuppers[box][0], j, k );
                                HYPRE(SStructVectorSetBoxValues) ( solver_struct->x, part, ll.get_pointer(), hh.get_pointer(), var, const_cast<double *> ( vals ) );
                            }
                    }

                if ( nparts > 1 )
                    HYPRE(SStructFACZeroAMRVectorData) ( solver_struct->x, plevels, prefinements );
            }

            if ( do_setup || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre solution vector assemble" );
                HYPRE(SStructVectorAssemble) ( solver_struct->x );
                solver_struct->created = true;
            }

#ifdef HYPRE_TIMING
            hypre_EndTiming ( tMatVecSetup_ );
#endif

            Timers::Simple solve_timer;
            solve_timer.start();

#ifdef HYPRE_TIMING
            hypre_BeginTiming ( tSolveOnly_ );
#endif

            /*-----------------------------------------------------------
            * Solve the system using FAC
            *-----------------------------------------------------------*/
            DOUT ( dbg_assembling, rank << "   hypre solve" );

            int   num_iterations = -1;
            double  final_res_norm = NAN;

            if ( do_setup || restart )
            {
                if ( nparts > 1 )
                    solver_struct->solver = scinew SolverImpl<FAC> ( comm, solver_struct, m_params, plevels, prefinements, m_guess_label );
                else if ( m_params->csolver_type == SysPFMG_PCG )
                    solver_struct->solver = scinew SolverImpl<PFMG_PCG> ( comm, solver_struct, m_params );
                else if ( m_params->csolver_type == SysPFMG )
                    solver_struct->solver = scinew SolverImpl<PFMG> ( comm, solver_struct, m_params, m_guess_label );
            }

            solver_struct->solver->solve ( solver_struct, num_iterations, final_res_norm );

            /*-----------------------------------------------------------
            * Gather the solution vector
            *-----------------------------------------------------------*/
            HYPRE(SStructVectorGather) ( solver_struct->x );

            /*-----------------------------------------------------------
            * Print the solution and other info
            *-----------------------------------------------------------*/

#ifdef PRINTSYSTEM
            //__________________________________
            //   Debugging
            std::ostringstream sname;
            std::vector<std::string> fname;
            sname << "." << d_myworld->myRank() << "." << timeStep;
            m_params->setOutputFileName ( sname.str() );
            m_params->getOutputFileName ( fname );

            for ( auto name : fname )
                DOUT ( dbg_doing, "   hypre print " << name );
            HYPRE_SStructMatrixPrint ( fname[0].c_str(), solver_struct->A, 0 );
            HYPRE_SStructVectorPrint ( fname[1].c_str(), solver_struct->b, 0 );
            HYPRE_SStructVectorPrint ( fname[2].c_str(), solver_struct->x, 0 );
#endif

            //__________________________________
            // Test for convergence
            if ( final_res_norm > m_params->tol || std::isfinite ( final_res_norm ) == 0 )
            {
                if ( m_params->getRecomputeTimeStepOnFailure() )
                {
                    if ( pg->myRank() == 0 )
                        std::cout << "HypreSolver not converged in " << num_iterations
                                  << "iterations, final residual= " << final_res_norm
                                  << ", requesting smaller timestep\n";
                    //new_dw->abortTimestep();
                    //new_dw->restartTimestep();
                }
                else
                {
                    throw ConvergenceFailure ( "HypreSolver variable: " + m_solution_label->getName() + ", solver: " + ( ( nparts > 1 ) ? "FAC" : ( m_params->csolver_type == SysPFMG_PCG ) ? "SysPFMG_PCG" : ( m_params->csolver_type == SysPFMG ) ? "SysPFMG" : "???" ),
                                               num_iterations, final_res_norm,
                                               m_params->tol, __FILE__, __LINE__ );
                }
            }

            solve_timer.stop();

#ifdef HYPRE_TIMING
            hypre_EndTiming ( tSolveOnly_ );
#endif

            //__________________________________
            // Push the solution into Uintah data structure
            for ( int p = 0; p < patches->size(); p++ )
                for ( int var = 0; var < nvars; ++var )
                {
                    const Patch * patch = patches->get ( p );

                    int part = patch->getLevel()->getIndex();
                    auto & pdata = solver_struct->pdatas[part];
                    int box = pdata->patch2box[patch->getID()];

                    CCVariable<double> solution_var;
                    if ( m_modifies_solution )
                        new_dw->getModifiable ( solution_var, m_solution_label, material, patch );
                    else
                        new_dw->allocateAndPut ( solution_var, m_solution_label, material, patch );

                    // Get the solution back from hypre
                    for ( int k = pdata->ilowers[box][2]; k <= pdata->iuppers[box][2]; k++ )
                        for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
                        {
                            double * vals = &solution_var[IntVector ( pdata->ilowers[box][0], j, k )];
                            IntVector ll ( pdata->ilowers[box][0], j, k );
                            IntVector hh ( pdata->iuppers[box][0], j, k );
                            HYPRE(SStructVectorGetBoxValues) ( solver_struct->x, part, ll.get_pointer(), hh.get_pointer(), var, vals );
                        }
                }

            //__________________________________
            // Push the solution into Uintah data structure
            for ( int p = 0; p < patches->size(); p++ )
                for ( int var = 0; var < nvars; ++var )
                {
                    const Patch * patch = patches->get ( p );

                    int part = patch->getLevel()->getIndex();
                    auto & pdata = solver_struct->pdatas[part];
                    int box = pdata->patch2box[patch->getID()];

                    CCVariable<double> rhs_var;
                    CCVariable<Stencil7> stencil_var;
                    rhs_dw->getModifiable ( rhs_var, m_rhs_label, material, patch );
                    matrix_dw->getModifiable ( stencil_var, m_stencil_entries_label, material, patch );

                    // Get the solution back from hypre
                    for ( int k = pdata->ilowers[box][2]; k <= pdata->iuppers[box][2]; k++ )
                        for ( int j = pdata->ilowers[box][1]; j <= pdata->iuppers[box][1]; ++j )
                            for ( int i = pdata->ilowers[box][0]; i <= pdata->iuppers[box][0]; ++i )
                            {
                                int index[3] = {i, j, k};
                                double * rhs = &rhs_var ( i, j, k );
                                double * Ap = &stencil_var ( i, j, k ).p;
                                double * Ae = &stencil_var ( i, j, k ).e;
                                double * Aw = &stencil_var ( i, j, k ).w;
                                double * An = &stencil_var ( i, j, k ).n;
                                double * As = &stencil_var ( i, j, k ).s;

                                HYPRE(SStructVectorGetValues) ( solver_struct->b, part, index, var, rhs );
                                HYPRE(SStructMatrixGetValues) ( solver_struct->A, part, index, var, 1, &stn0, Ap );
                                HYPRE(SStructMatrixGetValues) ( solver_struct->A, part, index, var, 1, &stnE, Ae );
                                HYPRE(SStructMatrixGetValues) ( solver_struct->A, part, index, var, 1, &stnW, Aw );
                                HYPRE(SStructMatrixGetValues) ( solver_struct->A, part, index, var, 1, &stnN, An );
                                HYPRE(SStructMatrixGetValues) ( solver_struct->A, part, index, var, 1, &stnS, As );

                            }
                }

            if ( nparts > 1 )
            {
                delete[] plevels;
                delete[] prefinements;
            }

#ifdef HYPRE_TIMING
            hypre_EndTiming ( tHypreAll_ );

            hypre_PrintTiming ( "Hypre Timings:", pg->getComm() );
            hypre_FinalizeTiming ( tMatVecSetup_ );
            hypre_FinalizeTiming ( tSolveOnly_ );
            hypre_FinalizeTiming ( tHypreAll_ );
            hypre_ClearTiming();
#endif

            timer.stop();

            if ( pg->myRank() == 0 )
            {
                std::cout << "Solve of " << m_solution_label->getName()
                          << " completed in " << timer().seconds()
                          << " s ( solve only: " << solve_timer().seconds() << " s, ";

                if ( timeStep > 2 )
                {
                    // alpha = 2/(N+1)
                    // averaging window is 10 timeSteps.
                    double alpha = 2.0 / ( std::min ( int ( timeStep ) - 2, 10 ) + 1 );
                    m_moving_average = alpha * solve_timer().seconds() + ( 1 - alpha ) * m_moving_average;

                    std::cout << "mean: " <<  m_moving_average << " s, ";
                }

                std::cout << num_iterations << " iterations, residual = "
                          << final_res_norm << " )." << std::endl;
            }

            timer.reset ( true );
        }
    }

    const VarLabel * hypre_global_data_label ()
    {
        return m_hypre_global_data_label;
    }

    const VarLabel * hypre_solver_struct_label ( const int mpi_rank )
    {
        auto res = m_hypre_solver_struct_label.find ( mpi_rank );
        if ( res != m_hypre_solver_struct_label.end() )
            return res->second;
        return m_hypre_solver_struct_label.insert ( std::make_pair ( mpi_rank, VarLabel::create ( "hypre_solver_struct_" + std::to_string ( mpi_rank ), SoleVariable<SolverStructP>::getTypeDescription() ) ) ).first->second;
    }
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreFAC_Solver_h
