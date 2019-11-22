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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Solver_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Solver_h

#include <CCA/Components/Solvers/HypreSStruct/Definitions.h>

#include <CCA/Components/Solvers/HypreSStruct/AdditionalEntriesP.h>
#include <CCA/Components/Solvers/HypreSStruct/AdditionalEntries.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructInterfaceP.h>
#include <CCA/Components/Solvers/HypreSStruct/GlobalData.h>
#include <CCA/Components/Solvers/HypreSStruct/PartData.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverFactory.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructInterface.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>
#include <CCA/Components/Solvers/SolverCommon.h>

#include <Core/Grid/Variables/PerPatch.h> // must be included after GlobalDataP/SolverStructP/AdditionalEntriesP where swapbytes override is defined
#include <Core/Grid/DbgOutput.h>
#include <Core/Util/Timers/Timers.hpp>
#include <Core/Exceptions/ConvergenceFailure.h>

/**
 *  @class  HypreSStruct::Solver
 *  @author Jon Matteo Church
 *  @brief  Uintah hypre solver interface.
 *  Allows the solution of a linear system of the form \[ \mathbf{A} \mathbf{x} = \mathbf{b}\] where \[\mathbf{A}\] is
 *  stencil7 matrix.
 */

namespace Uintah
{
namespace HypreSStruct
{

static constexpr bool dbg_doing = false;
static constexpr bool dbg_assembling = false;
extern DebugStream cout_doing;

template < int DIM >
class Solver
    : public SolverCommon
{
public:
    static std::string AdditionalEntriesSuffix;

private:
    static constexpr Patch::FaceType face_start = ( Patch::FaceType ) ( 0 );
    static constexpr Patch::FaceType face_end = ( Patch::FaceType ) ( 2 * DIM - 1 );

    const VarLabel * m_timeStepLabel;
    const VarLabel * m_hypre_global_data_label;
    std::map <int, const VarLabel *> m_hypre_sstruct_interface_label;

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

    const Ghost::GhostType gtype;
    const int gnum;

    SStructInterfaceFactory::FactoryMethod create_sstruct_interface;

public:
    Solver (
        const ProcessorGroup * myworld,
        const int & ghosts,
        SStructInterfaceFactory::FactoryMethod && sstruct_interface_creator
    )
        : SolverCommon ( myworld )
        , m_timeStepLabel ( VarLabel::create ( timeStep_name, timeStep_vartype::getTypeDescription() ) )
        , m_hypre_global_data_label ( VarLabel::create ( "hypre_global_data", SoleVariable<GlobalDataP>::getTypeDescription() ) )
        , m_hypre_sstruct_interface_label ()
        , m_params ( scinew SolverParams() )
        , m_moving_average ( 0. )
        , gtype ( ghosts ? Ghost::AroundCells : Ghost::None )
        , gnum ( ghosts )
        , create_sstruct_interface ( sstruct_interface_creator )
    {
    }

    virtual ~Solver()
    {
        VarLabel::destroy ( m_timeStepLabel );
        VarLabel::destroy ( m_hypre_global_data_label );
        for ( auto && l : m_hypre_sstruct_interface_label ) VarLabel::destroy ( l.second );
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
                    continue;

                int sFreq;
                int coefFreq;
                param_ps->getWithDefault ( "solveFrequency",      m_params->solveFrequency, 1 );
                param_ps->getWithDefault ( "setupFrequency",      sFreq,                    1 );
                param_ps->getWithDefault ( "updateCoefFrequency", coefFreq,                 1 );
                m_params->setSetupFrequency ( sFreq );
                m_params->setUpdateCoefFrequency ( coefFreq );

                int cSolverType;
                int relaxType;
                ASSERTMSG( param_ps->get ( "max_levels",    m_params->max_levels ),   "HypreSStruct::Solver ERROR. Missing parameter: max_level" );
                param_ps->getWithDefault ( "maxiterations", m_params->max_iter,       -1  );
                param_ps->getWithDefault ( "tolerance",     m_params->tol,            -1. );
                param_ps->getWithDefault ( "rel_change",    m_params->rel_change,     -1  );
                param_ps->getWithDefault ( "relax_type",    relaxType,                -1  );
                param_ps->getWithDefault ( "weight",        m_params->weight,         -1. );
                param_ps->getWithDefault ( "npre",          m_params->num_pre_relax,  -1  );
                param_ps->getWithDefault ( "npost",         m_params->num_post_relax, -1  );
                param_ps->getWithDefault ( "csolver_type",  cSolverType,              -1  );
                param_ps->getWithDefault ( "logging",       m_params->logging,        -1  );
                m_params->relax_type = ( RelaxType ) relaxType;
                m_params->csolver_type = ( CoarseSolverType ) cSolverType;

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
        }
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

        printSchedule ( level, cout_doing, "HypreSStruct::Solver:scheduleSolve" );

        auto global_data_label = hypre_global_data_label();
        auto sstruct_interface_label = hypre_sstruct_interface_label ( d_myworld->myRank() );
        auto grid = level->getGrid();

        Task * task = scinew Task ( "HypreSStruct::Solver::solve", this, &Solver::solve, global_data_label, sstruct_interface_label );
        task->requires ( Task::OldDW, sstruct_interface_label );
        task->requires ( Task::OldDW, global_data_label );

        task->computes ( global_data_label );
        task->computes ( sstruct_interface_label );

        task->requires ( matrix_dw, m_stencil_entries_label, gtype, gnum );
        task->requires ( matrix_dw, m_additional_entries_label, (MaterialSubset*)nullptr );

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
        auto sstruct_interface_label = hypre_sstruct_interface_label ( d_myworld->myRank() );

        printSchedule ( level, cout_doing, "HypreSStruct::Solver:scheduleInitialize" );

        Task * task = scinew Task ( "initialize_hypre", this, &Solver::initialize, level );

        task->computes ( global_data_label );
        task->computes ( sstruct_interface_label );

        task->setType ( Task::OncePerProc );

        ASSERT ( scheduler );
        ASSERT ( scheduler->getLoadBalancer() );
        auto * patches = scheduler->getLoadBalancer()->getPerProcessorPatchSet ( level );
        scheduler->addTask ( task, patches, matls );

        scheduler->overrideVariableBehavior ( global_data_label->getName(), false, false, false, true, true );
        scheduler->overrideVariableBehavior ( sstruct_interface_label->getName(), false, false, false, true, true );
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
        return "hypre_sstruct";
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
        const int rank = pg->myRank();
        const int nparts = level->getGrid()->numLevels();
        const int part = level->getIndex();
        const int nboxes = level->numPatches();
        const int npatches = patches->size();

        DOUT ( dbg_doing, "\nHypreSStruct::Solver::initialize ( rank " << rank << " level " << part << " ) " );

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

        SStructInterface * sstruct_interface;
        auto sstruct_interface_label = hypre_sstruct_interface_label ( rank );
        if ( new_dw->exists ( sstruct_interface_label ) )
        {
            SoleVariable<SStructInterfaceP> sstruct_interface_var;
            new_dw->get ( sstruct_interface_var, sstruct_interface_label );
            sstruct_interface = sstruct_interface_var.get().get_rep();
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
                global_data = scinew GlobalData ( nparts );
                SoleVariable<GlobalDataP> global_data_var;
                global_data_var.setData ( global_data );
                new_dw->put ( global_data_var, global_data_label );
            }

            sstruct_interface = dynamic_cast<SStructInterface *> ( create_sstruct_interface ( global_data ) );
            SoleVariable<SStructInterfaceP> sstruct_interface_var;
            sstruct_interface_var.setData ( sstruct_interface );
            new_dw->put ( sstruct_interface_var, sstruct_interface_label );
        }

        if ( !sstruct_interface->isPartSet ( part ) )
        {
            sstruct_interface->setPart ( part, npatches, nboxes,
                                         low.get_pointer(), high.get_pointer(),
                                         refinement.get_pointer(),
                                         periodic.get_pointer() );

            for ( int p = 0; p < patches->size(); p++ )
            {
                const Patch * patch = patches->get ( p );

                int id = patch->getID();
                IntVector ilower = patch->getCellLowIndex();
                IntVector iupper = patch->getCellHighIndex();

                bool interface[2 * HYPRE_MAXDIM] = { false, false, false, false, false, false };
                for ( Patch::FaceType f = face_start; f <= face_end; f = Patch::nextFace ( f ) )
                    if ( patch->getBCType ( f ) == Patch::Coarse )
                        interface[f] = true;

                sstruct_interface->addBox ( part, id, ilower.get_pointer(),
                                            iupper.get_pointer(), interface );
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
        const VarLabel * sstruct_interface_label
    )
    {
        auto comm = pg->getComm();
        auto rank = pg->myRank();

        DOUT ( dbg_doing, "\nHypreSStruct::Solver::solve (rank " << rank << ")" );

        struct SStructInterface * sstruct_interface = nullptr;

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

        ASSERTMSG ( old_dw->exists ( sstruct_interface_label ), "old hypre fac solver struct should exist" );
        ASSERTMSG ( !new_dw->exists ( sstruct_interface_label ), "new hypre fac solver struct should not exist" );
        SoleVariable<SStructInterfaceP> sstruct_interface_var;
        old_dw->get ( sstruct_interface_var, sstruct_interface_label );
        new_dw->put ( sstruct_interface_var, sstruct_interface_label );
        sstruct_interface = sstruct_interface_var.get().get_rep();

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

        bool restart = sstruct_interface->restart();

        for ( int m = 0; m < matls->size(); ++m )
        {
            int material = matls->get ( m );

#ifdef HYPRE_TIMING
            hypre_BeginTiming ( tMatVecSetup_ );
#endif

            GridP grd = nullptr;
            constCCVariable<Stencil7> ** stencil_entries = nullptr;
            AdditionalEntries ** * additional_entries = nullptr;
            constCCVariable<double> ** rhs = nullptr;
            constCCVariable<double> ** guess = nullptr;
            CCVariable<double> ** solution = nullptr;

            const GlobalDataP & gdata = sstruct_interface->globalData();
            const int & nparts = gdata->nParts();

            if ( patches->size() )
            {
                grd = patches->get ( 0 )->getLevel()->getGrid();

                stencil_entries = scinew constCCVariable<Stencil7> * [nparts];
                additional_entries = scinew AdditionalEntries ** [nparts];
                rhs = scinew constCCVariable<double> * [nparts];
                solution = scinew CCVariable<double> * [nparts];
                for ( int part = 0; part < nparts; ++part )
                {
                    const PartDataP & pdata = sstruct_interface->partData ( part );
                    const int & nboxes = gdata->nBoxes ( part );
                    stencil_entries[part] = scinew constCCVariable<Stencil7> [nboxes];
                    additional_entries[part] = scinew AdditionalEntries * [nboxes];
                    rhs[part] = scinew constCCVariable<double> [nboxes];
                    for ( int box = 0; box < nboxes; ++box )
                    {
                        const Patch * patch = grd->getPatchByID ( pdata->patch ( box ), part );
                        matrix_dw->get ( stencil_entries[part][box], m_stencil_entries_label, material, patch, gtype, gnum );
                        rhs_dw->get ( rhs[part][box], m_rhs_label, material, patch, Ghost::None, 0 );

                        PerPatch<AdditionalEntriesP> patch_additional_entries;
                        matrix_dw->get ( patch_additional_entries, m_additional_entries_label, material, patch );
                        additional_entries[part][box] = patch_additional_entries.get().get_rep();
                    }
                }
            }

            if ( m_guess_label )
            {
                guess = scinew constCCVariable<double> * [nparts];
                for ( int part = 0; part < nparts; ++part )
                {
                    const PartDataP & pdata = sstruct_interface->partData ( part );
                    const int & nboxes = gdata->nBoxes ( part );
                    guess[part] = scinew constCCVariable<double> [nboxes];
                    for ( int box = 0; box < nboxes; ++box )
                    {
                        const Patch * patch = grd->getPatchByID ( pdata->patch ( box ), part );
                        guess_dw->get ( guess[part][box], m_guess_label, material, patch, Ghost::None, 0 );
                    }
                }
            }

            if ( m_modifies_solution )
            {
                for ( int part = 0; part < nparts; ++part )
                {
                    const PartDataP & pdata = sstruct_interface->partData ( part );
                    const int & nboxes = gdata->nBoxes ( part );
                    solution[part] = scinew CCVariable<double> [nboxes];
                    for ( int box = 0; box < nboxes; ++box )
                    {
                        const Patch * patch = grd->getPatchByID ( pdata->patch ( box ), part );
                        new_dw->getModifiable ( solution[part][box], m_solution_label, material, patch );
                    }
                }
            }
            else
            {
                for ( int part = 0; part < nparts; ++part )
                {
                    const PartDataP & pdata = sstruct_interface->partData ( part );
                    const int & nboxes = gdata->nBoxes ( part );
                    solution[part] = scinew CCVariable<double> [nboxes];
                    for ( int box = 0; box < nboxes; ++box )
                    {
                        const Patch * patch = grd->getPatchByID ( pdata->patch ( box ), part );
                        new_dw->allocateAndPut ( solution[part][box], m_solution_label, material, patch );
                    }
                }
            }

            if ( do_setup || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre grid setup" );
                sstruct_interface->gridInitialize ( comm );

                DOUT ( dbg_assembling, rank << "   hypre stencil setup" );
                sstruct_interface->stencilInitialize ( comm );

                DOUT ( dbg_assembling, rank << "   hypre graph setup" );
                sstruct_interface->graphInitialize ( comm, grd, &additional_entries );

                DOUT ( dbg_assembling, rank << "   hypre matrix setup" );
                sstruct_interface->matrixInitialize ( comm );

                DOUT ( dbg_assembling, rank << "   hypre rhs vector setup" );
                sstruct_interface->rhsInitialize ( comm );

                DOUT ( dbg_assembling, rank << "   hypre solution vector setup" );
                sstruct_interface->solutionInitialize ( comm );

                DOUT ( dbg_assembling, rank << "   hypre solver initialize" );
                sstruct_interface->solverInitialize ( comm, m_params );
            }

            DOUT ( dbg_assembling, rank << "   hypre rhs vector update coefficients" );
            sstruct_interface->rhsUpdate ( &rhs );

            if ( do_update || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre matrix update coefficients" );
                sstruct_interface->matrixUpdate ( &stencil_entries, &additional_entries );
            }

            if ( m_guess_label )
            {
                DOUT ( dbg_assembling, rank << "   hypre solution update coefficients" );
                sstruct_interface->guessUpdate ( &guess );
            }

            if ( do_setup || restart )
            {
                sstruct_interface->assemble();
            }

            if ( do_update || restart )
            {
                DOUT ( dbg_assembling, rank << "   hypre solver update" );
                sstruct_interface->solverUpdate();
            }

#ifdef HYPRE_TIMING
            hypre_EndTiming ( tMatVecSetup_ );
#endif

            Timers::Simple solve_timer;
            solve_timer.start();

#ifdef HYPRE_TIMING
            hypre_BeginTiming ( tSolveOnly_ );
#endif
            //__________________________________
            //   Solve
            SolverOutput out;

            DOUT ( dbg_assembling, rank << "   hypre solve" );
            sstruct_interface->solve ( &out );

#ifdef PRINTSYSTEM
            //__________________________________
            //   Debugging
            std::ostringstream sname;
            std::vector<std::string> fname;
            sname << "." << rank << "." << timeStep;
            m_params->setOutputFileName ( sname.str() );
            m_params->getOutputFileName ( fname );

            DOUT ( dbg_doing, "   hypre print system" );
            sstruct_interface->printSystem ( fname.data() );
#endif

            //__________________________________
            // Test for convergence
            if ( out.res_norm > m_params->tol )
            {
                if ( m_params->getRecomputeTimeStepOnFailure() )
                {
                    if ( pg->myRank() == 0 )
                        std::cout << "HypreSolver not converged in " << out.num_iterations
                                  << "iterations, final residual= " << out.res_norm
                                  << ", requesting smaller timestep\n";
                    //new_dw->abortTimestep();
                    //new_dw->restartTimestep();
                }
                else
                {
                    throw ConvergenceFailure ( "HypreSolver variable: " + m_solution_label->getName() + ", solver: " + ( ( nparts > 1 ) ? "FAC" : ( m_params->csolver_type == SysPFMG_PCG ) ? "SysPFMG_PCG" : ( m_params->csolver_type == SysPFMG ) ? "SysPFMG" : "???" ),
                                               out.num_iterations, out.res_norm,
                                               m_params->tol, __FILE__, __LINE__ );
                }
            }

            solve_timer.stop();

#ifdef HYPRE_TIMING
            hypre_EndTiming ( tSolveOnly_ );
#endif

            //__________________________________
            // Push the solution into Uintah data structure
            sstruct_interface->getSolution ( &solution );

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

                std::cout << out.num_iterations << " iterations, residual = "
                          << out.res_norm << " )." << std::endl;
            }

            timer.reset ( true );
        }
    }

    const VarLabel * hypre_global_data_label ()
    {
        return m_hypre_global_data_label;
    }

    const VarLabel * hypre_sstruct_interface_label ( const int mpi_rank )
    {
        auto res = m_hypre_sstruct_interface_label.find ( mpi_rank );
        if ( res != m_hypre_sstruct_interface_label.end() )
            return res->second;
        return m_hypre_sstruct_interface_label.insert ( std::make_pair ( mpi_rank, VarLabel::create ( "hypre_sstruct_interface_" + std::to_string ( mpi_rank ), SoleVariable<SStructInterfaceP>::getTypeDescription() ) ) ).first->second;
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Solver_h

