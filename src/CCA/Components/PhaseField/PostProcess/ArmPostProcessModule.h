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
 * @file CCA/Components/PhaseField/PostProcess/ArmPostProcessModule.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/06
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_tip_ArmPostProcessModule_h
#define Packages_Uintah_CCA_Components_PhaseField_tip_ArmPostProcessModule_h

#include <CCA/Components/PostProcessUda/Module.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessorFactory.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessor.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWView.h>
#include <CCA/Ports/Scheduler.h>
#include <forward_list>

namespace Uintah
{
namespace PhaseField
{

/// Debugging stream for component schedulings
static constexpr bool dbg_arm_scheduling = false;

static constexpr bool dbg_arm_tip = false;

template<VarType VAR, DimType DIM, StnType STN, bool AMR> class PureMetal;

template<VarType VAR, DimType DIM, StnType STN, bool AMR = false>
class ArmPostProcessModule
    : public Module

{
    using PM = PureMetal<VAR, DIM, STN, AMR>;

    const Application< PureMetalProblem<VAR, STN>, AMR > * const d_app;
    const Regridder * const d_regridder;

    ProblemSpecP d_param;
    const VarLabel * const d_psi_label;

    ArmPostProcessor<VAR, DIM> * arm_postproc;
    const VarLabel * tip_position_label;
    const VarLabel * tip_velocity_label;
    const VarLabel * tip_curvatures_label;

    double delt;

public: // CONSTRUCTORS/DESTRUCTOR

    ArmPostProcessModule (
        const Application< PureMetalProblem<VAR, STN>, AMR > * app,
        const Regridder * regridder,
        const ProblemSpecP & params,
        const VarLabel * psi_label
    );

    /**
     * @brief Destructor
     */
    virtual ~ArmPostProcessModule();

    /// Prevent copy (and move) constructor
    ArmPostProcessModule ( const ArmPostProcessModule & ) = delete;

    /// Prevent copy (and move) assignment
    /// @return deleted
    ArmPostProcessModule & operator= ( const ArmPostProcessModule & ) = delete;

    virtual std::string
    getName()
    override
    {
        return "pure_metal_arm";
    };

public: // SETUP

    virtual void
    problemSetup ()
    override;

public: // SCHEDULINGS

    /**
     * @brief Schedule the initialization tasks
     *
     * Specify all tasks to be performed at initial timestep to initialize
     * variables in the DataWarehouse
     *
     * @params level grid level to be initialized
     * @params sched scheduler to manage the tasks
     */
    virtual void
    scheduleInitialize (
        SchedulerP & sched,
        const LevelP & level
    ) override;

    virtual void
    scheduleDoAnalysis (
        SchedulerP & sched,
        const LevelP & level
    ) override;

    virtual void
    scheduleDoAnalysis_preReloc (
        SchedulerP & sched,
        const LevelP & level
    ) {};

protected:

    template < bool MG >
    typename std::enable_if < !MG, void >::type
    scheduleDoAnalysis_tip (
        SchedulerP & sched,
        const LevelP & level
    );

    template < bool MG >
    typename std::enable_if < MG, void >::type
    scheduleDoAnalysis_tip (
        SchedulerP & sched,
        const LevelP & level
    );

protected: // TASKS

    void
    task_do_analysis_tip (
        const ProcessorGroup * myworld,
        const PatchSubset * patches,
        const MaterialSubset *,
        DataWarehouse * dw_old,
        DataWarehouse * dw_new
    );

}; // class ArmPostProcessModule

// CONSTRUCTORS/DESTRUCTOR

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
ArmPostProcessModule<VAR, DIM, STN, AMR>::ArmPostProcessModule (
    const Application< PureMetalProblem<VAR, STN>, AMR > * app,
    const Regridder * regridder,
    const ProblemSpecP & params,
    const VarLabel * psi_label
) : d_app ( app ),
    d_regridder ( regridder ),
    d_param ( params ),
    d_psi_label ( psi_label ),
    arm_postproc ( nullptr )
{
    tip_position_label = VarLabel::create ( "tip_position", max_vartype::getTypeDescription() );
    tip_velocity_label = VarLabel::create ( "tip_velocity", max_vartype::getTypeDescription() );
    tip_curvatures_label = VarLabel::create ( "tip_curvatures", maxvec_vartype::getTypeDescription() );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
ArmPostProcessModule<VAR, DIM, STN, AMR>::~ArmPostProcessModule()
{
    delete arm_postproc;

    VarLabel::destroy ( tip_position_label );
    VarLabel::destroy ( tip_velocity_label );
    VarLabel::destroy ( tip_curvatures_label );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
ArmPostProcessModule<VAR, DIM, STN, AMR>::problemSetup ()
{
    ProblemSpecP pure_metal = d_param->findBlock ( "PhaseField" );

    double epsilon;
    pure_metal->require ( "epsilon", epsilon );
    pure_metal->require ( "delt", delt );

    arm_postproc = ArmPostProcessorFactory<VAR, DIM>::create ( d_param->findBlock ( "ArmPostProcessor" ), epsilon, d_app->m_dbg_lvl3 );
}

// SCHEDULINGS

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
ArmPostProcessModule<VAR, DIM, STN, AMR>::scheduleInitialize (
    SchedulerP & sched,
    const LevelP & level
)
{
    scheduleDoAnalysis_tip<AMR> ( sched, level );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
ArmPostProcessModule<VAR, DIM, STN, AMR>::scheduleDoAnalysis (
    SchedulerP & sched,
    const LevelP & level
)
{
    scheduleDoAnalysis_tip<AMR> ( sched, level );
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < ! MG, void >::type
ArmPostProcessModule<VAR, DIM, STN, AMR>::scheduleDoAnalysis_tip (
    SchedulerP & sched,
    const LevelP & level
)
{
    if ( arm_postproc )
    {
        Task * task = scinew Task ( "ArmPostProcessModule::task_do_analysis_tip", this,
                                    &ArmPostProcessModule::task_do_analysis_tip );
        task->requires ( Task::NewDW, d_psi_label, get_var<VAR>::ghost_type, 1 );
        task->requires ( Task::OldDW, tip_position_label );
        task->computes ( tip_position_label );
        task->computes ( tip_velocity_label );
        task->computes ( tip_curvatures_label );
        task->setType ( Task::OncePerProc );
        task->usesMPI ( true );

        sched->addTask ( task, sched->getLoadBalancer()->getPerProcessorPatchSet ( level ), d_app->getMaterialManagerP()->allMaterials() );
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
template < bool MG >
typename std::enable_if < MG, void >::type
ArmPostProcessModule<VAR, DIM, STN, AMR>::scheduleDoAnalysis_tip (
    SchedulerP & sched,
    const LevelP & level
)
{
    if ( !level->hasFinerLevel() && arm_postproc )
    {
        Task * task = scinew Task ( "ArmPostProcessModule::task_do_analysis_tip", this,
                                    &ArmPostProcessModule::task_do_analysis_tip );
        task->requires ( Task::NewDW, d_psi_label, get_var<VAR>::ghost_type, 1 );
        task->requires ( Task::OldDW, tip_position_label );
        task->modifies ( d_regridder->getRefineFlagLabel(), d_regridder->refineFlagMaterials() );
        task->computes ( tip_position_label );
        task->computes ( tip_velocity_label );
        task->computes ( tip_curvatures_label );
        task->setType ( Task::OncePerProc );
        task->usesMPI ( true );

        sched->addTask ( task, sched->getLoadBalancer()->getPerProcessorPatchSet ( level ), d_app->getMaterialManagerP()->allMaterials() );
    }
}

template<VarType VAR, DimType DIM, StnType STN, bool AMR>
void
ArmPostProcessModule<VAR, DIM, STN, AMR>::task_do_analysis_tip (
    const ProcessorGroup * myworld,
    const PatchSubset * patches,
    const MaterialSubset *,
    DataWarehouse * dw_old,
    DataWarehouse * dw_new
)
{
    int myrank = myworld->myRank();
    DOUT ( d_app->m_dbg_lvl1,  myrank << "==== ArmPostProcessModule::task_do_analysis_tip ====" );

    const Level * level = dw_new->getGrid()->getLevel( dw_new->getGrid()->numLevels() - 1 ).get_rep();
    ASSERT ( patches->empty() || level == getLevel( patches ) );
    arm_postproc->setLevel ( level );

    arm_postproc->initializeLocations();

    std::forward_list< std::tuple<IntVector, IntVector, std::reference_wrapper< View < ScalarField<const double> > >, std::shared_ptr< View< ScalarField<int> > > > > list;
    auto iter = list.before_begin();

    for ( int p = 0; p < patches->size(); ++p )
    {
        const Patch * patch = patches->get ( p );
        DOUT ( d_app->m_dbg_lvl2,  myrank << "== Patch: " << *patch << " Level: " << patch->getLevel()->getIndex() );;

        SubProblems < PureMetalProblem<VAR, STN> > subproblems ( dw_new, d_app->getSubProblemsLabel(), PM::material, patch );
        std::shared_ptr< View< ScalarField<int> > > refine_flag { AMR ? scinew DWView < ScalarField<int>, CC, DIM > ( dw_new, d_regridder->getRefineFlagLabel(), PM::material, patch ) : nullptr };

        for ( const auto & p : subproblems )
        {
            DOUT ( d_app->m_dbg_lvl3,  myrank << "= Iterating over " << p );;

            const IntVector & low = p.get_low();
            const IntVector & high = p.get_high();
            const std::list<Patch::FaceType> & faces = p.get_faces();
            View < ScalarField<const double> > & psi = p.template get_fd_view<PM::PSI> ( dw_new );

            arm_postproc->setLocations ( patch, low, high, faces, psi );

            std::stringstream ss;
            arm_postproc->printLocations ( ss << myrank << ": " );
            DOUTR ( dbg_arm_tip, ss.str().c_str() );

            iter = list.emplace_after ( iter, low, high, psi, refine_flag );
        }
    }

    std::stringstream ss;

    DOUT ( d_app->m_dbg_lvl2,  myrank );

    if ( dbg_arm_tip )
    {
        arm_postproc->printLocations ( ss << myrank << ": " );
        DOUTR ( dbg_arm_tip, ss.str().c_str() );
        ss.str ( "" );
        ss.clear();
    }

    arm_postproc->reduceLocations ( myworld );

    if ( dbg_arm_tip )
    {
        DOUT ( d_app->m_dbg_lvl2,  myrank );
        arm_postproc->printLocations ( ss << myrank << ": " );
        DOUTR ( dbg_arm_tip, ss.str().c_str() );
        ss.str ( "" );
        ss.clear();
    }

    arm_postproc->initializeData();

    for ( const auto & it : list )
    {
        arm_postproc->setData ( std::get<0> ( it ), std::get<1> ( it ), std::get<2> ( it ), std::get<3> ( it ).get() );
    }

    if ( dbg_arm_tip )
    {
        arm_postproc->printData ( ss << myrank << ": " );
        DOUTR ( dbg_arm_tip, ss.str().c_str() );
        ss.str ( "" );
        ss.clear();
    }

    arm_postproc->reduceData ( myworld );

    if ( myrank == 0 )
    {
        if ( dbg_arm_tip )
        {
            arm_postproc->printData ( ss << myrank << ": " );
            DOUTR ( dbg_arm_tip, ss.str().c_str() );
        }

        double tip_position;
        double tip_curvatures[3];
        arm_postproc->computeTipInfo ( tip_position, tip_curvatures );

        max_vartype old_position ( tip_position );
        if ( dw_old )
            dw_old->get ( old_position, tip_position_label );
        else
            old_position.setData ( tip_position );

        dw_new->put ( max_vartype ( tip_position ), tip_position_label );
        dw_new->put ( max_vartype ( ( tip_position - old_position ) / delt ), tip_velocity_label );
        dw_new->put ( maxvec_vartype ( { tip_curvatures[0], tip_curvatures[1], tip_curvatures[2] } ), tip_curvatures_label );
    }

    DOUT ( d_app->m_dbg_lvl2,  myrank );
}

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_tip_ArmPostProcessModule_h
