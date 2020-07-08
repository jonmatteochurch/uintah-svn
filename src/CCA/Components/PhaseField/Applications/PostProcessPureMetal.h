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

#ifndef Packages_Uintah_CCA_Components_PhaseField_Applications_PostProcessPureMetal_h
#define Packages_Uintah_CCA_Components_PhaseField_Applications_PostProcessPureMetal_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/DataTypes/PureMetalProblem.h>
#include <CCA/Components/PhaseField/Applications/Application.h>
#include <CCA/Components/PostProcessUda/PostProcessUda.h>
#include <vector>

namespace Uintah
{
namespace PhaseField
{

template<VarType VAR, DimType DIM, StnType STN, bool AMR = false>
class PostProcessPureMetal
    : public Application< PureMetalProblem<VAR, STN>, AMR >
    , public Implementation < PostProcessPureMetal<VAR, DIM, STN, AMR>, UintahParallelComponent, const ProcessorGroup *, const MaterialManagerP, const std::string &, int >
{

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    const static FactoryString Name;

    PostProcessUda postproc;

public:
    PostProcessPureMetal (
        const ProcessorGroup * myworld,
        const MaterialManagerP material_manager,
        const std::string & uda,
        int verbosity
    ) : Application< PureMetalProblem<VAR, STN>, AMR > ( myworld, material_manager, verbosity ),
        postproc ( myworld, material_manager, uda )
    {};

    virtual ~PostProcessPureMetal() {};

    virtual inline void
    problemSetup (
        const ProblemSpecP & params,
        const ProblemSpecP & restart_prob_spec,
        GridP & grid
    ) override
    {
        postproc.problemSetup ( params, restart_prob_spec, grid );
    }

    virtual inline void
    scheduleInitialize (
        const LevelP & level,
        SchedulerP & sched
    ) override
    {
        postproc.scheduleInitialize ( level, sched );
    }

    virtual inline void
    scheduleRestartInitialize (
        const LevelP & level,
        SchedulerP & sched
    ) override
    {
        postproc.scheduleRestartInitialize ( level, sched );
    }

    virtual inline void
    scheduleComputeStableTimeStep (
        const LevelP & level,
        SchedulerP & sched
    ) override
    {
        postproc.scheduleComputeStableTimeStep ( level, sched );
    }
};

} // namespace PhaseField
} // namespace Uintah

#endif
