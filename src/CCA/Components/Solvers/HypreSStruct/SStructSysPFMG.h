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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructSysPFMG_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructSysPFMG_h

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>

#include <HYPRE_sstruct_mv.h>

namespace Uintah
{
namespace HypreSStruct
{

template<int DIM, int C2F>
class SStructSolver<S::SysPFMG, DIM, C2F>
    : public SStructImplementation<DIM, C2F>
    , public Implementation< SStructSolver<S::SysPFMG, DIM, C2F>, SStructInterface, const GlobalDataP & >
{
public:
    using SStruct = SStructImplementation<DIM, C2F>;

    using SStruct::solver;
    using SStruct::A;
    using SStruct::b;
    using SStruct::x;

    using SStruct::solver_initialized;
    using SStruct::guess_updated;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    static constexpr auto precond = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSolve;
    static constexpr auto precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSetup;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStruct ( gdata )
    {
    }

    virtual ~SStructSolver()
    {
        solverFinalize();
    }

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override
    {
        ASSERT ( !solver_initialized );

        HYPRE_SStructSysPFMGCreate ( comm, &solver );

        if ( params->tol > 0 )
            HYPRE_SStructSysPFMGSetTol ( solver, params->tol );
        if ( params->max_iter > 0 )
            HYPRE_SStructSysPFMGSetMaxIter ( solver, params->max_iter );
        if ( params->rel_change > -1 )
            HYPRE_SStructSysPFMGSetRelChange ( solver, params->rel_change );
        if ( !guess_updated )
            HYPRE_SStructSysPFMGSetZeroGuess ( solver );
        else
            HYPRE_SStructSysPFMGSetNonZeroGuess ( solver );
        if ( params->relax_type >  -1 )
            HYPRE_SStructSysPFMGSetRelaxType ( solver, params->relax_type );
        if ( params->relax_type == WeightedJacobi )
            HYPRE_SStructSysPFMGSetJacobiWeight ( solver, params->weight );
        if ( params->num_pre_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPreRelax ( solver, params->num_pre_relax );
        if ( params->num_post_relax >  -1 )
            HYPRE_SStructSysPFMGSetNumPostRelax ( solver, params->num_post_relax );
        if ( params->skip_relax >  -1 )
            HYPRE_SStructSysPFMGSetSkipRelax ( solver, params->skip_relax );
        if ( params->logging >  -1 )
            HYPRE_SStructSysPFMGSetLogging ( solver, params->logging );

        solver_initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        HYPRE_SStructSysPFMGSetup ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        HYPRE_SStructSysPFMGSolve ( solver, A, b, x );
        HYPRE_SStructSysPFMGGetNumIterations ( solver, &out->num_iterations );
        HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm ( solver, &out->res_norm );
        guess_updated = false;
    }

public: // exposing this for when used as precond

    virtual void
    solverFinalize()
    override
    {
        if ( solver_initialized ) HYPRE_SStructSysPFMGDestroy ( solver );
        solver_initialized = false;
    }

public: // required if precond

    operator HYPRE_Solver ()
    {
        return ( HYPRE_Solver ) solver;
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructSysPFMG_h


