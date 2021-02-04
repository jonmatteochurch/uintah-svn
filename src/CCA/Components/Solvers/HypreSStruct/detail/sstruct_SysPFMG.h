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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_SysPFMG_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_SysPFMG_h

#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_implementation.h>

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>

#include <HYPRE_sstruct_ls.h>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

template<int DIM, int C2F, bool precond>
class sstruct_solver< ( int ) S::SysPFMG, DIM, C2F, precond>
    : public sstruct_implementation<DIM, C2F>
{
protected:
    using sstruct = sstruct_implementation<DIM, C2F>;

    using sstruct::gdata;

    using sstruct::solver;
    using sstruct::A;
    using sstruct::b;
    using sstruct::x;

    using sstruct::solver_initialized;
    using sstruct::guess_updated;

public: // STATIC MEMBERS
    static const HYPRE_PtrToSolverFcn precond_solve;
    static const HYPRE_PtrToSolverFcn precond_setup;

    sstruct_solver (
        const GlobalDataP & gdata
    ) : sstruct ( gdata )
    {
    }

    virtual ~sstruct_solver()
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
        ASSERTMSG ( gdata->nParts() == 1, "SysPFMG can be used only for problems with one grid part." );

        IGNOREPARAM ( SStructSysPFMG, params, csolver_type );
        IGNOREPARAM ( SStructSysPFMG, params, max_levels );
        IGNOREPARAM ( SStructSysPFMG, params, ssolver );
        IGNOREPARAM ( SStructSysPFMG, params, abs_tol );
        IGNOREPARAM ( SStructSysPFMG, params, two_norm );
        IGNOREPARAM ( SStructSysPFMG, params, recompute_residual );
        IGNOREPARAM ( SStructSysPFMG, params, recompute_residual_p );

        HYPRE_SStructSysPFMGCreate ( comm, &solver );

        if ( params->tol != -1 )
            HYPRE_SStructSysPFMGSetTol ( solver, params->tol );
        if ( params->max_iter != -1 )
            HYPRE_SStructSysPFMGSetMaxIter ( solver, params->max_iter );
        if ( params->rel_change != -1 )
            HYPRE_SStructSysPFMGSetRelChange ( solver, params->rel_change );
        if ( precond || !guess_updated )
            HYPRE_SStructSysPFMGSetZeroGuess ( solver );
        else
            HYPRE_SStructSysPFMGSetNonZeroGuess ( solver );
        if ( params->relax_type != -1 )
            HYPRE_SStructSysPFMGSetRelaxType ( solver, params->relax_type );
        if ( params->relax_type == WeightedJacobi && params->weight != -1 )
            HYPRE_SStructSysPFMGSetJacobiWeight ( solver, params->weight );
        if ( params->num_pre_relax != -1 )
            HYPRE_SStructSysPFMGSetNumPreRelax ( solver, params->num_pre_relax );
        if ( params->num_post_relax !=  -1 )
            HYPRE_SStructSysPFMGSetNumPostRelax ( solver, params->num_post_relax );
        if ( params->skip_relax != -1 )
            HYPRE_SStructSysPFMGSetSkipRelax ( solver, params->skip_relax );
        if ( params->logging != -1 )
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

template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::SysPFMG, DIM, C2F, precond>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSolve;
template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::SysPFMG, DIM, C2F, precond>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSysPFMGSetup;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructSysPFMG_h
