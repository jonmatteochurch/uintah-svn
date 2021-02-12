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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_Split_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_Split_h

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

template<int DIM, int C2F, bool P>
class sstruct_solver< ( int ) S::Split, DIM, C2F, P>
    : virtual public sstruct_implementation<DIM, C2F>
{
protected:
    using sstruct = sstruct_implementation<DIM, C2F>;

    HYPRE_SStructSolver & solver, precond;
    using sstruct::A;
    using sstruct::b;
    using sstruct::x;

    bool & initialized, precond_initialized;
    using sstruct::guess_updated;

public: // STATIC MEMBERS
    static const HYPRE_PtrToSolverFcn precond_solve;
    static const HYPRE_PtrToSolverFcn precond_setup;

    sstruct_solver (
        const GlobalDataP & gdata
    ) : sstruct ( gdata ),
        solver ( P ? precond : sstruct::solver ),
        initialized ( P ? precond_initialized : sstruct::solver_initialized ),
        precond_initialized ( false )
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
        ASSERT ( !initialized );

        IGNOREPARAM ( SStructSplit, params, csolver_type );
        IGNOREPARAM ( SStructSplit, params, max_levels );
        IGNOREPARAM ( SStructSplit, params, num_post_relax );
        IGNOREPARAM ( SStructSplit, params, num_pre_relax );
        IGNOREPARAM ( SStructSplit, params, rel_change );
        IGNOREPARAM ( SStructSplit, params, relax_type );
        IGNOREPARAM ( SStructSplit, params, skip_relax );
        IGNOREPARAM ( SStructSplit, params, weight );
        IGNOREPARAM ( SStructSplit, params, abs_tol );
        IGNOREPARAM ( SStructSplit, params, two_norm );
        IGNOREPARAM ( SStructSplit, params, recompute_residual );
        IGNOREPARAM ( SStructSplit, params, recompute_residual_p );
        IGNOREPARAM ( SStructSplit, params, min_iter );
        IGNOREPARAM ( SStructSplit, params, logging );

        HYPRE_SStructSplitCreate ( comm, &solver );

        if ( params->tol != -1 )
            HYPRE_SStructSplitSetTol ( solver, params->tol );
        if ( params->max_iter != -1 )
            HYPRE_SStructSplitSetMaxIter ( solver, params->max_iter );
        if ( P || !guess_updated )
            HYPRE_SStructSplitSetZeroGuess ( solver );
        else
            HYPRE_SStructSplitSetNonZeroGuess ( solver );
        if ( params->ssolver != -1 )
            HYPRE_SStructSplitSetStructSolver ( solver, params->ssolver );

        initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        HYPRE_SStructSplitSetup ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        HYPRE_SStructSplitSolve ( solver, A, b, x );
        HYPRE_SStructSplitGetNumIterations ( solver, &out->num_iterations );
        HYPRE_SStructSplitGetFinalRelativeResidualNorm ( solver, &out->res_norm );
        guess_updated = false;
    }

    virtual void
    solverFinalize()
    override
    {
        if ( initialized ) HYPRE_SStructSplitDestroy ( solver );
        initialized = false;
    }
};

template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::Split, DIM, C2F, precond>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSplitSolve;
template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::Split, DIM, C2F, precond>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructSplitSetup;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_Split_h
