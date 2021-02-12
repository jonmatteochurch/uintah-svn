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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructPCG_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructPCG_h

#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_implementation.h>

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>

#include <HYPRE_krylov.h>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

template<int DIM, int C2F, bool P>
class sstruct_solver< ( int ) S::PCG, DIM, C2F, P>
    : virtual public sstruct_implementation<DIM, C2F>
{
protected:
    using sstruct = sstruct_implementation<DIM, C2F>;

    HYPRE_Solver * psolver, & solver, precond;
    HYPRE_Matrix * pA, & A;
    HYPRE_Vector * pb, & b;
    HYPRE_Vector * px, & x;

    bool & initialized, precond_initialized;
    using sstruct::guess_updated;

public: // STATIC MEMBERS
    static constexpr auto SetPrecond = HYPRE_PCGSetPrecond;

    static const HYPRE_PtrToSolverFcn precond_solve;
    static const HYPRE_PtrToSolverFcn precond_setup;

    sstruct_solver (
        const GlobalDataP & gdata
    ) : sstruct ( gdata ),
        psolver ( P ? & precond : reinterpret_cast<HYPRE_Solver *> ( & ( sstruct::solver ) ) ), solver ( *psolver ),
        pA ( reinterpret_cast<HYPRE_Matrix *> ( & ( sstruct::A ) ) ), A ( *pA ),
        pb ( reinterpret_cast<HYPRE_Vector *> ( & ( sstruct::b ) ) ), b ( *pb ),
        px ( reinterpret_cast<HYPRE_Vector *> ( & ( sstruct::x ) ) ), x ( *px ),
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

        IGNOREPARAM ( SStructPCG, params, csolver_type );
        IGNOREPARAM ( SStructPCG, params, max_levels );
        IGNOREPARAM ( SStructPCG, params, num_post_relax );
        IGNOREPARAM ( SStructPCG, params, num_pre_relax );
        IGNOREPARAM ( SStructPCG, params, relax_type );
        IGNOREPARAM ( SStructPCG, params, skip_relax );
        IGNOREPARAM ( SStructPCG, params, ssolver );
        IGNOREPARAM ( SStructPCG, params, weight );

        HYPRE_SStructPCGCreate ( comm, ( HYPRE_SStructSolver * ) psolver );

        if ( params->tol != -1 )
            HYPRE_PCGSetTol ( solver, params->tol );
        if ( params->abs_tol != -1 )
            HYPRE_PCGSetAbsoluteTol ( solver, params->abs_tol );
        if ( params->res_tol != -1 )
            HYPRE_PCGSetResidualTol ( solver, params->res_tol );
        if ( params->abstolf != -1 )
            HYPRE_PCGSetAbsoluteTolFactor ( solver, params->abstolf );
        if ( params->cf_tol != -1 )
            HYPRE_PCGSetConvergenceFactorTol ( solver, params->cf_tol );
        if ( params->max_iter != -1 )
            HYPRE_PCGSetMaxIter ( solver, params->max_iter );
        if ( params->two_norm != -1 )
            HYPRE_PCGSetTwoNorm ( solver, params->two_norm );
        if ( params->rel_change != -1 )
            HYPRE_PCGSetRelChange ( solver, params->rel_change );
        if ( params->recompute_residual != -1 )
            HYPRE_PCGSetRecomputeResidual ( solver, params->recompute_residual );
        if ( params->recompute_residual_p != -1 )
            HYPRE_PCGSetRecomputeResidualP ( solver, params->recompute_residual_p );
        if ( params->logging != -1 )
            HYPRE_PCGSetLogging ( solver, params->logging );

        initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        HYPRE_PCGSetup ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        HYPRE_PCGSolve ( solver, A, b, x );
        HYPRE_PCGGetNumIterations ( solver, &out->num_iterations );
        HYPRE_PCGGetFinalRelativeResidualNorm ( solver, &out->res_norm );
        guess_updated = false;
    }

    virtual void
    solverFinalize()
    override
    {
        if ( initialized ) HYPRE_SStructPCGDestroy ( ( HYPRE_SStructSolver ) solver  );
        initialized = false;
    }
};

template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::PCG, DIM, C2F, precond>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_PCGSolve;
template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::PCG, DIM, C2F, precond>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_PCGSetup;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructPCG_h
