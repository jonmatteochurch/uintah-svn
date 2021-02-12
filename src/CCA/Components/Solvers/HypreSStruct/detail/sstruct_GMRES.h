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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructGMRES_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructGMRES_h

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
class sstruct_solver< ( int ) S::GMRES, DIM, C2F, P>
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
    static constexpr auto SetPrecond = HYPRE_GMRESSetPrecond;

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

        IGNOREPARAM ( SStructGMRES, params, csolver_type );
        IGNOREPARAM ( SStructGMRES, params, max_levels );
        IGNOREPARAM ( SStructGMRES, params, num_post_relax );
        IGNOREPARAM ( SStructGMRES, params, num_pre_relax );
        IGNOREPARAM ( SStructGMRES, params, relax_type );
        IGNOREPARAM ( SStructGMRES, params, skip_relax );
        IGNOREPARAM ( SStructGMRES, params, ssolver );
        IGNOREPARAM ( SStructGMRES, params, weight );

        HYPRE_SStructGMRESCreate ( comm, ( HYPRE_SStructSolver * ) psolver );

        if ( params->tol != -1 )
            HYPRE_GMRESSetTol ( solver, params->tol );
        if ( params->abs_tol != -1 )
            HYPRE_GMRESSetAbsoluteTol ( solver, params->abs_tol );
        if ( params->cf_tol != -1 )
            HYPRE_GMRESSetConvergenceFactorTol ( solver, params->cf_tol );
        if ( params->min_iter != -1 )
            HYPRE_GMRESSetMinIter ( solver, params->min_iter );
        if ( params->max_iter != -1 )
            HYPRE_GMRESSetMaxIter ( solver, params->max_iter );
        if ( params->two_norm != -1 )
            HYPRE_GMRESSetKDim ( solver, params->k_dim );
        if ( params->rel_change != -1 )
            HYPRE_GMRESSetRelChange ( solver, params->rel_change );
        if ( params->skip_real_r_check != -1 )
            HYPRE_GMRESSetSkipRealResidualCheck ( solver, params->skip_real_r_check );
        if ( params->logging != -1 )
            HYPRE_GMRESSetLogging ( solver, params->logging );

        initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        HYPRE_GMRESSetup ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        HYPRE_GMRESSolve ( solver, A, b, x );
        HYPRE_GMRESGetNumIterations ( solver, &out->num_iterations );
        HYPRE_GMRESGetFinalRelativeResidualNorm ( solver, &out->res_norm );
        guess_updated = false;
    }

    virtual void
    solverFinalize()
    override
    {
        if ( initialized ) HYPRE_SStructGMRESDestroy ( ( HYPRE_SStructSolver ) solver );
        initialized = false;
    }
};

template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::GMRES, DIM, C2F, precond>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_GMRESSolve;
template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::GMRES, DIM, C2F, precond>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_GMRESSetup;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructGMRES_h
