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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStrucPCG_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStrucPCG_h

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructPSolver.h>

#include <HYPRE_krylov.h>

namespace Uintah
{
namespace HypreSStruct
{

template<int DIM, int C2F, bool final>
class SStructSolver<S::PCG, DIM, C2F, P::None, final>
    : public SStructImplementation<DIM, C2F>
    , public std::conditional < final,
    Implementation< SStructSolver<S::PCG, DIM, C2F, P::None, final>, SStructInterface, const GlobalDataP & >,
    SStructPSolver
    >::type
{
protected:

    using SStruct = SStructImplementation<DIM, C2F>;

    HYPRE_Solver solver;
    HYPRE_Matrix A;
    HYPRE_Vector b;
    HYPRE_Vector x;

    using SStruct::solver_initialized;
    using SStruct::guess_updated;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    static constexpr auto SetPrecond = HYPRE_PCGSetPrecond;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStruct ( gdata ),
        solver ( ( HYPRE_Solver ) SStruct::solver ),
        A ( ( HYPRE_Matrix ) SStruct::A ),
        b ( ( HYPRE_Vector ) SStruct::b ),
        x ( ( HYPRE_Vector ) SStruct::x )
    {
        std::cout << "construct" << __FILE__ << ":" << __LINE__ << std::endl;
    }

    virtual ~SStructSolver()
    {
        std::cout << "destruct" << __FILE__ << ":" << __LINE__ << std::endl;
        solverFinalize();
    }

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override
    {
        std::cout << "solverInitialize" << __FILE__ << ":" << __LINE__ << std::endl;
        ASSERT ( !solver_initialized );

        HYPRE_SStructPCGCreate ( comm, & ( SStruct::solver ) );

        if ( params->tol > 0 )
            HYPRE_PCGSetTol ( solver, params->tol );
        if ( params->abs_tol > 0 )
            HYPRE_PCGSetAbsoluteTol ( solver, params->abs_tol );
        if ( params->res_tol > 0 )
            HYPRE_PCGSetResidualTol ( solver, params->res_tol );
        if ( params->max_iter > 0 )
            HYPRE_PCGSetMaxIter ( solver, params->max_iter );
        if ( params->two_norm > -1 )
            HYPRE_PCGSetTwoNorm ( solver, params->two_norm );
        if ( params->rel_change > -1 )
            HYPRE_PCGSetRelChange ( solver, params->rel_change );
        if ( params->recompute_residual > -1 )
            HYPRE_PCGSetRecomputeResidual ( solver, params->recompute_residual );
        if ( params->recompute_residual_p > -1 )
            HYPRE_PCGSetRecomputeResidualP ( solver, params->recompute_residual_p );
        if ( params->logging >  -1 )
            HYPRE_PCGSetLogging ( solver, params->logging );

        solver_initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        std::cout << "solverUpdate" << __FILE__ << ":" << __LINE__ << std::endl;
        HYPRE_PCGSetup ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        std::cout << "solve" << __FILE__ << ":" << __LINE__ << std::endl;
        HYPRE_PCGSolve ( solver, A, b, x );
        HYPRE_PCGGetNumIterations ( solver, &out->num_iterations );
        HYPRE_PCGGetFinalRelativeResidualNorm ( solver, &out->res_norm );
        HYPRE_PCGGetConverged ( solver, &out->converged );
        HYPRE_PCGGetResidual ( solver, &out->residual );
        guess_updated = false;
    }

protected:

    virtual void
    solverFinalize()
    override
    {
        std::cout << "solverFinalize" << __FILE__ << ":" << __LINE__ << std::endl;
        if ( solver_initialized ) HYPRE_SStructPCGDestroy ( SStruct::solver );
        solver_initialized = false;
    }
};

template<int DIM, int C2F>
class SStructSolver<S::PCG, DIM, C2F, P::Diag>
    : public SStructSolver<S::PCG, DIM, C2F>
    , public Implementation< SStructSolver<S::PCG, DIM, C2F, P::Diag>, SStructInterface, const GlobalDataP & >
{
    using PSolver = SStructSolver<S::PCG, DIM, C2F, P::None, false>;

    HYPRE_Solver precond_solver;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStructSolver<S::PCG, DIM, C2F> ( gdata ),
        precond_solver ( nullptr )
    {
        std::cout << "construct" << __FILE__ << ":" << __LINE__ << std::endl;
    }

    virtual ~SStructSolver()
    {
        std::cout << "destruct" << __FILE__ << ":" << __LINE__ << std::endl;
    }

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override
    {
        std::cout << "solverInitialize" << __FILE__ << ":" << __LINE__ << std::endl;
        PSolver::solverInitialize();
        precondInitialize();
    }

protected:

    virtual void 
    precondInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override

    {
        std::cout << "precondFinalize" << __FILE__ << ":" << __LINE__ << std::endl;
        HYPRE_PCGSetPrecond ( PSolver::solver,
                              ( HYPRE_PtrToSolverFcn ) HYPRE_SStructDiagScale,
                              ( HYPRE_PtrToSolverFcn ) HYPRE_SStructDiagScaleSetup,
                              precond_solver );
    }

};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartDataP_h


