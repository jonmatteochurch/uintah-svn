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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructFAC_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructFAC_h

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>

#include <HYPRE_sstruct_mv.h>

#ifdef PRINTSYSTEM
#   include "/home/jonmatteochurch/Developer/hypre/fork/src/sstruct_ls/fac.h"
#endif

namespace Uintah
{
namespace HypreSStruct
{

template<int DIM, int C2F >
class SStructSolver<S::FAC, DIM, C2F>
    : public SStructImplementation<DIM, C2F>
    , public Implementation< SStructSolver<S::FAC, DIM, C2F>, SStructInterface, const GlobalDataP & >
{
protected:
    using SStruct = SStructImplementation<DIM, C2F>;

    using SStruct::gdata;
    using SStruct::plevels;
    using SStruct::prefinements;

    using SStruct::solver;
    using SStruct::A;
    using SStruct::b;
    using SStruct::x;

    using SStruct::solver_initialized;
    using SStruct::guess_updated;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    static constexpr auto precond = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructFACSolve3;
    static constexpr auto precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructFACSetup2;

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

        const auto & nparts = gdata->nParts();

        HYPRE ( SStructFACCreate ) ( comm, &solver );
        HYPRE ( SStructFACSetPLevels ) ( solver, nparts, plevels );
        HYPRE ( SStructFACSetPRefinements ) ( solver, nparts, prefinements );

        if ( params->max_levels > 0 )
            HYPRE ( SStructFACSetMaxLevels ) ( solver, params->max_levels );
        if ( params->tol > 0 )
            HYPRE ( SStructFACSetTol ) ( solver, params->tol );
        if ( params->max_iter > 0 )
            HYPRE ( SStructFACSetMaxIter ) ( solver, params->max_iter );
        if ( params->rel_change > -1 )
            HYPRE ( SStructFACSetRelChange ) ( solver, params->rel_change );
        if ( !guess_updated )
            HYPRE ( SStructFACSetZeroGuess ) ( solver );
        else
            HYPRE ( SStructFACSetNonZeroGuess ) ( solver );
        if ( params->relax_type >  -1 )
            HYPRE ( SStructFACSetRelaxType ) ( solver, params->relax_type );
        if ( params->relax_type == WeightedJacobi )
            HYPRE_SStructFACSetJacobiWeight ( solver, params->weight );
        if ( params->num_pre_relax >  -1 )
            HYPRE ( SStructFACSetNumPreRelax ) ( solver, params->num_pre_relax );
        if ( params->num_post_relax >  -1 )
            HYPRE ( SStructFACSetNumPostRelax ) ( solver, params->num_post_relax );
        if ( params->csolver_type >  -1 )
            HYPRE ( SStructFACSetCoarseSolverType ) ( solver, params->csolver_type );
        if ( params->logging >  -1 )
            HYPRE ( SStructFACSetLogging ) ( solver, params->logging );

        solver_initialized = true;
    }

    virtual void
    solverUpdate (
    ) override
    {
        HYPRE ( SStructFACSetup2 ) ( solver, A, b, x );
    }

    virtual void
    solve (
        SolverOutput * out
    ) override
    {
        HYPRE ( SStructFACSolve3 ) ( solver, A, b, x );
        HYPRE ( SStructFACGetNumIterations ) ( solver, &out->num_iterations );
        HYPRE ( SStructFACGetFinalRelativeResidualNorm ) ( solver, &out->res_norm );
        guess_updated = false;
    }

#ifdef PRINTSYSTEM
    virtual void
    printSystem (
        std::string * fname
    ) override
    {
        SStruct::printSystem ( fname );

        std::string name = "fac_" + fname[0];
        HYPRE_SStructMatrix fac_A = ( ( hypre_FACData * ) solver )->A_rap;
        HYPRE_SStructMatrixPrint ( name.c_str(), fac_A, 0 );
    }
#endif

public: // exposing this for when used as precond

    virtual void
    solverFinalize()
    override
    {
        if ( solver_initialized ) HYPRE ( SStructFACDestroy2 ) ( solver );
    }

public: // required if precond

    operator HYPRE_Solver ()
    {
        return ( HYPRE_Solver ) solver;
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructFAC_h


