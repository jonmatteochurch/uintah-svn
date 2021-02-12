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

#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_implementation.h>

#include <CCA/Components/Solvers/HypreSStruct/SolverParams.h>
#include <CCA/Components/Solvers/HypreSStruct/SolverOutput.h>

#include <Core/Exceptions/ProblemSetupException.h>

#include <HYPRE_sstruct_ls.h>

#ifdef PRINTSYSTEM
#   include "/home/jonmatteochurch/Developer/hypre/fork/src/sstruct_ls/_hypre_sstruct_ls.h"
#   include "/home/jonmatteochurch/Developer/hypre/fork/src/sstruct_ls/fac.h"
#endif

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

template<int DIM, int C2F, bool P>
class sstruct_solver< ( int ) S::FAC, DIM, C2F, P>
    : virtual public sstruct_implementation<DIM, C2F>
{
protected:
    using sstruct = sstruct_implementation<DIM, C2F>;

    using sstruct::gdata;
    using sstruct::plevels;
    using sstruct::prefinements;

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

        const auto & nparts = gdata->nParts();

        HYPRE ( SStructFACCreate ) ( comm, &solver );
        HYPRE ( SStructFACSetPLevels ) ( solver, nparts, plevels );
        HYPRE ( SStructFACSetPRefinements ) ( solver, nparts, prefinements );

        if ( params->max_levels > 0 )
            HYPRE ( SStructFACSetMaxLevels ) ( solver, params->max_levels );
        else
            SCI_THROW ( ProblemSetupException ( "HypreSStruct FAC solver: missing required parameter 'max_levels' cannot be deduced", __FILE__, __LINE__ ) );
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

        initialized = true;
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
        sstruct::printSystem ( fname );

        std::string name = "fac_" + fname[0];
        HYPRE_SStructMatrix fac_A = ( ( hypre_FACData * ) solver )->A_rap;
        HYPRE_SStructMatrixPrint ( name.c_str(), fac_A, 0 );
    }
#endif

    virtual void
    solverFinalize()
    override
    {
        if ( initialized ) HYPRE ( SStructFACDestroy2 ) ( solver );
        initialized = false;
    }
};

template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::FAC, DIM, C2F, precond>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructFACSolve3;
template<int DIM, int C2F, bool precond> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) S::FAC, DIM, C2F, precond>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructFACSetup2;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructFAC_h
