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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructSolver_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructSolver_h

#include <CCA/Components/Solvers/HypreSStruct/Definitions.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructInterface.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_solver.h>

#include <Core/Util/Factory/Implementation.h>


namespace Uintah
{
namespace HypreSStruct
{

template < S SLV, int DIM, int C2F, P PCN = P::None> class SStructSolver;

template < S SLV, int DIM, int C2F>
class SStructSolver<SLV, DIM, C2F, P::None> final
    : public detail::sstruct_ssolver<SLV, DIM, C2F>
    , Implementation< SStructSolver<SLV, DIM, C2F, P::None>, SStructInterface, const GlobalDataP & >
{
    using SStruct = detail::sstruct_implementation<DIM, C2F>;
    using SSolver = detail::sstruct_ssolver<SLV, DIM, C2F>;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStruct ( gdata ),
        SSolver ( gdata )
    {}

    virtual ~SStructSolver() = default;
};

template < S SLV, int DIM, int C2F, P PCN >
class SStructSolver final
    : public detail::sstruct_ssolver<SLV, DIM, C2F>
    , public detail::sstruct_precond<PCN, DIM, C2F>
    , Implementation< SStructSolver<SLV, DIM, C2F, PCN>, SStructInterface, const GlobalDataP & >
{
    static_assert ( PCN != P::None, "Generic preconditioned SStructSolver implementation initialized for un-preconditioned SStructSolver" );
    static_assert ( PCN != P::Diagonal, "Generic preconditioned SStructSolver implementation initialized for un-preconditioned SStructSolver" );
    static_assert ( (int) SLV > 0, "Generic preconditioned SStructSolver implementation initialized for SStructSolver which cannot be preconditioned" );

    using SStruct = detail::sstruct_implementation<DIM, C2F>;
    using SSolver = detail::sstruct_ssolver<SLV, DIM, C2F>;
    using Precond = detail::sstruct_precond<PCN, DIM, C2F>;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStruct ( gdata ),
        SSolver ( gdata ),
        Precond ( gdata )
    {
    }

    virtual ~SStructSolver()
    {
    }

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) final
    {
        Precond::solverInitialize ( comm, params );
        SSolver::solverInitialize ( comm, params );
        SSolver::SetPrecond ( SSolver::solver, Precond::precond_solve, Precond::precond_setup, (HYPRE_Solver) Precond::solver );
    }

    virtual void
    solverUpdate (
    ) final
    {
        return SSolver::solverUpdate();
    }

    virtual void
    solve (
        SolverOutput * out
    ) final
    {
        return SSolver::solve( out );
    }

    virtual void
    solverFinalize (
    ) final
    {
        Precond::solverFinalize();
        SSolver::solverFinalize();
    }
};

template < S SLV, int DIM, int C2F >
class SStructSolver<SLV, DIM, C2F, P::Diagonal> final
    : public detail::sstruct_ssolver<SLV, DIM, C2F>
    , Implementation< SStructSolver<SLV, DIM, C2F, P::Diagonal>, SStructInterface, const GlobalDataP & >
{
    static_assert ( (int) SLV > 0, "Generic preconditioned SStructSolver implementation initialized for solver type  that cannot be preconditioned" );

    using SStruct = detail::sstruct_implementation<DIM, C2F>;
    using SSolver = detail::sstruct_ssolver<SLV, DIM, C2F>;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    SStructSolver (
        const GlobalDataP & gdata
    ) : SStruct ( gdata ),
        SSolver ( gdata )
    {
    }

    virtual ~SStructSolver()
    {
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartDataP_h
