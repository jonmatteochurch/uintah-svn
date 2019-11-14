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

#include <Core/Util/Factory/Implementation.h>

namespace Uintah
{
namespace HypreSStruct
{

template < S, int DIM, int C2F, P PCN = P::None, bool final = true > class SStructSolver;

template < S SLV, int DIM, int C2F, P PCN, bool final >
class SStructSolver
    : public SStructSolver<SLV, DIM, C2F, P::None, false>
    , public Implementation< SStructSolver<SLV, DIM, C2F, PCN>, SStructInterface, const GlobalDataP & >
{
    static_assert ( PCN != P::None, "Generic preconditioned SStructSolver implementation initialized for un-preconditioned SStructSolver" );
    static_assert ( final, "Generic preconditioned SStructSolver implementation must be final" );

    using PSolver = SStructSolver<SLV, DIM, C2F, P::None, false>;
    using Precond = SStructSolver< ( S ) PCN, DIM, C2F>;

    Precond precond_solver;

public: // STATIC MEMBERS

    /// Class name as used by ApplicationFactory
    static const std::string Name;

    SStructSolver (
        const GlobalDataP & gdata
    ) : PSolver ( gdata ),
        precond_solver ( gdata )
    {
    }

    virtual ~SStructSolver()
    {
        precondFinalize();
    }

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override
    {
        PSolver::solverInitialize ( comm, params );
        precondInitialize( comm, params );
    }

protected:

    virtual void
    solverFinalize (
    ) override
    {
        PSolver::solverFinalize();
        precondFinalize();
    }

protected:

    virtual void
    precondInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) override
    {
        precond_solver.solverInitialize ( comm, params );
        PSolver::SetPrecond ( PSolver::solver, Precond::precond, Precond::precond_setup, precond_solver );
    }

    virtual void
    precondFinalize (
    ) override
    {
        precond_solver.solverFinalize();
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartDataP_h


