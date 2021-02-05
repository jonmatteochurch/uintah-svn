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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructDiag_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructDiag_h

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

template<int DIM, int C2F>
class sstruct_solver< ( int ) P::Diagonal, DIM, C2F, true>
{
public: // STATIC MEMBERS
    static const HYPRE_PtrToSolverFcn precond_solve;
    static const HYPRE_PtrToSolverFcn precond_setup;

    sstruct_solver (
        const GlobalDataP & /*gdata*/
    )
    {}

    void
    solverInitialize (
        const MPI_Comm & /*comm*/,
        const SolverParams * /*params*/
    )
    {}

    void
    solverFinalize (
    )
    {}

    operator HYPRE_SStructSolver() { return nullptr; }
    operator HYPRE_Solver() { return nullptr; }
};

template<int DIM, int C2F> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) P::Diagonal, DIM, C2F, true>::precond_solve = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructDiagScale;
template<int DIM, int C2F> const HYPRE_PtrToSolverFcn sstruct_solver< ( int ) P::Diagonal, DIM, C2F, true>::precond_setup = ( HYPRE_PtrToSolverFcn ) HYPRE_SStructDiagScaleSetup;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstructDiag_h
