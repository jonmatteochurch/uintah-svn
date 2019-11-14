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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructInterface_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructInterface_h

#include <CCA/Components/Solvers/HypreSStruct/GlobalDataP.h>
#include <CCA/Components/Solvers/HypreSStruct/PartDataP.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Grid/Variables/CCVariable.h>

namespace Uintah
{
namespace HypreSStruct
{

class SolverParams;
class SolverOutput;
class AdditionalEntries;

class SStructInterface
    : public RefCounted
{
public:
    virtual const GlobalDataP &
    globalData (
    ) const = 0;

    virtual const PartDataP &
    partData (
        const int & part
    ) const = 0;

    virtual bool
    isPartSet (
        const int & part
    ) const = 0;

    virtual void
    setPart (
        const int & part,
        const int & npatches,
        const int & nboxes,
        const int low[3],
        const int high[3],
        const int refinement[3],
        const int periodic[3]
    ) = 0;

    virtual void
    addBox (
        const int & part,
        const int & id,
        const int ilower[3],
        const int iupper[3],
        const bool interface[6]
    ) = 0;

    virtual void
    gridInitialize (
        const MPI_Comm & comm
    ) = 0;

    virtual void
    stencilInitialize (
        const MPI_Comm & comm
    ) = 0;

    virtual void
    graphInitialize (
        const MPI_Comm & comm,
        const GridP & grd,
        AdditionalEntries **** additional_entries
    ) = 0;

    virtual void
    matrixInitialize (
        const MPI_Comm & comm
    ) = 0;;

    virtual void
    rhsInitialize (
        const MPI_Comm & comm
    ) = 0;

    virtual void
    solutionInitialize (
        const MPI_Comm & comm
    ) = 0;

    virtual void
    solverInitialize (
        const MPI_Comm & comm,
        const SolverParams * params
    ) = 0;

    virtual void
    matrixUpdate (
        constCCVariable<Stencil7> *** stencil_entries,
        AdditionalEntries **** additional_entries
    ) = 0;

    virtual void
    rhsUpdate (
        constCCVariable<double> *** rhs
    ) = 0;

    virtual void
    guessUpdate (
        constCCVariable<double> *** guess
    ) = 0;

    virtual void
    assemble (
    ) = 0;

    virtual void
    solverUpdate (
    ) = 0;

#ifdef PRINTSYSTEM
    virtual void
    printSystem (
        std::string * fname
    ) = 0;
#endif

    virtual void
    solve (
        SolverOutput * out
    ) = 0;

    virtual void
    getSolution (
        CCVariable<double> *** solution
    ) = 0;

    virtual void
    finalize (
    ) = 0;

    virtual bool
    restart (
    ) = 0;
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructInterface_h




