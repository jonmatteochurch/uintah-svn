/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h
#define Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h

//#define HYPRE_TIMING

// #include <CCA/Components/Solvers/SolverCommon.h>
#include <CCA/Components/Solvers/HypreFAC/PartDataP.h>
#include <CCA/Components/Solvers/HypreFAC/GlobalDataP.h>
//
// #include <Core/Exceptions/InternalError.h>
// #include <Core/Grid/SimulationState.h>
// #include <Core/Util/Handle.h>
#include <Core/Util/RefCounted.h>
#include <Core/Util/DebugStream.h>
// #include <Core/Util/Timers/Timers.hpp>
// #include <Core/Grid/Variables/PerPatch.h> // must be included after ProblemsP/AdditionalEntriesP where swapbytes override is defined
//
#include <HYPRE_sstruct_ls.h>
//
// #include <iostream>

namespace Uintah
{
namespace HypreFAC
{

extern DebugStream cout_doing;

struct SolverStruct : public RefCounted
{
    GlobalDataP  data;   // global data
    PartDataP  * pdatas; // local data per part/level

    bool                    created;
    bool                    restart;
    HYPRE_SStructGrid   *   grid;
    HYPRE_SStructStencil  * stencil;
    HYPRE_SStructGraph   *  graph;
    HYPRE_SStructMatrix  *  A;
    HYPRE_SStructVector  *  b;
    HYPRE_SStructVector  *  x;

    SolverStruct()
        : RefCounted()
        , data ( nullptr )
        , pdatas ( nullptr )
        , created ( false )
        , restart ( true )
        , grid ( nullptr )
        , stencil ( nullptr )
        , graph ( nullptr )
        , A ( nullptr )
        , b ( nullptr )
        , x ( nullptr )
    {
        cout_doing << " ### CREATED HypreFAC::SolverStruct " << std::endl;
    };

    virtual ~SolverStruct()
    {
        if ( created )
        {
            HYPRE_SStructVectorDestroy ( *x );
            HYPRE_SStructVectorDestroy ( *b );
            HYPRE_SStructMatrixDestroy ( *A );
            HYPRE_SStructGraphDestroy ( *graph );
            HYPRE_SStructStencilDestroy ( *stencil );
            HYPRE_SStructGridDestroy ( *grid );
        }

//         if ( data )
//         {
//             delete data;
//             data = nullptr;
//         }
        if ( pdatas )
        {
            delete[] pdatas;
            data = nullptr;
        }
        if ( grid )
        {
            delete grid;
            grid = nullptr;
        }
        if ( stencil )
        {
            delete stencil;
            stencil = nullptr;
        }
        if ( graph )
        {
            delete graph;
            graph = nullptr;
        }
        if ( A )
        {
            delete A;
            A = nullptr;
        }
        if ( b )
        {
            delete b;
            b = nullptr;
        }
        if ( x )
        {
            delete x;
            x = nullptr;
        }

        cout_doing << " ### DELETED hypre_fac_solver_struct " << std::endl;
    };
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h



