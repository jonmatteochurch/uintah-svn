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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h
#define Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h

#include <CCA/Components/Solvers/HypreFAC/PartDataP.h>
#include <CCA/Components/Solvers/HypreFAC/GlobalDataP.h>
#include <Core/Util/RefCounted.h>
#include <HYPRE_sstruct_ls.h>

namespace Uintah
{
namespace HypreFAC
{

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

        if ( pdatas )
        {
            delete[] pdatas;
            pdatas = nullptr;
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
    };
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreFACSolverStruct_h



