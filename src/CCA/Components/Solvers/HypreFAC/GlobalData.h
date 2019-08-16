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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreFAC_GlobalData_h
#define Packages_Uintah_CCA_Components_Solvers_HypreFAC_GlobalData_h

#include <Core/Util/RefCounted.h>
#include <Core/Util/DebugStream.h>
#include <vector>

//#define HYPRE_TIMING

// #include <CCA/Components/Solvers/SolverCommon.h>
// #include <CCA/Components/Solvers/HypreFACSolverParams.h>
//
// #include <Core/Exceptions/InternalError.h>
// #include <Core/Grid/SimulationState.h>
// #include <Core/Util/Handle.h>
// #include <Core/Util/DebugStream.h>
// #include <Core/Util/Timers/Timers.hpp>
// #include <Core/Grid/Variables/PerPatch.h> // must be included after ProblemsP/AdditionalEntriesP where swapbytes override is defined
//
// #include <HYPRE_sstruct_ls.h>
//

namespace Uintah
{
namespace HypreFAC
{

extern DebugStream cout_doing;

struct GlobalData : public RefCounted
{
    int nparts; // # of total parts/levels
    std::vector<std::vector<int>> nboxes; // # of boxes/patches per part/level per mpi rank

    int                                 nvars;    // tot # of variables on this part/level
    std::vector<HYPRE_SStructVariable>  vartypes; // variable types on this part/level

    GlobalData ()
        : RefCounted()
        , nparts ( 0 )
        , nboxes ( 0 )
        , nvars ( 0 )
        , vartypes ( 0 )
    {
        cout_doing << " ### CREATED HypreFAC::GlobalData " << std::endl;
    };

    virtual ~GlobalData()
    {
        cout_doing << " ### DELETED HypreFAC::GlobalData " << std::endl;
    };
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreFAC_GlobalData_h


