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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_GlobalData_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_GlobalData_h

#include <Core/Util/RefCounted.h>
#include <Core/Malloc/Allocator.h>
// #include <vector>

#include <HYPRE_sstruct_mv.h>

namespace Uintah
{
namespace HypreSStruct
{

class GlobalData
    : public RefCounted
{
    const int               m_nparts;   // # of total parts/levels
    int *                   m_nboxes;   // # of boxes/patches per part/level
    const int               m_nvars;    // tot # of variables on this part/level
    HYPRE_SStructVariable * m_vartypes; // variable types on this part/level

public:
    GlobalData (
        const int & nparts
    ) : RefCounted(),
        m_nparts ( nparts ),
        m_nboxes ( scinew int[nparts] ),
        m_nvars ( 1 ), // TODO generalize for systems
        m_vartypes ( scinew HYPRE_SStructVariable[nparts] )
    {
    };

    virtual ~GlobalData()
    {
        delete[] m_nboxes;
        delete[] m_vartypes;
    };

    const int &
    nParts (
    ) const
    {
        return m_nparts;
    };

    const int &
    nBoxes (
        const int & part
    ) const
    {
        return m_nboxes[part];
    };

    const int &
    nVars (
    ) const
    {
        return m_nvars;
    };

    const HYPRE_SStructVariable *
    varTypes (
    ) const
    {
        return m_vartypes;
    }

    void
    setPart (
        const int & part,
        const int & npatches
    )
    {
        m_nboxes[part] = npatches;
        m_vartypes[part] = HYPRE_SSTRUCT_VARIABLE_CELL;
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_GlobalData_h


