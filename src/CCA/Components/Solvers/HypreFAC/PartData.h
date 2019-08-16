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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreFAC_PartData_h
#define Packages_Uintah_CCA_Components_Solvers_HypreFAC_PartData_h

#include <Core/Util/DebugStream.h>
#include <Core/Util/RefCounted.h>

namespace Uintah
{
namespace HypreFAC
{
    
struct PartData : public RefCounted
{
    /* for GridSetExtents */
    int                             nboxes;     // tot # of boxes/patches on this part/level
    std::vector<IntVector>          ilowers;    // lower cell index per my_rank owned box on this part/level
    std::vector<IntVector>          iuppers;    // upper cell index per my_rank owned box on this part/level
    std::vector<int>                boxsizes;   // tot # of cells per my_rank owned box on this part/level
    std::vector<std::array<bool, 6>> interfaces; //

    std::vector<int>                box2patch;
    std::map<int, int>               patch2box;

    /* for amr structure */
    int plevel;         // depth of this part/level
    int prefinement[3]; // refinement factors of this part/level with respect to the immidiately coarser level

    PartData()
        : RefCounted()
        , nboxes ( 0 )
        , ilowers ( 0 )
        , iuppers ( 0 )
        , boxsizes ( 0 )
        , interfaces ( 0 )
        , box2patch ( 0 )
        , patch2box()
        , plevel ( -1 )
        , prefinement { -1, -1, 1}
    {
        cout_doing << " ### CREATED HypreFAC::PartData " << std::endl;
    };

    virtual ~PartData()
    {
        cout_doing << " ### DELETED HypreFAC::PartData " << std::endl;
    };
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreFAC_PartData_h



