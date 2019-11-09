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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartData_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartData_h

#include <Core/Util/RefCounted.h>
#include <Core/Malloc/Allocator.h>

#include <map>
#include <vector>

namespace Uintah
{
namespace HypreSStruct
{

class PartData
    : public RefCounted
{
    int ( * m_ilowers ) [HYPRE_MAXDIM];    // lower cell index per my_rank owned box on this part/level
    int ( * m_iuppers ) [HYPRE_MAXDIM];    // upper cell index per my_rank owned box on this part/level
    int * m_boxsizes;   // tot # of cells per my_rank owned box on this part/level
    bool ( * m_interfaces ) [2 * HYPRE_MAXDIM]; //

    std::vector<int> m_box2patch;
    std::map<int, int> m_patch2box;

    const int m_low[HYPRE_MAXDIM];
    const int m_high[HYPRE_MAXDIM];
    const int m_periodic[HYPRE_MAXDIM];
    const double m_prefinement_weight;

public:
    PartData (
        const int & nboxes,
        const int low[HYPRE_MAXDIM],
        const int high[HYPRE_MAXDIM],
        const int refinement[HYPRE_MAXDIM],
        const int periodic[HYPRE_MAXDIM]
    ) : RefCounted(),
        m_ilowers ( scinew int[nboxes][HYPRE_MAXDIM] ),
        m_iuppers ( scinew int[nboxes][HYPRE_MAXDIM] ),
        m_boxsizes ( scinew int [nboxes] ),
        m_interfaces ( scinew bool[nboxes][2 * HYPRE_MAXDIM] ),
        m_low
        {
            low[0] * refinement[0],
            low[1] * refinement[1],
            low[2] * refinement[2]
        },
        m_high
        {
            high[0] * refinement[0] - 1,
            high[1] * refinement[1] - 1,
            high[2] * refinement[2] - 1
        },
        m_periodic
        {
            periodic[0] * ( m_high[0] - m_low[0] + 1 ),
            periodic[1] * ( m_high[1] - m_low[1] + 1 ),
            periodic[2] * ( m_high[2] - m_low[2] + 1 )
        },
        m_prefinement_weight { 1. / ( refinement[0] * refinement[1] * refinement[2] ) }
    {
        m_box2patch.reserve ( nboxes );
    }

    virtual ~PartData()
    {
    }

    const int *
    iLower (
        const int & box
    ) const
    {
        return m_ilowers[box];
    }

    const int &
    iLowerD (
        const int & box,
        const int & d
    ) const
    {
        return m_ilowers[box][d];
    }

    const int *
    iUpper (
        const int & box
    ) const
    {
        return m_iuppers[box];
    }

    const int &
    iUpperD (
        const int & box,
        const int & d
    ) const
    {
        return m_iuppers[box][d];
    }

    const int &
    boxSize (
        const int & box
    ) const
    {
        return m_boxsizes[box];
    }

    const int &
    patch (
        const int & box
    ) const
    {
        return m_box2patch[box];
    }

    const int &
    box (
        const int & patch
    ) const
    {
        return m_patch2box.at ( patch );
    }

    const int *
    periodic (
    ) const
    {
        return m_periodic;
    }

    const double &
    pRefinementWeight (
    ) const
    {
        return m_prefinement_weight;
    }

    void
    addBox (
        const int & id,
        const int ilower[HYPRE_MAXDIM],
        const int iupper[HYPRE_MAXDIM],
        const bool interface[2 * HYPRE_MAXDIM]
    )
    {
        int box = m_box2patch.size();
        m_boxsizes[box] = 1;
        for ( int i = 0; i < HYPRE_MAXDIM; i++ )
        {
            m_ilowers[box][i] = ilower[i];
            m_iuppers[box][i] = iupper[i] - 1;
            m_boxsizes[box] *= ( iupper[i] - ilower[i] );
            m_interfaces[box][i] = interface[i];
            m_interfaces[box][i + HYPRE_MAXDIM] = interface[i + HYPRE_MAXDIM];
        }

        // std::map::operator[] leaks memory!
        if ( m_patch2box.find ( id ) == m_patch2box.end() )
            m_patch2box.insert ( std::make_pair ( id, m_box2patch.size() ) );
        m_box2patch.emplace_back ( id );
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartData_h



