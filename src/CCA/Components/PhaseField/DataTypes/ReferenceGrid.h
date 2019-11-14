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

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_ReferenceGrid_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_ReferenceGrid_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <Core/Grid/Grid.h>

namespace Uintah
{
namespace PhaseField
{

template<DimType DIM>
class ReferenceGrid
    : public Grid
{
public:
    ReferenceGrid (
        const Grid * grid,
        const int & index
    ) : Grid()
    {
        d_levels.reserve ( index + 2 );
        for ( int i = 0; i <= index; ++i )
            d_levels.emplace_back ( grid->getLevel ( i ) );
    }

    virtual ~ReferenceGrid () = default;

    void
    mapPatchToFinest (
        const Patch * patch,
        const int & index,
        IntVector & low,
        IntVector & high
    ) const
    {

        IntVector r_ratio  = IntVector ( 1, 1, 1 );
        for ( size_t i = index + 1; i < d_levels.size(); ++i )
            r_ratio  = r_ratio * d_levels[i]->getRefinementRatio();

        low = patch->getCellLowIndex() * r_ratio;
        high  = patch->getCellHighIndex() * r_ratio;


        for ( size_t D = 0; D < DIM; ++D )
            if ( r_ratio[D] > 1 )
            {
                if ( low[D] < 0 ) low[D] += 1;
                if ( high[D] < 0 ) high[D] += 1;
            }
    }

    Level *
    addFinestLevel ( 
        int k
    )
    {
        Point anchor ( d_levels.back()->getAnchor() );
        Vector dcell ( d_levels.back()->dCell() );

        int r = 1 << k;
        for ( size_t D = 0; D < DIM; ++D )
            dcell[D] /= r;
        return addLevel ( anchor, dcell );
    }

    Patch *
    addFinestPatch (
        const Patch * patch,
        const int & index 
    )
    {
        IntVector low, high;
        mapPatchToFinest ( patch, index, low, high );
        return d_levels.back()->addPatch ( low, high, low, high, this, -999 );
    }
};

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_Entry_h


