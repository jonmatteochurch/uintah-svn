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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraValues_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraValues_h

#include <CCA/Components/Solvers/HypreSStruct/ExtraIndices.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructStencil.h>
#include <Core/Grid/Variables/CCVariable.h>

namespace Uintah
{
namespace HypreSStruct
{

class ExtraValue
{
public:
    virtual ~ExtraValue() = default;

    virtual double
    value (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries
    ) const = 0;
};

template<int DIM>
class StencilEntryValue:
    public ExtraValue
{
    int box;
    IntVector index;
    int entry;
    double weight;

public:
    StencilEntryValue (
        const StnIndex & index,
        const double & weight
    ) : box ( index.box() ),
        index ( index.index() ),
        entry ( SStructStencil<DIM>::entry[index.rank()] ),
        weight ( weight )
    {
    }

    virtual double
    value (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** /*additional_entries*/
    ) const override
    {
        return weight * stencil_entries[box][index][entry];
    }
};

class AdditionalEntryValue:
    public ExtraValue
{
    int box;
    MatrixIndex index;
    double weight;

public:
    AdditionalEntryValue (
        const AddIndex & index,
        const double & weight
    ) : box ( index.box() ),
        index ( index.index(), index.toPart(), index.toIndex() ),
        weight ( weight )
    {}

    virtual
    double
    value (
        constCCVariable<Stencil7> * /*stencil_entries*/,
        AdditionalEntries ** additional_entries
    ) const override
    {
        return weight * additional_entries[box]->at ( index );
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraValues_h
