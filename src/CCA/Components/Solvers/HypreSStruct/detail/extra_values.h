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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_extra_values_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_extra_values_h

#include <Core/Grid/Variables/CCVariable.h>

#include <CCA/Components/Solvers/HypreSStruct/detail/extra_indices.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_stencil.h>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

class extra_value
{
public:
    virtual ~extra_value() = default;

    virtual double
    value (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries
    ) const = 0;

    virtual std::ostream &
    print (
        std::ostream & os
    ) const = 0;
};

inline std::ostream & operator<< ( std::ostream & os, const extra_value & value )
{
    return value.print(os);
}

template<int DIM>
class stencil_entry_value:
    public extra_value
{
    int m_box;
    IntVector m_index;
    int m_entry;
    double m_weight;

public:
    stencil_entry_value (
        const stn_index & index,
        const double & weight
    ) : m_box ( index.box() ),
        m_index ( index.index() ),
        m_entry ( sstruct_stencil<DIM>::entry[index.rank()] ),
        m_weight ( weight )
    {
    }

    virtual double
    value (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** /*additional_entries*/
    ) const override
    {
        return m_weight * stencil_entries[m_box][m_index][m_entry];
    }

    virtual std::ostream &
    print (
        std::ostream & os
    ) const override
    {
         return os << "{stencil," << m_box << "," << m_index << "," << m_entry << "," << m_weight << "}";
    }
};

class additional_entry_value:
    public extra_value
{
    int m_box;
    MatrixIndex m_index;
    double m_weight;

public:
    additional_entry_value (
        const add_index & index,
        const double & weight
    ) : m_box ( index.box() ),
        m_index ( index.index(), index.to_part(), index.to_index() ),
        m_weight ( weight )
    {}

    virtual
    double
    value (
        constCCVariable<Stencil7> * /*stencil_entries*/,
        AdditionalEntries ** additional_entries
    ) const override
    {
        return m_weight * additional_entries[m_box]->at ( m_index );
    }

    virtual std::ostream &
    print (
        std::ostream & os
    ) const override
    {
         return os << "{additional," << m_box << "," << m_index << "," << m_weight << "}";
    }
};

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_extra_values_h
