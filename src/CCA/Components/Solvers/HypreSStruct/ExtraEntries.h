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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraEntries_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraEntries_h

#include <CCA/Components/Solvers/HypreSStruct/ExtraValues.h>

#include <vector>

namespace Uintah
{
namespace HypreSStruct
{

class ExtraEntry
{
protected:
    const IntVector * m_index;
    std::vector<int> m_entries;
    std::vector<double> m_values;

public :
    int * index() const
    {
        return const_cast<IntVector *> ( m_index )->get_pointer();
    }

    int * entries() const
    {
        return const_cast<int *> ( m_entries.data() );
    }

    double * values() const
    {
        return const_cast<double *> ( m_values.data() );
    }

    int nValues() const
    {
        return m_values.size();
    }
};

template<typename Index, int DIM>
class ExtraIterator
    : protected ExtraEntry
{
    using map = std::map< Index, std::list<ExtraValue *> >;
    using cit = typename map::const_iterator;

    constCCVariable<Stencil7> * m_stencil_entries;
    AdditionalEntries ** m_additional_entries;

    const cit m_end;
    cit m_it;
    cit m_next;

    void
    update (
        const StnIndex & stn
    )
    {
        if ( m_end == m_it ) return;

        this->m_index = &std::get<1> ( stn );
        this->m_entries.clear();
        this->m_values.clear();

        for ( m_next = m_it; m_next != m_end; ++m_next )
        {
            if ( std::get<1> ( m_next->first ) != *this->m_index )
                break;

            double value = 0.;
            for ( const auto & extra : m_next->second )
                value += extra->value ( m_stencil_entries, m_additional_entries );

            this->m_entries.emplace_back ( std::get<2> ( m_next->first ) );
            this->m_values.emplace_back ( value );
        }
    }

    void
    update (
        const AddIndex & add
    )
    {
        if ( m_end == m_it ) return;

        this->m_index = &std::get<1> ( add );
        this->m_entries.clear();
        this->m_values.clear();

        int entry = SStructStencil<DIM>::size;
        for ( m_next = m_it; m_next != m_end; ++m_next )
        {
            if ( std::get<1> ( m_next->first ) != *this->m_index )
                break;

            double value = 0.;
            for ( const auto & extra : m_next->second )
                value += extra->value ( m_stencil_entries, m_additional_entries );

            this->m_entries.emplace_back ( entry++ );
            this->m_values.emplace_back ( value );
        }
    }

public:
    ExtraIterator
    (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries,
        const int & capacity,
        const cit & end,
        const cit & it
    ) : m_stencil_entries ( stencil_entries ),
        m_additional_entries ( additional_entries ),
        m_end ( end ),
        m_it ( it )
    {
        this->m_index = nullptr;
        this->m_entries.reserve ( capacity );
        this->m_values.reserve ( capacity );
        update ( m_it->first );
    }

    bool
    operator != (
        const ExtraIterator & it
    ) const
    {
        return it.m_it != m_it;
    };

    ExtraIterator &
    operator ++ ()
    {
        m_it = m_next;
        update ( m_it->first );
        return *this;
    };

    ExtraEntry
    operator * ()
    {
        return *this;
    }
};

template<typename Index, int DIM>
class ExtraRange
{
    using map = std::map< Index, std::list<ExtraValue *> >;
    using cit = typename map::const_iterator;

    constCCVariable<Stencil7> * m_stencil_entries;
    AdditionalEntries ** m_additional_entries;
    const cit m_begin;
    const cit m_end;
    const int m_capacity;

public:
    ExtraRange (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries,
        const map * extra_entries
    ) : m_stencil_entries ( stencil_entries ),
        m_additional_entries ( additional_entries ),
        m_begin ( extra_entries->begin() ),
        m_end ( extra_entries->end() ),
        m_capacity ( extra_entries->size() )
    {}

    ExtraIterator<Index, DIM>
    begin() const
    {
        return ExtraIterator<Index, DIM> ( m_stencil_entries, m_additional_entries, m_capacity, m_end, m_begin );
    }

    ExtraIterator<Index, DIM>
    end() const
    {
        return ExtraIterator<Index, DIM> ( m_stencil_entries, m_additional_entries, 0, m_end, m_end );
    }
};

template<typename Index, int DIM>
class ExtraEntries
{
    std::map< Index, std::list<ExtraValue *> > m_extra_entries;

public:
    ~ExtraEntries()
    {
        for ( const std::pair< Index, std::list<ExtraValue *> > & entry: m_extra_entries )
            for ( ExtraValue * value: entry.second )
                delete value;
    }

    bool
    emplace_back (
        const Index & ind,
        const StnIndex & stn,
        const double & weight
    )
    {
        ExtraValue * value = scinew StencilEntryValue<DIM> { stn, weight };
        auto res = m_extra_entries.emplace ( ind, std::list<ExtraValue *> ( 1, value ) );
        if ( !res.second ) // entry already in map
            res.first->second.emplace_back ( value );
        return res.second;
    }

    bool
    emplace_back (
        const Index & ind,
        const AddIndex & add,
        const double & weight
    )
    {
        ExtraValue * value = scinew AdditionalEntryValue { add, weight };
        auto res = m_extra_entries.emplace ( ind, std::list<ExtraValue *> ( 1, value ) );
        if ( !res.second ) // entry already in map
            res.first->second.emplace_back ( value );
        return res.second;
    }

    ExtraRange<Index, DIM>
    operator() (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries
    )
    {
        return ExtraRange<Index, DIM> ( stencil_entries, additional_entries, &m_extra_entries );
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraEntries_h

