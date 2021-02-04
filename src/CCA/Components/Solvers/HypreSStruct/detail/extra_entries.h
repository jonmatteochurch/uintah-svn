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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_extra_entries_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_extra_entries_h

#include <CCA/Components/Solvers/HypreSStruct/detail/extra_values.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_stencil.h>

#include <vector>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

enum extra_operation { NO_OP, ADDTO, SET };

class extra_entry
{
protected:
    const IntVector * m_index;
    std::vector<int> m_set_entries;
    std::vector<double> m_set_values;
    std::vector<int> m_addto_entries;
    std::vector<double> m_addto_values;

public :
    bool get_set_values (
        int *& index,
        int & nentries,
        int *& entries,
        double *& values
    ) const
    {
        index = const_cast<IntVector *> ( m_index )->get_pointer();
        nentries = m_set_entries.size();
        entries = const_cast<int *> ( m_set_entries.data() );
        values = const_cast<double *> ( m_set_values.data() );
        return nentries;
    }

    bool get_addto_values (
        int *& index,
        int & nentries,
        int *& entries,
        double *& values
    ) const
    {
        index = const_cast<IntVector *> ( m_index )->get_pointer();
        nentries = m_addto_entries.size();
        entries = const_cast<int *> ( m_addto_entries.data() );
        values = const_cast<double *> ( m_addto_values.data() );
        return nentries;
    }
};

template<typename Index, int DIM>
class extra_iterator
    : protected extra_entry
{
    using extra_values = std::list<extra_value *>;
    using extra_pair = std::pair<extra_operation, extra_values>;
    using pair_map = std::map<Index, extra_pair>;
    using pair_cit = typename pair_map::const_iterator;

    constCCVariable<Stencil7> * m_stencil_entries;
    AdditionalEntries ** m_additional_entries;

    const pair_cit m_end;
    pair_cit m_it;
    pair_cit m_next;

    void
    update (
        const stn_index & stn
    )
    {
        if ( m_end == m_it ) return;

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\textra_iterator update stn " << stn << ": ";
#endif

        this->m_index = &std::get<1> ( stn );
        this->m_set_entries.clear();
        this->m_set_values.clear();
        this->m_addto_entries.clear();
        this->m_addto_values.clear();

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\nm_it\t" << m_it->first <<  " -> " << ( m_it->second.first == ADDTO ? "ADDTO " : ( m_it->second.first == SET ? "SET " : "NO_OP " ) );
        for ( const auto & pe : m_it->second.second ) std::cout << *pe << " ";
#endif

        for ( m_next = m_it; m_next != m_end; ++m_next )
        {
#ifdef EXTRA_ENTRIES_DBG
            std::cout << "\nm_next\t" << m_next->first <<  " -> " << ( m_next->second.first == ADDTO ? "ADDTO " : ( m_next->second.first == SET ? "SET " : "NO_OP " ) );
            for ( const auto & pe : m_next->second.second ) std::cout << *pe << " ";
#endif

            if ( std::get<1> ( m_next->first ) != *this->m_index )
                break;

            const extra_pair & pair ( m_next->second );
            const extra_operation & operation ( pair.first );
            const extra_values & values ( pair.second );

            double value = 0.;
            for ( const auto & extra : values )
                value += extra->value ( m_stencil_entries, m_additional_entries );

            switch ( operation )
            {
            case SET:
                ASSERTEQ ( values.size(), 1 );
                this->m_set_entries.emplace_back ( std::get<2> ( m_next->first ) );
                this->m_set_values.emplace_back ( value );
                break;
            case ADDTO:
                this->m_addto_entries.emplace_back ( std::get<2> ( m_next->first ) );
                this->m_addto_values.emplace_back ( value );
                break;
            case NO_OP:
            default:
                ASSERTFAIL ( "operation should be set" );
            }
        }

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\n";
#endif

        ASSERTMSG ( std::none_of ( m_set_entries.begin(), m_set_entries.end(), [this] ( int e )
        {
            return std::count ( m_addto_entries.begin(), m_addto_entries.end(), e );
        } ),
        "only one operation is permitted for each entry" );

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "set_entries=[ ";
        std::copy ( this->m_set_entries.begin(), this->m_set_entries.end(), std::ostream_iterator<int> ( std::cout, " " ) );
        std::cout << "] set_values=[ ";
        std::copy ( this->m_set_values.begin(), this->m_set_values.end(), std::ostream_iterator<double> ( std::cout, " " ) );
        std::cout << "] addto_entries=[ ";
        std::copy ( this->m_addto_entries.begin(), this->m_addto_entries.end(), std::ostream_iterator<int> ( std::cout, " " ) );
        std::cout << "] addto_values=[ ";
        std::copy ( this->m_addto_values.begin(), this->m_addto_values.end(), std::ostream_iterator<double > ( std::cout, " " ) );
        std::cout << "]" << std::endl;
#endif
    }

    void
    update (
        const add_index & add
    )
    {
        if ( m_end == m_it ) return;

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\textra_iterator update add " << add << ": ";
#endif

        this->m_index = &std::get<1> ( add );
        this->m_set_entries.clear();
        this->m_set_values.clear();
        this->m_addto_entries.clear();
        this->m_addto_values.clear();

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\nm_it\t" << m_it->first <<  " -> " << ( m_it->second.first == ADDTO ? "ADDTO " : ( m_it->second.first == SET ? "SET " : "NO_OP " ) );
        for ( const auto & pe : m_it->second.second ) std::cout << *pe << " ";
#endif

        for ( m_next = m_it; m_next != m_end; ++m_next )
        {
#ifdef EXTRA_ENTRIES_DBG
            std::cout << "\nm_next\t" << m_next->first <<  " -> " << ( m_next->second.first == ADDTO ? "ADDTO " : ( m_next->second.first == SET ? "SET " : "NO_OP " ) );
            for ( const auto & pe : m_next->second.second ) std::cout << *pe << " ";
#endif

            if ( std::get<1> ( m_next->first ) != *this->m_index )
                break;

            const extra_pair & pair ( m_next->second );
            const extra_operation & operation ( pair.first );
            const extra_values & values ( pair.second );

            double value = 0.;
            for ( const auto & extra : values )
                value += extra->value ( m_stencil_entries, m_additional_entries );

            switch ( operation )
            {
            case SET:
                ASSERTEQ ( values.size(), 1 );
                this->m_set_entries.emplace_back ( std::get<2> ( m_next->first ) );
                this->m_set_values.emplace_back ( value );
                break;
            case ADDTO:
                this->m_addto_entries.emplace_back ( std::get<2> ( m_next->first ) );
                this->m_addto_values.emplace_back ( value );
                break;
            case NO_OP:
            default:
                ASSERTFAIL ( "operation should be set" );
            }
        }
#ifdef EXTRA_ENTRIES_DBG
        std::cout << "\n";
#endif

        ASSERTMSG ( std::none_of ( m_set_entries.begin(), m_set_entries.end(), [this] ( int e )
        {
            return std::count ( m_addto_entries.begin(), m_addto_entries.end(), e );
        } ),
        "only one operation is permitted for each entry" );

#ifdef EXTRA_ENTRIES_DBG
        std::cout << "set_entries=[ ";
        std::copy ( this->m_set_entries.begin(), this->m_set_entries.end(), std::ostream_iterator<int> ( std::cout, " " ) );
        std::cout << "] set_values=[ ";
        std::copy ( this->m_set_values.begin(), this->m_set_values.end(), std::ostream_iterator<double> ( std::cout, " " ) );
        std::cout << "] addto_entries=[ ";
        std::copy ( this->m_addto_entries.begin(), this->m_addto_entries.end(), std::ostream_iterator<int> ( std::cout, " " ) );
        std::cout << "] addto_values=[ ";
        std::copy ( this->m_addto_values.begin(), this->m_addto_values.end(), std::ostream_iterator<double> ( std::cout, " " ) );
        std::cout << "]" << std::endl;
#endif
    }

public:
    extra_iterator
    (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries,
        const int & capacity,
        const pair_cit & end,
        const pair_cit & it
    ) : m_stencil_entries ( stencil_entries ),
        m_additional_entries ( additional_entries ),
        m_end ( end ),
        m_it ( it )
    {
        this->m_index = nullptr;
        this->m_set_entries.reserve ( capacity );
        this->m_set_values.reserve ( capacity );
        this->m_addto_entries.reserve ( capacity );
        this->m_addto_values.reserve ( capacity );
        update ( m_it->first );
    }

    bool
    operator != (
        const extra_iterator & it
    ) const
    {
        return it.m_it != m_it;
    };

    extra_iterator &
    operator ++ ()
    {
        m_it = m_next;
        update ( m_it->first );
        return *this;
    };

    extra_entry &
    operator * ()
    {
        return *this;
    }

    const extra_entry &
    operator * ()
    const
    {
        return *this;
    }
};

template<typename Index, int DIM>
class extra_range
{
    using extra_values = std::list<extra_value *>;
    using extra_pair = std::pair<extra_operation, extra_values>;
    using pair_map = std::map<Index, extra_pair>;
    using pair_cit = typename pair_map::const_iterator;

    constCCVariable<Stencil7> * m_stencil_entries;
    AdditionalEntries ** m_additional_entries;
    const pair_cit m_begin;
    const pair_cit m_end;
    const int m_capacity;

public:
    extra_range (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries,
        const pair_map * extra_entries
    ) : m_stencil_entries ( stencil_entries ),
        m_additional_entries ( additional_entries ),
        m_begin ( extra_entries->begin() ),
        m_end ( extra_entries->end() ),
        m_capacity ( extra_entries->size() )
    {}

    extra_iterator<Index, DIM>
    begin() const
    {
        return extra_iterator<Index, DIM> ( m_stencil_entries, m_additional_entries, m_capacity, m_end, m_begin );
    }

    extra_iterator<Index, DIM>
    end() const
    {
        return extra_iterator<Index, DIM> ( m_stencil_entries, m_additional_entries, 0, m_end, m_end );
    }
};

template<typename Index, int DIM>
class extra_entries
{
    using extra_values = std::list<extra_value *>;
    using extra_pair = std::pair<extra_operation, extra_values>;
    using pair_map = std::map<Index, extra_pair>;
    using pair_it = typename pair_map::iterator;

    pair_map m_extra_entries;

    template<typename I>
    typename std::enable_if< std::is_same<typename std::remove_const<I>::type, stn_index>::value, bool>::type
    insert (
        const extra_operation & op,
        const I & ind,
        extra_value * value
    )
    {
        bool inserted;
        pair_it it;
        std::tie ( it, inserted ) = m_extra_entries.emplace ( ind, std::make_pair ( op, extra_values ( 1, value ) ) );

        if ( !inserted ) // entry already in map
        {
            extra_pair & pair = it->second;
            ASSERTEQ ( pair.first, ADDTO );
            ASSERTEQ ( op, ADDTO );
            pair.second.emplace_back ( value );
        }
        return inserted;
    }

    template<typename I>
    typename std::enable_if< std::is_same<typename std::remove_const<I>::type, add_index>::value, bool>::type
    insert (
        const extra_operation & op,
        I & ind,
        extra_value * value
    )
    {
        std::get<2> ( ind ) = -1;
        pair_it it = m_extra_entries.lower_bound ( ind );

        std::get<2> ( ind ) = sstruct_stencil<DIM>::size;
        for ( ; it->first.index() == ind.index(); ++it, ++std::get<2> ( ind ) )
        {
            if ( it->first == ind )
            {
                extra_pair & pair = it->second;
                ASSERTEQ ( pair.first, ADDTO );
                ASSERTEQ ( op, ADDTO );
                pair.second.emplace_back ( value );
                return false;
            }
        }
        m_extra_entries.emplace_hint ( it, ind, std::make_pair ( op, extra_values ( 1, value ) ) );
        return true;
    }

public:
    ~extra_entries()
    {
        for ( const auto & it : m_extra_entries )
            for ( extra_value * value : it.second.second )
                delete value;
    }

    bool
    emplace_back (
        const extra_operation & op,
        typename std::conditional<std::is_same<Index, add_index>::value, Index &, const Index &>::type ind,
        const stn_index & stn,
        const double & weight
    )
    {
#ifdef EXTRA_ENTRIES_DBG
        std::cout << "extra_entries::emplace_back: " << ( op == ADDTO ? "ADDTO " : ( op == SET ? "SET " : "NO_OP " ) ) << "at " << ind <<  " from " << weight << " x stn " << stn << std::endl;
#endif
        extra_value * value = scinew stencil_entry_value<DIM> { stn, weight };
        return insert ( op, ind, value );
    }

    bool
    emplace_back (
        const extra_operation & op,
        typename std::conditional<std::is_same<Index, add_index>::value, Index &, const Index &>::type ind,
        const add_index & add,
        const double & weight
    )
    {
#ifdef EXTRA_ENTRIES_DBG
        std::cout << "extra_entries::emplace_back: " << ( op == ADDTO ? "ADDTO " : ( op == SET ? "SET " : "NO_OP " ) ) << "at " << ind <<  " from " << weight << " x add " << add << std::endl;
#endif
        extra_value * value = scinew additional_entry_value { add, weight };
        return insert ( op, ind, value );
    }

    inline void
    clear()
    {
#ifdef EXTRA_ENTRIES_DBG
        std::cout << "extra_entries::clear" << std::endl;
#endif
        m_extra_entries.clear();
    }

    extra_range<Index, DIM>
    operator() (
        constCCVariable<Stencil7> * stencil_entries,
        AdditionalEntries ** additional_entries
    )
    {
        return extra_range<Index, DIM> ( stencil_entries, additional_entries, &m_extra_entries );
    }

    virtual std::ostream &
    print (
        std::ostream & os
    ) const
    {
        os << "{ extra entries:";
        for ( const auto & it : m_extra_entries )
        {
            os << "\n\t" << it.first <<  " -> " << ( it.second.first == ADDTO ? "ADDTO " : ( it.second.first == SET ? "SET " : "NO_OP " ) );
            for ( const auto & pe : it.second.second ) os << *pe << " ";
        }
        return os << "\n}";
    }
};

template<typename Index, int DIM>
inline std::ostream & operator<< ( std::ostream & os, const extra_entries<Index, DIM> & entries )
{
    return entries.print ( os );
}

} // namespace detailextra_operation
} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_extra_entries_h
