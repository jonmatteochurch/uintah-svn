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

/**
 * @file CCA/Components/PhaseField/DataTypes/Entries.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2019/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataType_Entries_h
#define Packages_Uintah_CCA_Components_PhaseField_DataType_Entries_h

#include <CCA/Components/PhaseField/DataTypes/Entry.h>

namespace Uintah
{
namespace PhaseField
{

/**
 * @brief System row entries
 *
 * Structure containing a list of matrix entries and rhs value on a
 * (not given) row of a linear system
 */
template <typename T>
struct Entries
{
public: // MEMBERS

    /// Matrix row entries
    std::list<Entry<T> > values;

    /// Rhs entry
    T rhs;

public: // CONSTRUCTORS

    /**
     * @brief Constructor
     *
     * Populates members
     *
     * @tparam E... matrix entry type
     * @param rhs rhs entry
     * @param values matrix row entries
     */
    template<typename... E>
    Entries (
        T rhs,
        E... values
    ) : values { values... },
        rhs ( rhs )
    {}

    /**
     * @brief Constructor
     *
     * Populates members (homogeneous rhs)
     *
     * @tparam E... matrix entry type
     * @param values matrix row entries
     */
    template<typename... E>
    Entries (
        E... values
    ) : values { values... },
        rhs ( 0. )
    {}

    /// Default move constructor
    Entries ( Entries<T> && copy ) = default;

public: // MODIFIERS

    /**
     * @brief Matrix row entries sum
     *
     * Add to the system row another row
     *
     * @param entries system row entries to be added
     */
    void
    add (
        Entries<T> && entries
    )
    {
        values.splice ( values.end(), entries.values );
        rhs += entries.rhs;
    }

    /**
     * @brief Matrix row entries weighted sum
     *
     * Add to the system row another row multiplied by a weight
     *
     * @param entries system row entries to be added (rvalue)
     * @param weight field value multiplier
     */
    void
    add (
        Entries<T> && entries,
        T weight = 1.
    )
    {
        for ( Entry<T> & entry : entries.values )
            entry.weight *= weight;
        values.splice ( values.end(), entries.values );
        rhs += weight * entries.rhs;
    }

    /**
     * @brief Simplify row entries
     *
     * Check if multiple row entries points to the same matrix element
     * and combine them
     */
    void
    simplify()
    {
        for ( auto it1 = values.begin(); it1 != values.end(); it1++ )
        {
            auto it2 = std::next ( it1 );
            while ( it2 != values.end() )
            {
                if ( it1->level == it2->level &&
                     it1->index[0] == it2->index[0] &&
                     it1->index[1] == it2->index[1] )
                {
                    it1->weight += it2->weight;
                    it2 = values.erase ( it2 );
                }
                else
                    it2++;
            }
        }
    }
};

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataType_Entries_h


