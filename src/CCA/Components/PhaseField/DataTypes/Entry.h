/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
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
 * @file CCA/Components/PhaseField/DataTypes/Entry.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2019/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_Entry_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_Entry_h

#include <Core/Geometry/IntVector.h>

namespace Uintah
{
namespace PhaseField
{

/**
 * @brief Matrix row entry
 *
 * Structure containing all information for defining the matrix
 * entry on a (not given) row at the column identified by the
 * level/index key pair.
 */
template <typename T>
struct Entry
{
public: // MEMBERS

    /// Column level index
    int level;

    /// Column grid entry index
    IntVector index;

    /// Matrix entry value
    T weight;

public: // CONSTRUCTORS

    /**
     * @brief Constructor
     *
     * Populates members
     *
     * @param level column level index
     * @param index column grid entry index
     * @param weight field value multiplier
     */
    Entry (
        int level,
        IntVector index,
        T weight = 1
    ) : level ( level ),
        index ( index ),
        weight ( weight )
    {}

};

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_Entry_h


