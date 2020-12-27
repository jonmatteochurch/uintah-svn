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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraIndices_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraIndices_h

#include <Core/Geometry/IntVector.h>
// 
#include <tuple>

namespace Uintah
{
namespace HypreSStruct
{

class AddIndex
    : public std::tuple<int, IntVector, int, IntVector>
{
public:
    using std::tuple<int, IntVector, int, IntVector>::tuple;

    const int &
    box() const
    {
        return std::get<0> ( *this );
    }

    const IntVector &
    index() const
    {
        return std::get<1> ( *this );
    }

    const int &
    toPart() const
    {
        return std::get<2> ( *this );
    }

    const IntVector &
    toIndex() const
    {
        return std::get<3> ( *this );
    }
};

class StnIndex
    : public std::tuple<int, IntVector, int>
{
public:
    using std::tuple<int, IntVector, int>::tuple;

    const int &
    box() const
    {
        return std::get<0> ( *this );
    }

    const IntVector &
    index() const
    {
        return std::get<1> ( *this );
    }

    const int &
    rank() const
    {
        return std::get<2> ( *this );
    }
};

inline std::ostream & operator<< ( std::ostream & os, const AddIndex & add )
{
    os << "{" << add.box() << "," << add.index() << "," << add.toPart() << "," << add.toIndex() << "}";
    return os;
}

inline std::ostream & operator<< ( std::ostream & os, const StnIndex & add )
{
    os << "{" << add.box() << "," << add.index() << "," << add.rank() << "}";
    return os;
}

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_ExtraIndices_h



