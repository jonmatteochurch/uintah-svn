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

#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_stencil.h>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

template<> int sstruct_stencil<2>::offsets[5][2] = {
    {0,0},
    {-1,0},
    {1,0},
    {0,-1},
    {0,1}
};

template<> const int sstruct_stencil<2>::entry[5] = { 6, 0, 1, 2, 3 };

template<> int sstruct_stencil<3>::offsets[7][3] = {
    {0,0,0},
    {-1,0,0},
    {1,0,0},
    {0,-1,0},
    {0,1,0},
    {0,0,-1},
    {0,0,1},
};

template<> const int sstruct_stencil<3>::entry[7] = { 6, 0, 1, 2, 3, 4, 5 };

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah
