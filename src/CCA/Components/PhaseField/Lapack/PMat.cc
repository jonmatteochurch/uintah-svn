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
 * @file CCA/Components/PhaseField/Lapack/PMat.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#include <CCA/Components/PhaseField/Lapack/PMat.h>

using namespace Uintah;
using namespace PhaseField;
using namespace Lapack;

#ifdef HAVE_LAPACK

#include <Core/Malloc/Allocator.h>

#include <cstring>
#include <numeric>

PMat::PMat ( int m ) :
    raw ( scinew int[m] ),
    m ( m ),
    n ( m )
{
    std::iota ( raw, raw + m, 1 );
}

PMat::PMat ( const PMat & copy ) :
    raw ( scinew int[copy.m] ),
    m ( copy.m ),
    n ( copy.m )
{
    std::memcpy ( raw, copy.raw, m * sizeof ( double ) );
}

// PMat::PMat ( std::initializer_list<int> l ) :
//     raw ( ( int * ) std::malloc ( l.size() * sizeof ( int ) ) ),
//     m ( l.size() ),
//     n ( l.size() )
// {
//     std::copy ( l.begin(), l.end(), raw );
// }

PMat::~PMat()
{
    delete[] raw;
}

int * PMat::data() const
{
    return raw;
}

int PMat::operator() ( int i, int j )
{
    return ( i != j ) * raw[i];
}

#endif
