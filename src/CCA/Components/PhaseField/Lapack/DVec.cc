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
 * @file CCA/Components/PhaseField/Lapack/DVec.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#include <CCA/Components/PhaseField/Lapack/DVec.h>

#ifdef HAVE_LAPACK

#include <sci_defs/uintah_defs.h> // For FIX_NAME

#define DORMQR FIX_NAME(dormqr)
#define DTRTRS FIX_NAME(dtrtrs)
#define DLAPMR FIX_NAME(dlapmr)

extern "C" void DORMQR ( const char * side, const char * trans, const int * m, const int * n, const int * k, const double * a, const int * lda, const double * tau, double * c, const int * ldc, double * work, int * lwork, int * info );
extern "C" void DTRTRS ( const char * uplo, const char * trans, const char * diag, const int * n, const int * nrhs, const double * a, const int * lda, double * b, const int * ldb, int * info );
extern "C" void DLAPMR ( const int * forward, const int * m, const int * n, double * x, const int * ldx, int * k );

#include <Core/Malloc/Allocator.h>
#include <Core/Util/Assert.h>
#include <CCA/Components/PhaseField/Lapack/DMat.h>

#include <algorithm>
#include <cstring>

using namespace Uintah;
using namespace PhaseField;
using namespace Lapack;

DVec::DVec ( int m )
    : raw ( ( double * ) std::malloc ( m * sizeof ( double ) ) ),
      m ( m ),
      n ( 1 )
{
}

DVec::DVec ( int m, double v )
    : raw ( ( double * ) std::malloc ( m * sizeof ( double ) ) ),
      m ( m ),
      n ( 1 )
{
    std::fill ( raw, raw + m, v );
}

DVec::DVec ( int m, const double * v )
    : raw ( ( double * ) std::malloc ( m * sizeof ( double ) ) ),
      m ( m ),
      n ( 1 )
{
    std::memcpy ( raw, v, m * sizeof ( double ) );
}

// DVec::DVec ( std::initializer_list<double> l )
//     : raw ( ( double * ) std::malloc ( l.size() * sizeof ( double ) ) ),
//       m ( l.size() ),
//       n ( 1 )
// {
//     std::copy ( l.begin(), l.end(), raw );
// }

DVec::DVec ( const DVec & copy )
    : raw ( ( double * ) std::malloc ( copy.m * sizeof ( double ) ) ),
      m ( copy.m ),
      n ( copy.n )
{
    std::memcpy ( raw, copy.raw, m * sizeof ( double ) );
}

DVec::~DVec()
{
    std::free ( raw );
}

double *
DVec::data() const
{
    return raw;
}

double &
DVec::operator() ( int i )
{
    return raw[i];
}

const double & 
DVec::operator() ( int i ) const
{
    return raw[i];
}

void
DVec::premult_trsp_q ( const DMat & qrp )
{
    ASSERT ( m == qrp.m );

    double * work = ( double * ) std::malloc ( sizeof ( double ) );
    int lwork = -1;
    int info;

    dormqr_ ( "L", "T", &m, &n, &qrp.k, qrp.data(), &qrp.m, qrp.tau().data(), raw, &m, work, &lwork, &info );

    lwork = work[0];
    work = ( double * ) std::realloc ( work, lwork * sizeof ( double ) );

    dormqr_ ( "L", "T", &m, &n, &qrp.k, qrp.data(), &qrp.m, qrp.tau().data(), raw, &m, work, &lwork, &info );

    std::free ( work );

    m = qrp.k;
    raw = ( double * ) std::realloc ( raw, m * sizeof ( double ) );
}

void
DVec::premult_inv_r ( const DMat & qrp )
{
    ASSERT ( m == qrp.k );

    int info;
    dtrtrs_ ( "U", "N", "N", &qrp.k, &n, qrp.data(), &qrp.m, raw, &m, &info );
}

void
DVec::premult_inv_perm ( const PMat & p )
{
    ASSERT ( m == p.m );

    int forward = 0;
    dlapmr_ ( &forward, &m, &n, raw, &m, p.data() );
}

#endif
