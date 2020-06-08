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
 * @file CCA/Components/PhaseField/Lapack/DMat.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#include <CCA/Components/PhaseField/Lapack/DMat.h>
#include <Core/Util/Assert.h>

#ifdef HAVE_LAPACK

#include <sci_defs/uintah_defs.h> // For FIX_NAME

#define DGEQP3 FIX_NAME(dgeqp3)
#define DGESVD FIX_NAME(dgesvd)

extern "C" void DGEQP3 ( const int * m, const int * n, double * a, const int * lda, int * jpvt, const double * tau, double * work, const int * lwork, int * info );
extern "C" void DGESVD ( const char * jobu, const char * jobvt, const int * m, const int * n, double * a, const int * lda, double * s, double * u, const int * ldu, double * vt, const int * ldvy, double * work, const int * lwork, int * info );

#include <Core/Malloc/Allocator.h>

#include <algorithm>
#include <cstring>

using namespace Uintah;
using namespace PhaseField;
using namespace Lapack;

DMat::DMat ( int m, int n ) :
    raw ( scinew double[m * n] ),
    pp ( nullptr ),
    ptau ( nullptr ),
    ps ( nullptr ),
    pu ( nullptr ),
    pvt ( nullptr ),
    m ( m ),
    n ( n ),
    k ( std::min ( m, n ) )
{}

DMat::DMat ( int m, int n, double v ) :
    raw ( scinew double[m * n] ),
    pp ( nullptr ),
    ptau ( nullptr ),
    ps ( nullptr ),
    pu ( nullptr ),
    pvt ( nullptr ),
    m ( m ),
    n ( n ),
    k ( std::min ( m, n ) )
{
    std::fill ( raw, raw + m * n, v );
}

DMat::DMat ( const DMat & copy ) :
    raw ( scinew double[copy.m * copy.n] ),
    pp ( copy.pp ? scinew PMat ( *copy.pp ) : nullptr ),
    ptau ( copy.ptau ? scinew DVec ( *copy.ptau ) : nullptr ),
    ps ( copy.ps ? scinew DVec ( *copy.ps ) : nullptr ),
    pu ( copy.pu ? scinew DMat ( *copy.pu ) : nullptr ),
    pvt ( copy.pvt ? scinew DMat ( *copy.pvt ) : nullptr ),
    m ( copy.m ),
    n ( copy.n ),
    k ( copy.k )
{
}

DMat::~DMat()
{
    delete pp;
    delete ptau;
    delete ps;
    delete pu;
    delete pvt;
    delete[] raw;
}

double &
DMat::operator() ( int i, int j )
{
    return raw[j * m + i];
}

const double &
DMat::operator() ( int i, int j ) const
{
    return raw[j * m + i];
}

double *
DMat::data () const
{
    return raw;
}

double *
DMat::cdata ( int j ) const
{
    return raw + j * m;
}

const PMat &
DMat::p() const
{
    ASSERTMSG ( pp, "accessing qr decomposition before computing it" );
    return *pp;
}

const DVec &
DMat::tau () const
{
    ASSERTMSG ( ptau, "accessing qr decomposition before computing it" );
    return *ptau;
}

const DVec &
DMat::s() const
{
    ASSERTMSG ( ps, "accessing svd decomposition before computing it" );
    return *ps;
}

const DMat &
DMat::u() const
{
    ASSERTMSG ( pu, "accessing svd decomposition before computing it" );
    return *pu;
}

const DMat &
DMat::vt() const
{
    ASSERTMSG ( pvt, "accessing svd decomposition before computing it" );
    return *pvt;
}

void
DMat::qr ( bool preserve )
{
    if ( !pp ) pp = scinew PMat ( n );
    if ( !ptau ) ptau = scinew DVec ( k );

    double * tmp;
    if ( preserve )
    {
        tmp = scinew double[m * n];
        std::memcpy ( tmp, raw, m * n * sizeof ( double ) );
    }

    double * work = ( double * ) std::malloc ( sizeof ( double ) );
    int lwork = -1;
    int info;

    // retrieve optimal work size
    DGEQP3 ( &m, &n, raw, &m, pp->data(), ptau->data(), work, &lwork, &info );

    lwork = work[0];
    work = ( double * ) std::realloc ( work, lwork * sizeof ( double ) );

    DGEQP3 ( &m, &n, raw, &m, pp->data(), ptau->data(), work, &lwork, &info );

    std::free ( work );

    if ( preserve )
    {
        delete[] raw;
        raw = tmp;
    }
}

void
DMat::svd ( bool preserve )
{
    if ( !ps ) ps = scinew DVec ( k );
    if ( !pu ) pu = scinew DMat ( m, k );
    if ( !pvt ) pvt = scinew DMat ( k, n );

    double * tmp;
    if ( preserve )
    {
        tmp = scinew double[m * n];
        std::memcpy ( tmp, raw, m * n * sizeof ( double ) );
    }

    double * work = ( double * ) std::malloc ( sizeof ( double ) );
    int lwork = -1;
    int info;

    // retrieve optimal work size
    DGESVD ( "S", "S", &m, &n, raw, &m, ps->data(), pu->data(), &m, pvt->data(), &k, work, &lwork, &info );

    lwork = work[0];
    work = ( double * ) std::realloc ( work, lwork * sizeof ( double ) );

    DGESVD ( "S", "S", &m, &n, raw, &m, ps->data(), pu->data(), &m, pvt->data(), &k, work, &lwork, &info );

    std::free ( work );

    if ( preserve )
    {
        delete[] raw;
        raw = tmp;
    }
}

DMat
DMat::Vander ( int m, int n, const double * c )
{
    DMat v ( m, n );
    DVec d ( m, 1 );

    for ( int j = n - 1; j >= 0; --j )
    {
        std::memcpy ( v.cdata ( j ), d.data(), m * sizeof ( double ) );
        for ( int i = 0; i < m; ++i )
            d ( i ) *= c[i];
    }
    return v;
}

#endif
