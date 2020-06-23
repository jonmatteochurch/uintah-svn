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
 * @file CCA/Components/PhaseField/Lapack/Poly.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#include <CCA/Components/PhaseField/Lapack/Poly.h>

#ifdef HAVE_LAPACK

#include <sci_defs/uintah_defs.h> // For FIX_NAME

#define DGEEV FIX_NAME(dgeev)

extern "C" void DGEEV ( const char * jobvl, const char * jobvr, const int * n, double * a, const int * lda, double * wr, double * wi, double * vl, int * ldvl, double * vr, int * ldvr, double * work, const int * lwork, int * info );

#include <CCA/Components/PhaseField/Lapack/DMat.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>

#include <cstring>
#include <cmath>
#include <cstdlib>

#include <iomanip>
#define _p16 std::scientific << std::setprecision(17) << std::setw(19) << std::right <<

using namespace Uintah;
using namespace PhaseField;
using namespace Lapack;

static constexpr bool dbg_fit = false;

Poly::Poly ( int d ) :
    raw ( scinew double[d + 1] ),
    rr ( scinew double[d] ),
    ri ( scinew double[d] ),
    d ( d ),
    nr ( -1 )
{
}

Poly::Poly ( const Poly & copy ) :
    raw ( scinew double[copy.d + 1] ),
    rr ( scinew double[copy.d] ),
    ri ( scinew double[copy.d] ),
    d ( copy.d ),
    nr ( copy.nr )
{
    std::memcpy ( raw, copy.raw, ( d + 1 ) * sizeof ( double ) );
    std::memcpy ( raw, copy.raw, d * sizeof ( double ) );
    std::memcpy ( raw, copy.raw, d * sizeof ( double ) );
}

Poly::~Poly()
{
    delete[] raw;
    delete[] rr;
    delete[] ri;
}

double &
Poly::cfx ( int i )
{
    return raw[d - i];
}

double
Poly::cfx ( int i ) const
{
    return raw[i];
}

bool
Poly::root_in_range ( double min, double max, double & r ) const
{
    for ( int i = 0; i < nr; ++i )
        if ( ri[i] == 0 && min <= rr[i] && rr[i] <= max )
        {
            r = rr[i];
            return true;
        }

    // if almost tangent to y=0 return 0
    if ( abs ( raw[d] ) < 1e-3 )
    {
        r = 0.;
        return min <= r && r <= max;
    }

    DOUT ( dbg_fit, "range: [" << min << "," << max << "]" );
    for ( int i = 0; i < nr; ++i )
        DOUT ( dbg_fit, "root i: " << rr[i] << " + " << ri[i] << "i" );

    return false;
}

double
Poly::operator() ( double x0 ) const
{
    double b = raw[0];
    for ( int i = 1; i <= d; ++i )
        b = raw[i] + b * x0;
    return b;
}

double
Poly::dx ( double x0 ) const
{
    double b = d * raw[0];
    for ( int i = 1; i <= d - 1; ++i )
        b = ( d - i ) * raw[i] + b * x0;
    return b;
}

double
Poly::ddx ( double x0 ) const
{
    double b = d * ( d - 1 ) * raw[0];
    for ( int i = 1; i <= d - 2; ++i )
        b = ( d - i ) * ( d - i - 1 ) * raw[i] + b * x0;
    return b;
}

void
Poly::fit ( int m, const double * x, const double * y )
{
    if ( dbg_fit )
    {
        std::stringstream ss;

        ss << "x   =[";
        for ( int i = 0; i < m - 1; ++i ) ss << _p16 x[i] << ",";
        ss << _p16 x[m - 1] << "];";
        DOUT ( true, ss.str() );

        ss.str ( "" );
        ss.clear();

        ss << "y   =[";
        for ( int i = 0; i < m - 1; ++i ) ss << _p16 y[i] << ",";
        ss << _p16 y[m - 1] << "];";
        DOUT ( true, ss.str() );

        DOUT ( true, "p = polyfit(x, y, " << d << ")" );
    }

    // Construct the Vandermonde matrix.
    DMat a = DMat::Vander ( m, d + 1, x );
    DVec b = DVec ( m, y );

    // Solve by QR decomposition.
    a.qr();
    b.premult_trsp_q ( a );
    b.premult_inv_r ( a );
    b.premult_inv_perm ( a.p() );
    std::memcpy ( raw, b.data(), ( d + 1 ) * sizeof ( double ) );

    if ( dbg_fit )
    {
        std::stringstream ss;

        ss << "p   =[";
        for ( int i = 0; i < d; ++i ) ss << _p16 raw[i] << ",";
        ss << _p16 raw[d] << "];";
        DOUT ( true, ss.str() );
        DOUT ( true, "" );
    }
}

void
Poly::roots()
{
    nr = 0;

    double v_max = 0., tmp;
    for ( int i = 0; i <= d; ++i )
        if ( ( tmp = std::abs ( raw[i] ) ) > v_max )
            v_max = tmp;

    if ( v_max == 0 ) return;

    double * v = scinew double[d + 1];
    int * f = scinew int[d + 1];

    int k = -1; // leading zero cfx
    while ( raw[++k] == 0 );

    int l = d + 1; // trailing zero cfx
    while ( raw[--l] == 0 );

    int m = 0;
    for ( int i = k; i <= l; ++i )
    {
        v[m] = raw[i] / v_max;
        f[m] = i;
        ++m;
    }

    if ( m > 1 )
    {
        nr = m - 1;
        DMat a ( nr, nr, 0 );
        for ( int i = 1; i < nr; ++i )
            a ( i, i - 1 ) = 1;
        for ( int i = 0; i < nr; ++i )
            a ( 0, i ) = -v[i + 1] / v[0];

        double * work = ( double * ) std::malloc ( sizeof ( double ) );
        int lwork = -1;
        int ldv = 1;
        int info;

        DGEEV ( "N", "N", &a.m, a.data(), &a.m, rr, ri, nullptr, &ldv, nullptr, &ldv, work, &lwork, &info );

        lwork = work[0];
        work = ( double * ) std::realloc ( work, lwork * sizeof ( double ) );

        DGEEV ( "N", "N", &a.m, a.data(), &a.m, rr, ri, nullptr, &ldv, nullptr, &ldv, work, &lwork, &info );

        for ( int i = f[m - 1]; i < d; ++i )
        {
            rr[i] = 0.;
            ri[i] = 0.;
            ++nr;
        }

        std::free ( work );
    }
    else
    {
        nr = d - f[m - 1];
        for ( int i = 0; i < nr; ++i )
        {
            rr[i] = 0.;
            ri[i] = 0.;
        }
    }

    delete[] f;
    delete[] v;
}

#endif
