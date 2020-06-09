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
 * @file CCA/Components/PhaseField/Lapack/Tanh1.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/04
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Lapack_Tanh1_h
#define Packages_Uintah_CCA_Components_PhaseField_Lapack_Tanh1_h

#include <sci_defs/lapack_defs.h>

#ifdef HAVE_LAPACK

#include <CCA/Components/PhaseField/Lapack/DMat.h>

#include <string>
#include <iomanip>
#include <numeric>

#define _p16 std::scientific << std::setprecision(17) << std::setw(19) << std::right <<
#define _w15 std::scientific << std::setprecision(4) << std::setw(15) << std::right <<
#define _e4 std::scientific << std::setprecision(4) <<

namespace Uintah
{
namespace PhaseField
{
namespace Lapack
{

struct TrustRegionSetup
{
    double ftol, xtol, gtol, trtol;
    int max_nfev, max_triter;
};


template < int N >
int solve_lsq_trust_region ( const Lapack::DMat & J, const double Uf[N], const double SUf[N], const double SUf2[N], const double & Delta2, double ( &d ) [N], double & alpha, const double & max_triter, const double & trtol2 )
{
    const DVec & S = J.s();
    const DMat & VT = J.vt();

    int iter = 0;
    double tmp0, tmp1, pnorm, phi, phi_prime, ratio;

    double Delta = std::sqrt ( Delta2 );
    double alpha_upper = std::sqrt ( std::accumulate( SUf2, SUf2+N, 0. ) / Delta2 );
    double alpha_lower = 0.;

    // Check if J has full rank and try Gauss-Newton step.
    bool full_rank = ( J.m >= N ) && ( S ( N - 1 ) > DBL_EPSILON * J.m * S ( 0 ) );

    if ( full_rank )
    {
        std::fill ( d, d + N, 0. );
        for ( int i = 0; i < N; ++i )
            for ( int j = 0; j < N; ++j )
                d[j] -= VT ( i, j ) * Uf[i] / S ( i );

        double step_norm2 = 0.;
        for ( int j = 0; j < N; ++j )
            step_norm2 += d[j] * d[j];

        if ( step_norm2 <= Delta2 ) return iter;

        pnorm = phi_prime = 0.;
        for ( int i = 0; i < N; ++i )
        {
            tmp0 = S ( i ) * S ( i );
            tmp1 = SUf2[i] / ( tmp0 * tmp0 );
            pnorm += tmp1;
            phi_prime -= tmp1 / tmp0;
        }
        pnorm = std::sqrt ( pnorm );
        phi = pnorm - Delta;
        phi_prime /= pnorm;
        alpha_lower = -phi / phi_prime;
    }
    else
    {
        alpha = std::max ( 0.001 * alpha_upper, std::sqrt ( alpha_lower * alpha_upper ) );
    }

    while ( ++iter <= max_triter )
    {
        if ( alpha < alpha_lower || alpha > alpha_upper )
            alpha = std::max ( 0.001 * alpha_upper, std::sqrt ( alpha_lower * alpha_upper ) );

        pnorm = phi_prime = 0.;
        for ( int i = 0; i < N; ++i )
        {
            tmp0 = S ( i ) * S ( i ) + alpha;
            tmp1 = SUf2[i] / ( tmp0 * tmp0 );
            pnorm += tmp1;
            phi_prime -= tmp1 / tmp0;
        }
        pnorm = std::sqrt ( pnorm );
        phi = pnorm - Delta;
        phi_prime /= pnorm;
        ratio = phi / phi_prime;

        if ( phi < 0 ) alpha_upper = alpha;
        if ( alpha - ratio  > alpha_lower ) alpha_lower = alpha - ratio;

        alpha = alpha - ( phi + Delta ) * ratio / Delta;
        if ( phi * phi < trtol2 * Delta2 ) break;
    }

    std::fill ( d, d + N, 0. );
    for ( int i = 0; i < N; ++i )
    {
        tmp0 = S ( i ) * S ( i ) + alpha;
        for ( int j = 0; j < N; ++j )
            d[j] -= VT ( i, j ) * SUf[i] / tmp0;
    }

    double step_norm2 = 0.;
    for ( int j = 0; j < N; ++j )
        step_norm2 += d[j] * d[j];

    /* Make the norm of p equal to Delta, p is changed only slightly
     * during this. It is done to prevent p lie outside the trust
     * region (which can cause problems later).
     */
    double adj = std::sqrt ( Delta2 / step_norm2 );

    for ( int j = 0; j < N; ++j )
        d[j] *= adj;

    return iter;
}

class Tanh1
{
    static constexpr unsigned NOT_CONVERGED = 0x00, GTOL_SATISFIED = 0x01, S2TOL_SATISFIED = 0x02, BETATOL_SATISFIED = 0x04;
    static constexpr int A = 0, B = 1, C = 2;
    static constexpr int N = 3;

    static constexpr bool dbg_fit = false;
    static constexpr int verbosity = 2;

    static std::string
    termination_message ( unsigned status )
    {
        switch ( status )
        {
        case NOT_CONVERGED:
            return "The maximum number of function evaluations is exceeded.";
        case GTOL_SATISFIED:
            return "`gtol` termination condition is satisfied.";
        case S2TOL_SATISFIED:
            return "`ftol` termination condition is satisfied.";
        case BETATOL_SATISFIED:
            return  "`xtol` termination condition is satisfied.";
        case GTOL_SATISFIED | S2TOL_SATISFIED:
            return  "Both `gtol` and `ftol` termination conditions are satisfied.";
        case GTOL_SATISFIED | BETATOL_SATISFIED:
            return  "Both `gtol` and `xtol` termination conditions are satisfied.";
        case S2TOL_SATISFIED | BETATOL_SATISFIED:
            return  "Both `ftol` and `xtol` termination conditions are satisfied.";
        case GTOL_SATISFIED | S2TOL_SATISFIED | BETATOL_SATISFIED:
            return  "All `gtol`, `ftol` and `xtol` termination conditions are satisfied.";
        default:
            std::string res ( 55, ' ' );
            std::sprintf ( &res[0], "Improper status (%#04x) returned from `leastsq`", ( unsigned ) status );
            return res;
        }
    }

    double m_S2tol, m_betatol2, m_gtol, m_trtol2;
    int m_max_nfev, m_max_triter;
    double m_a0, m_a, m_b, m_c;

public:
    Tanh1 (
        const double & gamma0,
        const double & r0,
        const TrustRegionSetup & setup = { 1e-8, 1e-8, 1e-8, 1e-2, 200, 10 }
    ) : m_S2tol ( setup.ftol ),
        m_betatol2 ( setup.xtol * setup.xtol ),
        m_gtol ( setup.gtol ),
        m_trtol2 ( setup.trtol * setup.trtol ),
        m_max_nfev ( setup.max_nfev ),
        m_max_triter ( setup.max_triter ),
        m_a0 ( 0 ),
        m_a ( -r0 * r0 ),
        m_b ( 0. ),
        m_c ( gamma0 )
    {}

    double operator() ( const double & x ) const
    {
        return -std::tanh ( m_a + m_b * x + m_c * x * x  );
    }

    double zero () const
    {
        return ( std::sqrt ( m_b * m_b - 4. * m_a * m_c ) - m_b ) / ( 2. * m_c );
    }

    bool fit ( int m, double * x, double * z, double a0 )
    {
        m_a += a0 - m_a0;
        m_a0 = a0;

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
            for ( int i = 0; i < m - 1; ++i ) ss << _p16 z[i] << ",";
            ss << _p16 z[m - 1] << "];";
            DOUT ( true, ss.str() );

            DOUT ( true, "beta=[" << _p16 m_a << "," << _p16 m_b << "," << _p16 m_c << "];" );
            DOUT ( true, "% beta = lsq1(x,y,z,beta);" );
        }

        double * zh = scinew double[m];
        double * dz = scinew double[m];
        double * zh_old = scinew double[m];
        double * dz_old = scinew double[m];

        Lapack::DMat J ( m, N );

        int iter = 0;
        int nfev = 0;
        int njev = 0;

        double S2, S20, g[N], tmp;

        S2 = 0.;
        std::fill ( g, g + N, 0. );
        for ( int i = 0; i < m; ++i )
        {
            zh[i] = ( *this ) ( x[i] );
            dz[i] = zh[i] - z[i];

            tmp = zh[i] * zh[i] - 1.;
            J ( i, A ) = tmp;
            J ( i, B ) = x[i] * tmp;
            J ( i, C ) = x[i] * x[i] * tmp;

            for ( int j = 0; j < N; ++j )
                g[j] += J ( i, j ) * dz[i];

            S2 += dz[i] * dz[i];
        }
        ++nfev;
        ++njev;

        S20 = S2;

        double Delta2 = m_a * m_a + m_b * m_b + m_c * m_c;
        if ( Delta2 == 0. ) Delta2 = 1.;
        double alpha = 0.;
        unsigned converged = NOT_CONVERGED;

        DOUT ( verbosity > 1, ">> " << _w15 "Iteration" << _w15 "Total nfev"
               << _w15 "S2" << _w15 "S2 reduction" << _w15 "Step norm"
               << _w15 "Optimality" );

        double g_norm = std::abs ( g[0] );
        for ( int j = 1; j < N; ++j )
            if ( g_norm < ( tmp = std::abs ( g[j] ) ) )
                g_norm = tmp;

        if ( g_norm < m_gtol ) converged = GTOL_SATISFIED;

        DOUT ( verbosity > 1, ">> " << _w15 iter << _w15 nfev << _w15 S2
               << _w15 "" << _w15 "" << _w15 g_norm );

        double udz[N], sudz[N], sudz2[N], q, S2_new, dS2, dSh2, a_old, b_old, c_old, d[N], step_norm2;
        while ( !converged && nfev < m_max_nfev )
        {
            J.svd ( true );
            const DMat & U = J.u();
            const DVec & S = J.s();

            std::fill ( udz, udz + N, 0. );
            for ( int i = 0; i < m; ++i )
                for ( int j = 0; j < N; ++j )
                    udz[j] += U ( i, j ) * dz[i];

            if ( m >= N )
                for ( int j = 0; j < N; ++j )
                    sudz[j] = S ( j ) * udz[j];
            else
                for ( int j = 0; j < N; ++j )
                    sudz[j] = S ( 0 ) * udz[j];

            for ( int j = 0; j < N; ++j )
                sudz2[j] = sudz[j] * sudz[j];

            dS2 = -1;
            a_old = m_a;
            b_old = m_b;
            c_old = m_c;
            std::swap ( zh_old, zh );
            std::swap ( dz_old, dz );

            while ( !converged && dS2 <= 0 && nfev < m_max_nfev )
            {
                solve_lsq_trust_region ( J, udz, sudz, sudz2, Delta2, d, alpha, m_max_triter, m_trtol2 );

                m_a = a_old + d[A];
                m_b = b_old + d[B];
                m_c = c_old + d[C];

                q = S2_new = 0.;
                for ( int i = 0; i < m; ++i )
                {
                    tmp = 0.;
                    for ( int j = 0; j < N; ++j )
                        tmp += J ( i, j ) * d[j];
                    q += tmp * tmp;

                    zh[i] = ( *this ) ( x[i] );
                    dz[i] = zh[i] - z[i];

                    S2_new += dz[i] * dz[i];
                }
                ++nfev;

                step_norm2 = 0.;
                dSh2 = -q;
                for ( int j = 0; j < N; ++j )
                {
                    step_norm2 += d[j] * d[j];
                    dSh2 -=  2. * g[j] * d[j];
                }
                dS2 = S2 - S2_new;

                bool bound_hit = step_norm2 > 0.9025 * Delta2;

                double ratio = 0.;
                if ( dSh2 > 0 ) ratio = dS2 / dSh2;
                else if ( dSh2 == 0 && dS2 == 0 ) ratio = 1.;

                if ( ratio < 0.25 )
                {
                    alpha *= 4. * std::sqrt ( Delta2 / step_norm2 );
                    Delta2 = step_norm2 / 16.;
                }
                else if ( ratio > 0.75 && bound_hit )
                {
                    alpha *= 0.5;
                    Delta2 *= 4.;
                }

                if ( dS2 < m_S2tol * S2 && ratio > 0.25 )
                    converged |= S2TOL_SATISFIED;

                if ( step_norm2 < m_betatol2 * ( m_a * m_a + m_b * m_b + m_c * m_c ) )
                    converged |= BETATOL_SATISFIED;
            }

            if ( dS2 > 0 )
            {
                S2 = S2_new;

                std::fill ( g, g + N, 0. );
                for ( int i = 0; i < m; ++i )
                {
                    tmp = zh[i] * zh[i] - 1.;
                    J ( i, A ) = tmp;
                    J ( i, B ) = x[i] * tmp;
                    J ( i, C ) = x[i] * x[i] * tmp;

                    for ( int j = 0; j < N; ++j )
                        g[j] += J ( i, j ) * dz[i];
                }
                ++njev;
            }
            else
            {
                m_a = a_old;
                m_b = b_old;
                m_c = c_old;
                std::swap ( zh, zh_old );
                std::swap ( dz, dz_old );
                step_norm2 = 0;
                dS2 = 0;
            }

            g_norm = std::abs ( g[0] );
            for ( int j = 1; j < N; ++j )
                if ( g_norm < ( tmp = std::abs ( g[j] ) ) )
                    g_norm = tmp;

            if ( g_norm < m_gtol )
                converged |= GTOL_SATISFIED;

            DOUT ( verbosity > 1, ">> " << _w15 ++iter << _w15 nfev << _w15 S2
                   << _w15 dS2 << _w15 std::sqrt ( step_norm2 )
                   << _w15 g_norm );
        }

        delete[] zh;
        delete[] dz;
        delete[] zh_old;
        delete[] dz_old;

        DOUT ( verbosity, ">> " << termination_message ( converged ).c_str()
               << "\n>> Function evaluations " << nfev
               << ", initial cost " << _e4 S20
               << ", final cost " << _e4 S2
               << ", first-order optimality " << _e4 g_norm  << "." );

        if ( dbg_fit )
        {
            std::stringstream ss;

            DOUT ( true, "beta=[" << _p16 m_a << "," << _p16 m_b << "," << _p16 m_c << "];" );
            DOUT ( true, ss.str() );
            DOUT ( true, "" );
        }

        return converged;
    }
};

} // namespace Lapack
} // namespace PhaseField
} // namespace Uintah

#undef _w15
#undef _p16
#undef _e4

#endif

#endif // Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostprocessorFactory_h







