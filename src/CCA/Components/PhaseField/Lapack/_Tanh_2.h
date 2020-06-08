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
 * @file CCA/Components/PhaseField/Lapack/Tanh.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/04
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Lapack_Tanh_h
#define Packages_Uintah_CCA_Components_PhaseField_Lapack_Tanh_h

#include <sci_defs/lapack_defs.h>

#ifdef HAVE_LAPACK

#include <CCA/Components/PhaseField/Lapack/DMat.h>

#include <string>
#include <iomanip>

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

class Tanh1
{
    static constexpr unsigned NOT_CONVERGED = 0x00, GTOL_SATISFIED = 0x01, S2TOL_SATISFIED = 0x02, BETATOL_SATISFIED = 0x04;

    static constexpr bool dbg_fit = false;
    static constexpr int verbosity = 0;

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
    double m_r2, m_gamma/*, m_b*/;

    int solve_lsq_trust_region ( const Lapack::DMat & J, const double & Uf, const double & SUf, const double & SUf2, const double & Delta2, double & dr2, double & alpha )
    {
        int iter = 0;
        double d, tmp, pnorm, phi, phi_prime, ratio;

        double Delta = std::sqrt ( Delta2 );
        double alpha_upper = std::sqrt ( SUf2 / Delta2 );
        double alpha_lower = 0.;

        alpha = std::max ( 0.001 * alpha_upper, std::sqrt ( alpha_lower * alpha_upper ) );

        while ( ++iter <= m_max_triter )
        {
            if ( alpha < alpha_lower || alpha > alpha_upper )
                alpha = std::max ( 0.001 * alpha_upper, std::sqrt ( alpha_lower * alpha_upper ) );

            d = J.s() ( 0 ) * J.s() ( 0 ) + alpha;
            tmp = SUf2 / ( d * d );
            pnorm = std::sqrt ( tmp );
            phi = pnorm - Delta;
            phi_prime = - ( tmp / d ) / pnorm;
            ratio = phi / phi_prime;

            if ( phi < 0 ) alpha_upper = alpha;
            if ( alpha - ratio  > alpha_lower ) alpha_lower = alpha - ratio;

            alpha = alpha - ( phi + Delta ) * ratio / Delta;
            if ( phi * phi < m_trtol2 * Delta2 ) break;
        }

        d = J.s() ( 0 ) * J.s() ( 0 ) + alpha;
        dr2 = - J.vt() ( 0, 0 ) * SUf / d;

        /* Make the norm of p equal to Delta, p is changed only slightly
         * during this. It is done to prevent p lie outside the trust
         * region (which can % cause problems later).
         */
        double adj = std::sqrt ( Delta2 / ( dr2 * dr2 ) );
        dr2 *= adj;

        return iter;
    }

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
        m_r2 ( r0 * r0 ),
        m_gamma ( gamma0 )//,
//         m_b ( m_a * m_r2 )
    {}

//     void reset ( const double & gamma0, const double & r0 )
//     {
//         m_r2 = r0 * r0;
//         m_gamma = gamma0;
// //         m_b = m_a * m_r2;
//     }

//     void y0 ( const double & y0 )
//     {
// //         m_b = m_a * ( m_r2 - y0 * y0 );
//     }
    double r0() const
    {
        return std::sqrt(m_r2);
    }

    double operator() ( const double & x, const double & y ) const
    {
        return -std::tanh ( m_gamma * ( x * x + y * y - m_r2 ) );
    }

    double zero( const double & y ) const
    {
        return std::sqrt ( m_r2 - y * y );
    }

    double dx ( const double & x, const double & y ) const
    {
        double tmp = std::tanh ( m_gamma * ( x * x + y * y - m_r2 ) );
        return -2. * m_gamma * x * ( 1. - tmp * tmp );
    }

    bool fit ( int m, double * x, double * y, double * z )
    {
        if ( dbg_fit )
        {
            std::stringstream ss;
            ss << "PY:x   =np.array([";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 x[i] << ",";
            ss << _p16 x[m-1] << "]);";
            DOUT ( true, ss.str() );

            ss.str ( "" );
            ss.clear();

            ss << "PY:y   =np.array([";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 y[i] << ",";
            ss << _p16 y[m-1] << "]);";
            DOUT ( true, ss.str() );

            ss.str ( "" );
            ss.clear();

            ss << "PY:z   =np.array([";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 z[i] << ",";
            ss << _p16 z[m-1] << "]);";
            DOUT ( true, ss.str() );

            DOUT ( true, "PY:beta=np.array([" << _p16 m_r2 << "]);" );
            DOUT ( true, "PY:beta = least_squares( fun, beta, jac=jac, verbose=2, x_scale=beta_scale, ftol=ftol, xtol=xtol, gtol=gtol );" );
            DOUT ( true, "PY:");

            ss.str("");
            ss.clear();

            ss << "ML:x   =[";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 x[i] << ",";
            ss << _p16 x[m-1] << "];";
            DOUT ( true, ss.str() );

            ss.str ( "" );
            ss.clear();

            ss << "ML:y   =[";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 y[i] << ",";
            ss << _p16 y[m-1] << "];";
            DOUT ( true, ss.str() );

            ss.str ( "" );
            ss.clear();

            ss << "ML:z   =[";
            for ( int i = 0; i < m-1; ++i ) ss << _p16 z[i] << ",";
            ss << _p16 z[m-1] << "];";
            DOUT ( true, ss.str() );

            DOUT ( true, "ML:beta=[" << _p16 m_r2 << "];" );
            DOUT ( true, "ML:beta_scale=[ " << _p16 m_r2 << "," << _p16 1 << " ];" );
            DOUT ( true, "ML:beta = least_squares( x', y', beta', 'verbose', 2, 'beta_scale', beta_scale, 'ftol', ftol, 'xtol', xtol, 'gtol', gtol );" );
            DOUT ( true, "ML:");
        }

        double * zh = scinew double[m];
        double * dz = scinew double[m];
        double * zh_old = scinew double[m];
        double * dz_old = scinew double[m];

        Lapack::DMat J ( m, 1 );
        double * j = J.cdata ( 0 );

        int iter = 0;
        int nfev = 0;
        int njev = 0;

        double S2, S20, g;
        S2 = g = 0.;

        for ( int i = 0; i < m; ++i )
        {
            zh[i] = ( *this ) ( x[i], y[i] );
            dz[i] = zh[i] - z[i];
//
            j[i] = 2. * m_gamma * ( 1. - zh[i] * zh[i] );

            g += j[i] * dz[i];
            S2 += dz[i] * dz[i];
        }
        ++nfev;
        ++njev;

        S20 = S2;

        double Delta2 = m_r2 * m_r2;
        if ( Delta2 == 0. ) Delta2 = 1.;
        double alpha = 0.;
        unsigned converged = NOT_CONVERGED;

        DOUT ( verbosity > 1, ">> " << _w15 "Iteration" << _w15 "Total nfev"
               << _w15 "S2" << _w15 "S2 reduction" << _w15 "Step norm"
               << _w15 "Optimality" );

        double g_norm = std::abs ( g );
        if ( g_norm < m_gtol ) converged = GTOL_SATISFIED;

        DOUT ( verbosity > 1, ">> " << _w15 iter << _w15 nfev << _w15 S2
               << _w15 "" << _w15 "" << _w15 g_norm );

        double udz, sudz, sudz2, q, S2_new, dS2, dSh2, r2_old, dr2, step_norm2, tmp;
        while ( !converged && nfev < m_max_nfev )
        {
            J.svd ( true );
            j = J.cdata ( 0 );

            udz = 0.;
            for ( int i = 0; i < m; ++i )
                udz += J.u() ( i, 0 ) * dz[i];

            sudz = J.s() ( 0 ) * udz;
            sudz2 = sudz * sudz;

            dS2 = -1;
            r2_old = m_r2;
            std::swap ( zh_old, zh );
            std::swap ( dz_old, dz );

            while ( !converged && dS2 <= 0 && nfev < m_max_nfev )
            {
                solve_lsq_trust_region ( J, udz, sudz, sudz2, Delta2, dr2, alpha );

                m_r2 = r2_old + dr2;

                q = 0.;
                S2_new = 0;
                for ( int i = 0; i < m; ++i )
                {
                    tmp = J ( i, 0 ) * dr2;
                    q += tmp * tmp;

                    zh[i] = ( *this ) ( x[i], y[i] );
                    dz[i] = zh[i] - z[i];

                    S2_new += dz[i] * dz[i];
                }
                ++nfev;

                step_norm2 = dr2 * dr2;
                dSh2 = - q - 2. * ( g * dr2);
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

                step_norm2 = dr2 * dr2;
                if ( step_norm2 < m_betatol2 * ( m_r2 * m_r2 ) )
                    converged |= BETATOL_SATISFIED;
            }

            if ( dS2 > 0 )
            {
                S2 = S2_new;

                g = 0.;
                for ( int i = 0; i < m; ++i )
                {
                    j[i] = 2. * m_gamma * ( 1. - zh[i] * zh[i] );

                    g += j[i] * dz[i];
                }
                ++njev;
            }
            else
            {
                m_r2 = r2_old;
                std::swap ( zh, zh_old );
                std::swap ( dz, dz_old );
                step_norm2 = 0;
                dS2 = 0;
            }

            g_norm = std::abs ( g );
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







