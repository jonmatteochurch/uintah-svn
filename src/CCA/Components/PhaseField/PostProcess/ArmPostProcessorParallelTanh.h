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
 * @file CCA/Components/PhaseField/PostProcess/ArmPostProcessorParallelTanh.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/06
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorParallelTanh_h
#define Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorParallelTanh_h

#include <sci_defs/lapack_defs.h>

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessor.h>
#include <Core/Exceptions/TypeMismatchException.h>

#ifdef HAVE_LAPACK
#   include <CCA/Components/PhaseField/Lapack/Poly.h>
#   include <CCA/Components/PhaseField/Lapack/Tanh1.h>
#   include <CCA/Components/PhaseField/Lapack/Tanh2.h>
#endif

#include <numeric>
#include <functional>
#include <CCA/Components/PhaseField/Views/View.h>

#define DBG_PRINT 0

namespace Uintah
{

extern Dout g_mpi_dbg;

namespace PhaseField
{

template<VarType VAR, DimType DIM>
class ArmPostProcessorParallelTanh;

template<VarType VAR>
class ArmPostProcessorParallelTanh < VAR, D2 >
    : public ArmPostProcessor < VAR, D2 >
{
private: // STATIC MEMBERS

    static constexpr TypeDescription::Type var_td = VAR == CC ? TypeDescription::CCVariable : TypeDescription::NCVariable;
    static const IntVector en;
    static const IntVector et;

protected: // MEMBERS

    const bool m_dbg;

    const int m_n0; // no. of psi values for computation of 0-level, and psi_n

    const int m_n0l;
    const int m_n0h;

    const int m_nn;  // size of tip data
    const int m_nt;  // size of tip data

    const int m_nnl; // no of pts left of tip
    const int m_nnh; // no of pts right of tip

    Lapack::TrustRegionSetup m_psetup;

    IntVector m_origin;
    IntVector m_low;
    IntVector m_high;
    double m_dn;
    double m_dt;

    int m_locations_size;
    int m_location_n0;
    int * m_locations;

    int m_data_size;
    int m_data_n0;

    double * m_data_t; // [*,n0]
    double * m_data_z; // [*,n0]

    int m_tip_n;
    double * m_tip_z; // [nt,nn]

private: // MEMBERS

    int n_ind ( const IntVector & id ) const
    {
        return id[0] - m_origin[0];
    }

    int n_ind ( const int & x, const int & /*y*/ ) const
    {
        return x - m_origin[0];
    }

    int t_ind ( const IntVector & id ) const
    {
        return id[1] - m_origin[1];
    }

    int t_ind ( const int & /*x*/, const int & y ) const
    {
        return y - m_origin[1];
    }

    int x_ind ( const IntVector & id ) const
    {
        return id[0] + m_origin[0];
    }

    int x_ind ( const int & n, const int & /*t*/ ) const
    {
        return n + m_origin[0];
    }

    int y_ind ( const IntVector & id ) const
    {
        return id[1] + m_origin[1];
    }

    int y_ind ( const int & /*n*/, const int & t ) const
    {
        return t + m_origin[1];
    }

    template < VarType V = VAR, typename std::enable_if<V == CC, int>::type = 0 >
    double n_coord ( const int & n )
    {
        return ( n + .5 ) * m_dn;
    }

    template < VarType V = VAR, typename std::enable_if<V == NC, int>::type = 0 >
    double n_coord ( const int & n )
    {
        return n * m_dn;
    }

    template < VarType V = VAR, typename std::enable_if<V == CC, int>::type = 0 >
    double t_coord ( const int & t )
    {
        return ( t + .5 ) * m_dt;
    }

    template < VarType V = VAR, typename std::enable_if<V == NC, int>::type = 0 >
    double t_coord ( const int & t )
    {
        return t * m_dt;
    }

    inline double & data_t ( int i, int j )
    {
        return m_data_t[i * m_n0 + j];    // [*,n0]
    }
    inline double & data_z ( int i, int j )
    {
        return m_data_z[i * m_n0 + j];    // [*,n0]
    }
    inline double & tip_z ( int i, int j )
    {
        return m_tip_z[i * m_nn + j];    // [nt,nn]
    }

    inline const double & data_t ( int i, int j ) const
    {
        return m_data_t[i * m_n0 + j];    // [*,n0]
    }
    inline const double & data_z ( int i, int j ) const
    {
        return m_data_z[i * m_n0 + j];    // [*,n0]
    }
    inline const double & tip_z ( int i, int j ) const
    {
        return m_tip_z[i * m_nn + j];    // [nt,nn]
    }

    inline double * data_t ( int i )
    {
        return m_data_t + i * m_n0;    // [*,n0]
    }
    inline double * data_z ( int i )
    {
        return m_data_z + i * m_n0;    // [*,n0]
    }
    inline double * tip_z ( int i )
    {
        return m_tip_z + i * m_nn;  // [nn,nn]
    }

public: // CONSTRUCTORS/DESTRUCTOR

    ArmPostProcessorParallelTanh (
        const Lapack::TrustRegionSetup psetup,
        int n0,
        int nn,
        int nt,
        bool dbg
    ) : m_dbg ( dbg ),
        m_n0 ( n0 ),
        m_n0l ( ( m_n0 - 1 ) / 2 ),
        m_n0h ( m_n0 - m_n0l ),
        m_nn ( nn ),
        m_nt ( nt ),
        m_nnl ( ( m_nn - 1 ) / 2 ),
        m_nnh ( m_nn - m_nnl ),
        m_psetup ( psetup ),
        m_locations_size ( 0 ),
        m_locations ( nullptr ),
        m_data_t ( nullptr ),
        m_data_z ( nullptr ),
        m_tip_z ( nullptr )
    {}

    /**
    * @brief Destructor
    */
    virtual ~ArmPostProcessorParallelTanh ()
    {
        delete[] m_locations;
        delete[] m_data_t;
        delete[] m_data_z;
        delete[] m_tip_z;
    }

public:

    virtual void setLevel (
        const Level * level
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelTanh::setLevel: " << level->getIndex() << " " );

        m_origin = level->getCellIndex ( {0., 0., 0.} );
        level->computeVariableExtents ( var_td, m_low, m_high );
        m_dn = level->dCell() [0];
        m_dt = level->dCell() [1];
    }

    virtual void
    initializeLocations (
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelTanh::initializeLocations " );

        m_location_n0 = std::max ( 0, n_ind ( m_low ) );
        m_locations_size = n_ind ( m_high ) - m_location_n0;
        m_locations = scinew int[m_locations_size];
        std::fill ( m_locations, m_locations + m_locations_size, INT_MAX );
    };

    virtual void
    setLocations (
        IntVector const & low,
        IntVector const & high,
        View < ScalarField<const double> > const & psi
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelTanh::setLocations: " << low << high << " " );

        // first entry with given n
        IntVector start = Max ( low, m_origin );
        IntVector id {0, 0, 0};

        // move along increasing n
        for ( id[0] = start[0]; id[0] < high[0]; ++id[0] )
        {
            // move along increasing t
            int end1 = ( VAR == CC ) ? std::min ( high[1] - 1, id[0] ) : std::min ( high[1] - 1, id[0] + 1 );
            for ( id[1] = start[1]; id[1] < end1; ++id[1] )
                if ( psi[id]*psi[id + et] <= 0. )
                {
                    int ind = n_ind ( id ) - m_location_n0;
                    if ( t_ind ( id ) < m_locations[ind] )
                        m_locations[ind] = t_ind ( id );
                    break;
                }
        }
    }

    virtual void
    reduceLocations (
        const ProcessorGroup * myworld
    ) override
    {
        if ( myworld->nRanks() <= 1 ) return;

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessor::reduceMPI " );

        int error = Uintah::MPI::Allreduce ( MPI_IN_PLACE, m_locations, m_locations_size, MPI_INT, MPI_MIN, myworld->getComm() );

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessor::reduceMPI, done " );

        if ( error )
        {
            DOUT ( true, "ArmPostProcessor::reduceMPI: Uintah::MPI::Allreduce error: " << error );
            SCI_THROW ( InternalError ( "ArmPostProcessor::reduceMPI: MPI error", __FILE__, __LINE__ ) );
        }
    }

    virtual void
    printLocations (
        std::ostream & out
    ) const override
    {
        for ( int i = 0; i < m_locations_size; ++i )
            if ( m_locations[i] == INT_MAX )
                out << "___ ";
            else
                out << std::setw ( 3 ) << m_locations[i] << " ";
    }

    virtual void
    initializeData ()
    override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelPoly::initializeData " );

        int i = m_locations_size;
        for ( i = m_locations_size; i > 0 && m_locations[i - 1] == INT_MAX; --i );
        m_data_size = i--;

        for ( ; i > 1 && m_locations[i - 2] < INT_MAX && m_locations[i - 2] >= m_locations[i]; --i );

        m_data_n0 = i - 1;
        m_data_size -= m_data_n0;

        m_data_t = scinew double[m_n0 * m_data_size];
        m_data_z = scinew double[m_n0 * m_data_size];

        std::fill ( m_data_t, m_data_t + m_n0 * m_data_size, -DBL_MAX );
        std::fill ( m_data_z, m_data_z + m_n0 * m_data_size, -DBL_MAX );

        m_tip_n = -1;

        m_tip_z = scinew double[m_nn * m_nt];
        std::fill ( m_tip_z, m_tip_z + m_nn * m_nt, -1. );
    };

    virtual void
    setData (
        const IntVector & low,
        const IntVector & high,
        View < ScalarField<const double> > const & psi
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelPoly::setData: " << low << high << " " );

        // arm contour not here
        if ( n_ind ( high ) < m_data_n0 - m_nnl ) return;

        IntVector id, sym;
        int it0, it1;
        int in, it;

        int n_ = -1, t_; // tip candidates

        for ( in = std::max ( m_data_n0, n_ind ( low ) ); in < n_ind ( high ); ++in )
        {
            t_ = m_locations[ in - m_location_n0 ];
            if ( t_ < INT_MAX )
            {
                it = t_ - m_n0l;

                it0 = std::max ( it, t_ind ( low ) );
                it1 = std::min ( it + m_n0, t_ind ( high ) );

                id[0] = x_ind ( in, it );
                id[1] = y_ind ( in, it );
                id[2] = 0;

                n_ = in - m_data_n0;
                t_ = 0;

                // symmetry on y axis
                sym = { id[0], m_origin[1] - id[1], 0 };
                if ( VAR == CC ) --sym[1];

                for ( ; id[1] < m_origin[1]; ++it, ++id[1], --sym[1], ++t_ )
                    if ( low[1] <= sym[1] && sym[1] < high[1] )
                    {
                        data_t ( n_, t_ ) = t_coord ( it );
                        data_z ( n_, t_ ) = psi[sym];
                    }

                // skip non patch region
                for ( ; it < it0; ++it, ++id[1], ++t_ );

                // set patch data
                for ( ; it < it1; ++it, ++id[1], ++t_ )
                {
                    data_t ( n_, t_ ) = t_coord ( it );
                    data_z ( n_, t_ ) = psi[id];
                }

                // extend psi to -1 out computational boundary
                if ( high[1] == m_high[1] )
                {
                    for ( t_ = t_ind ( high ); t_ < m_n0; ++t_ )
                    {
                        data_t ( n_, t_ ) = data_t ( n_, t_ - 1 ) + m_dt;
                        data_z ( n_, t_ ) = -1.;
                    }
                }
            }
        }

        // check for tip
        if ( ( t_ind ( 0, high[1] - 1 ) + m_nt ) * ( t_ind ( low ) - m_nt ) <= 0 )
        {
            bool is_tip = false;
            if ( m_tip_n < 0 )
            {
                if ( n_ < 0 ) // move back from low[0] at most m_nnh steps
                {
                    in = n_ind ( low ) - m_location_n0;
                    int dn0 = std::max ( m_nnh, -in );
                    for ( int dn = - 1; dn >= dn0; --dn )
                    {
                        if ( m_locations[in  + dn] < INT_MAX )
                        {
                            n_ = in + dn - m_data_n0;
                            break;
                        }
                    }
                }

                if ( n_ < 0 ) // move forward from high[0] at most m_nnl + 1 steps
                {
                    in = n_ind ( high - 1 ) - m_location_n0;
                    int dn0 = std::min ( m_nnl + 1, m_locations_size - 1 - in );
                    for ( int dn = 1; dn <= dn0; ++dn )
                    {
                        if ( m_locations[in  + dn] < INT_MAX )
                        {
                            n_ = in + dn - m_data_n0;
                            break;
                        }
                    }
                }

                if ( n_ >= 0 ) // i may be tip position
                {
                    is_tip = true;
                    n_ += m_data_n0;
                    int imax = std::min ( n_ + m_nnh + 1, m_locations_size );
                    int i = n_ + 1;
                    for ( ; i < imax; ++i )
                        if ( m_locations[i] < INT_MAX )
                            n_ = i;
                    for ( ; i < m_locations_size; ++i )
                        if ( m_locations[i] < INT_MAX )
                        {
                            is_tip = false;
                            break;
                        }
                }
            }
            else
            {
                is_tip = ( n_ind ( low ) <= m_tip_n + m_nnh && m_tip_n - m_nnl - 1 < n_ind ( high ) );
                n_ = m_tip_n;
            }

            if ( is_tip )
            {
                ASSERT ( m_tip_n < 0 || m_tip_n == n_ );

                m_tip_n = n_;

                in = n_ - m_nnl;
                it = 0;

                id = { x_ind ( in, it ), y_ind ( in, it ), 0 };

                int in0 = std::max ( 0, n_ind ( low ) - in );
                int in1 = std::min ( m_nn, n_ind ( high ) - in );
                int it0 = std::max ( 0, low[1] - it );
                int it1 = std::min ( m_nt, high[1] - it );

                // simmetry on x axis
                if ( id[0] < m_origin[0] )
                {
                    sym = { m_origin[0] - id[0], 0, 0 };
                    if ( VAR == CC ) sym[0] -= 1;
                    for ( n_ = 0; n_ < in0; ++n_, ++id[0], --sym[0] )
                        if ( low[0] <= sym[0] && sym[0] < high[0] )
                        {
                            sym[1] = m_origin[1] + it0;
                            for ( t_ = it0; t_ < it1; ++t_, ++sym[1] )
                                tip_z ( t_, n_ ) = psi[sym];
                        }
                }
                else
                    id[0] += in0;

                for ( n_ = in0; n_ < in1; ++n_, ++id[0] )
                {
                    id[1] = m_origin[1] + it0;
                    for ( t_ = it0; t_ < it1; ++t_, ++id[1] )
                        tip_z ( t_, n_ ) = psi[id];
                }
            }
        }
    }

    virtual void
    reduceData (
        const ProcessorGroup * myworld
    ) override
    {
        if ( myworld->nRanks() <= 1 ) return;

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessor::reduceMPI " );

        int error;
        if ( myworld->myRank() == 0 )
            error = Uintah::MPI::Reduce ( MPI_IN_PLACE, &m_tip_n, 1, MPI_INT, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( MPI_IN_PLACE, m_data_t, m_n0 * m_data_size, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( MPI_IN_PLACE, m_data_z, m_n0 * m_data_size, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( MPI_IN_PLACE, m_tip_z, m_nn * m_nt, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() );
        else
            error = Uintah::MPI::Reduce ( &m_tip_n, nullptr, 1, MPI_INT, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( m_data_t, nullptr, m_n0 * m_data_size, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( m_data_z, nullptr, m_n0 * m_data_size, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() ) ||
                    Uintah::MPI::Reduce ( m_tip_z, nullptr, m_nn * m_nt, MPI_DOUBLE, MPI_MAX, 0, myworld->getComm() );

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessor::reduceMPI, done " );

        if ( error )
        {
            DOUT ( true, "ArmPostProcessor::reduceMPI: Uintah::MPI::Allreduce error: " << error );
            SCI_THROW ( InternalError ( "ArmPostProcessor::reduceMPI: MPI error", __FILE__, __LINE__ ) );
        }
    }

    virtual void
    printData (
        std::ostream & out
    ) const override
    {
        out << "n_tip: " << m_tip_n << "\n";
        for ( int j = m_nt - 1; j >= 0; --j )
        {
            out << "tip_data: ";
            for ( int i = 0; i < m_nn; ++i )
                if ( tip_z ( j, i ) == -DBL_MAX )
                    out << "_________ ";
                else
                    out << std::setw ( 9 ) << tip_z ( j, i ) << " ";
            out << "\n";
        }

        for ( int i = m_n0 - 1; i >= 0; --i )
        {
            out << "t: ";
            for ( int j = 0; j < m_data_size; ++j )
                if ( data_t ( j, i ) == -DBL_MAX )
                    out << "___ ";
                else
                    out << std::setw ( 3 ) << data_t ( j, i ) << " ";
            out << "\n";
        }
        for ( int i = m_n0 - 1; i >= 0; --i )
        {
            out << "z: ";
            for ( int j = 0; j < m_data_size; ++j )
                if ( data_z ( j, i ) == -DBL_MAX )
                    out << "_________ ";
                else
                    out << std::setw ( 9 ) << data_z ( j, i ) << " ";
            out << "\n";
        }
    }

    virtual void
    computeTipInfo (
        double & tip_position,
        double tip_curvatures[3],
        const double & gamma_psi = 0
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorParallelTanh::computeTipInfo " );

        // Containers for 0-level
        double * arm_n = scinew double[m_data_size + m_nt];
        double * arm_t2 = scinew double[m_data_size + m_nt];

        // 1. Compute 0-level using data_t/data_z (first m_data_size points)

        Lapack::Tanh1 tanh1 ( gamma_psi, tip_position, m_psetup );
        double n, t;
        for ( int i = 0; i < m_data_size; ++i )
        {
            n = n_coord ( m_data_n0 + i );

            tanh1.fit ( m_n0, data_t ( i ), data_z ( i ), gamma_psi * n * n );
            t = tanh1.zero();

            arm_n[i] = n;
            arm_t2[i] = t * t;

            DOUT ( DBG_PRINT, "plot3(" << n << "+0*x,x,y,'-ok')" );
            DOUT ( DBG_PRINT && i == 0, "hold on" );
            DOUT ( DBG_PRINT, "plot3(" << n << "," << t << ",0,'*k')\n" );
        }

        // 2. 2D Interpolation at tip for position and curvature

        double * tip_n = scinew double[m_nn * m_nt];
        double * tip_t = scinew double[m_nn * m_nt];

        int n0 =  m_tip_n - m_nnl;

        for ( int i = 0; i < m_nt; ++i )
        {
            for ( int j = 0; j < m_nn; ++j )
                tip_n[m_nn * i + j] = n_coord ( n0 + j );
            std::fill ( tip_t + m_nn * i, tip_t + m_nn * ( i + 1 ), t_coord ( i ) );
        }

        Lapack::Tanh2 tanh2 ( gamma_psi, tip_position, m_psetup );
        tanh2.fit ( m_nn * m_nt, tip_n, tip_t, tip_z ( 0 ) );

        DOUT ( DBG_PRINT, "nn=" << m_nn << "; nt=" <<  m_nt << ";" );
        DOUT ( DBG_PRINT, "mesh(reshape(x,nn,nt),reshape(y,nn,nt),reshape(z,nn,nt),'FaceColor','none')" );

        for ( int i = 0; i < m_nt; ++i )
        {
            t = t_coord ( i );
            n = tanh2.zero ( t );
            arm_n[m_data_size + m_nt - i - 1] = n;
            arm_t2[m_data_size + m_nt - i - 1] = t * t;

            DOUT ( DBG_PRINT, "plot3(" << n << "," << t << ",0,'ro')" );
        }

        tip_position = tanh2.zero ( 0. );
        tip_curvatures[0] = 1. / tanh2.r0();

        DOUT ( DBG_PRINT, "plot3(" << tip_position << ",0,0,'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)" );

        delete[] tip_n;
        delete[] tip_t;

        // 3. Evaluate parabolic curvatue using full 0-level

        int skip = 0;
        while ( arm_t2[skip] <= arm_t2[skip + 1] && skip < m_data_size - 1 ) ++skip;

        Lapack::Poly parabola ( 1 );
        double kn, kd;
        parabola.fit ( m_data_size + m_nt - skip, arm_t2 + skip, arm_n + skip );
        kn = 2.* parabola.cfx ( 1 );
        kd = 1 - 4.* parabola.cfx ( 1 ) * parabola.cfx ( 2 );
        tip_curvatures[1] = std::sqrt ( ( kn * kn ) / ( kd * kd * kd ) );

        DOUT ( DBG_PRINT, "n=y;" );
        DOUT ( DBG_PRINT, "t1=sqrt(x);" );
        DOUT ( DBG_PRINT, "t2=-t1;" );
        DOUT ( DBG_PRINT, "tt=linspace(t2(1),t1(1));" );
        DOUT ( DBG_PRINT, "plot3(n,t1,0*n,'k.');" );
        DOUT ( DBG_PRINT, "plot3(n,t2,0*n,'k.');" )
        DOUT ( DBG_PRINT, "plot3(polyval(p,tt.^2),tt,0*tt,'r-');\n" );

        // 4. Evaluate parabolic curvatue excluding tip neighbour

        int arm_size = m_data_size;
        for ( ; arm_size > 0; --arm_size )
        {
            const double & t2 = arm_t2[arm_size - 1];
            const double & tn = tip_position - arm_n[arm_size - 1];
            if ( t2 < tn * tn ) break;
        }

        if ( arm_size - skip < 2 )
            tip_curvatures[2] = NAN;
        else
        {
            parabola.fit ( arm_size - skip, arm_t2 + skip, arm_n + skip );
            kn = 2.* parabola.cfx ( 1 );
            kd = 1 - 4.* parabola.cfx ( 1 ) * parabola.cfx ( 2 );
            tip_curvatures[2] = std::sqrt ( ( kn * kn ) / ( kd * kd * kd ) );

            DOUT ( DBG_PRINT, "n=y;" );
            DOUT ( DBG_PRINT, "t1=sqrt(x);" );
            DOUT ( DBG_PRINT, "t2=-t1;" );
            DOUT ( DBG_PRINT, "tt=linspace(t2(1),t1(1));" );
            DOUT ( DBG_PRINT, "plot3(n,t1,0*n,'b.');" );
            DOUT ( DBG_PRINT, "plot3(n,t2,0*n,'b.');" )
            DOUT ( DBG_PRINT, "plot3(polyval(p,tt.^2),tt,0*tt,'m-');\n" );
        }

        delete[] arm_n;
        delete[] arm_t2;
    }
}; // class ArmPostProcessorParallelTanhD2

template<VarType VAR> const IntVector ArmPostProcessorParallelTanh < VAR, D2 >::en { 1, 0, 0 };
template<VarType VAR> const IntVector ArmPostProcessorParallelTanh < VAR, D2 >::et { 0, 1, 0 };

template<VarType VAR>
class ArmPostProcessorParallelTanh < VAR, D3 >
    : public ArmPostProcessor < VAR, D3 >
{
public:
    ArmPostProcessorParallelTanh ( const Lapack::TrustRegionSetup & psetup, int npts, int ntip, int, bool dbg ) {}
    virtual ~ArmPostProcessorParallelTanh() {};
    virtual void setLevel ( const Level * level ) {}
    virtual void initializeLocations () {};
    virtual void setLocations ( IntVector const & low, IntVector const & high, View < ScalarField<const double> > const & psi ) {};
    virtual void reduceLocations ( const ProcessorGroup * myworld ) {};
    virtual void printLocations ( std::ostream & out ) const {};
    virtual void initializeData () {};
    virtual void setData ( IntVector const & low, IntVector const & high, View < ScalarField<const double> > const & psi ) {};
    virtual void reduceData ( const ProcessorGroup * myworld ) {};
    virtual void printData ( std::ostream & out ) const {};
    virtual void computeTipInfo ( double & tip_position, double tip_curvatures[3], const double & gamma_psi ) {};
}; // class ArmPostProcessorParallelTanhD3

} // namespace PhaseField
} // namespace Uintah

#undef DBG_PRINT

#endif // Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorParallelTanh_h
