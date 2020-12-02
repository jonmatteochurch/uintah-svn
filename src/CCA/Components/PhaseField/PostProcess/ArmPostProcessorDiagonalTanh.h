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
 * @file CCA/Components/PhaseField/PostProcess/ArmPostProcessorDiagonalTanh.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/06
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorDiagonalTanh_h
#define Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorDiagonalTanh_h

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

template<VarType VAR>
class ArmPostProcessorDiagonalTanh
    : public ArmPostProcessor < VAR, D2 >
    , public ArmPostProcessor < VAR, D3 >
{
private: // STATIC MEMBERS

    static constexpr TypeDescription::Type var_td = VAR == CC ? TypeDescription::CCVariable : TypeDescription::NCVariable;
    static const IntVector en;
    static const IntVector et;

protected: // MEMBERS

    const bool m_dbg;

    const double m_alpha;

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
        return ( id[0] - m_origin[0] ) + ( id[1] - m_origin[1] );
    }

    int n_ind ( const int & x, const int & y ) const
    {
        return ( x - m_origin[0] ) + ( y - m_origin[1] );
    }

    int t_ind ( const IntVector & id ) const
    {
        return ( id[0] - m_origin[0] ) - ( id[1] - m_origin[1] );
    }

    int t_ind ( const int & x, const int & y ) const
    {
        return ( x - m_origin[0] ) - ( y - m_origin[1] );
    }

    int x_ind ( const IntVector & id ) const
    {
        return ( id[0] + id[1] ) / 2 + m_origin[0];
    }

    int x_ind ( const int & n, const int & t ) const
    {
        return ( n + t ) / 2 + m_origin[0];
    }

    int y_ind ( const IntVector & id ) const
    {
        return ( id[0] - id[1] ) / 2 + m_origin[1];
    }

    int y_ind ( const int & n, const int & t ) const
    {
        return ( n - t ) / 2 + m_origin[1];
    }

    template < VarType V = VAR, typename std::enable_if<V == CC, int>::type = 0 >
    double n_coord ( const int & n )
    {
        return ( n + 1. ) * m_dn;
    }

    template < VarType V = VAR, typename std::enable_if<V == NC, int>::type = 0 >
    double n_coord ( const int & n )
    {
        return n * m_dn;
    }

    template < VarType V = VAR, typename std::enable_if<V == CC, int>::type = 0 >
    double t_coord ( const int & t )
    {
        return t * m_dt;
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

    ArmPostProcessorDiagonalTanh (
        const Lapack::TrustRegionSetup psetup,
        int n0,
        int nn,
        int nt,
        double alpha,
        bool dbg
    ) : m_dbg ( dbg ),
        m_alpha ( alpha ),
        m_n0 ( n0 ),
        m_n0l ( ( m_n0 - 1 ) / 2 ),
        m_n0h ( m_n0 - m_n0l ),
        m_nn ( nn ),
        m_nt ( nt ),
        m_nnl ( nn ),
        m_nnh ( 2 * m_nn - m_nnl ),
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
    virtual ~ArmPostProcessorDiagonalTanh ()
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
        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::setLevel: " << level->getIndex() << " " );

        m_origin = level->getCellIndex ( {0., 0., 0.} );
        level->computeVariableExtents ( var_td, m_low, m_high );
        m_dn = m_dt = 0.5 * std::sqrt ( level->dCell() [0] * level->dCell() [0] + level->dCell() [1] * level->dCell() [1] );
    }

    virtual void
    initializeLocations (
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::initializeLocations " );

        m_location_n0 = std::max ( 0, n_ind ( m_low ) );
        m_locations_size = n_ind ( m_high ) - m_location_n0;
        m_locations = scinew int[m_locations_size];
        std::fill ( m_locations, m_locations + m_locations_size, INT_MAX );
    };

    virtual void
    setLocations (
        const Patch * patch,
        IntVector const & low,
        IntVector const & high,
        const std::list<Patch::FaceType> & faces,
        View < ScalarField<const double> > const & psi
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::setLocations: " << low << high << " " );

        bool fc_yminus = patch->getBCType ( Patch::yminus ) == Patch::Coarse &&
                         std::find ( faces.begin(), faces.end(), Patch::yminus ) != faces.end();
        bool fc_xplus = patch->getBCType ( Patch::xplus ) == Patch::Coarse &&
                        std::find ( faces.begin(), faces.end(), Patch::xplus ) != faces.end();

        IntVector dl {0, fc_yminus ? 1 : 0, 0};
        IntVector dh {fc_xplus ? 1 : 0, 0, 0};
        IntVector inf = Max ( low + dl, m_origin );
        IntVector sup = Min ( high - dh, m_high - 1 );

        if ( inf[0] < inf[1] )
        {
            inf[0] = inf[1];
        }

        int iy0 = m_origin[1]; // symmetric
        int ix1 = m_high[0] - 1; // extend

        // move along increasing n
        bool found;
        int it1;
        IntVector id;
        int in, it;
        for ( in = n_ind ( inf ); in < n_ind ( sup ); ++in )
        {
            found = false;
            it1 = std::min (
            {
                in - 2 * low[1], 2 * ( sup[0] - 1 ) - in,
                in - 2 * ( iy0 + 1 ), 2 * ( ix1 - 1 ) - in,
            } );

            id = inf;
            // move along increasing t
            for ( it = t_ind ( inf ); it <= it1; it += 2, id += et )
                if ( ( found = psi[id] * psi[id + et] <= 0. ) )
                    if ( it < m_locations[in - m_location_n0] )
                    {
                        // if convave edge id+et may not be refined
                        if ( patch->getLevel()->containsCell ( id + et ) )
                            m_locations[in - m_location_n0] = it;
                        break;
                    }

            if ( !found )
            {
                if ( id[0] == ix1 )
                {
                    // extends to -1
                    if ( psi[id] >= 0. )
                        if ( it < m_locations[in - m_location_n0] )
                        {
                            m_locations[in - m_location_n0] = it;
                        }
                }
                else if ( id[1] == iy0 )
                {
                    // symmetry on x
                    IntVector sym { id[0] + 1, id[1] + ( VAR == NC ? 1 : 0 ), 0 };
                    if ( sym[0] < sup[0] && psi[id] * psi[sym] <= 0. )
                        if ( it < m_locations[in - m_location_n0] )
                        {
                            m_locations[in - m_location_n0] = it;
                        }
                }
            }

            // increment n
            // - first along left edge till I bisector or top edge
            // - then if I bisector intersect patch along it
            // - then along top edge
            if ( inf[1] < inf[0] && inf[1] < sup[1] - 1 )
            {
                ++inf[1];
            }
            else
            {
                ++inf[0];
            }
        }
        while ( inf[0] < sup[0] );
    }

    virtual void
    reduceLocations (
        const ProcessorGroup * myworld
    ) override
    {
        if ( myworld->nRanks() <= 1 )
        {
            return;
        }

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessorDiagonalTanh::reduceMPI " );

        int error = Uintah::MPI::Allreduce ( MPI_IN_PLACE, m_locations, m_locations_size, MPI_INT, MPI_MIN, myworld->getComm() );

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessorDiagonalTanh::reduceMPI, done " );

        if ( error )
        {
            DOUT ( true, "ArmPostProcessorDiagonalTanh::reduceMPI: Uintah::MPI::Allreduce error: " << error );
            SCI_THROW ( InternalError ( "ArmPostProcessorDiagonalTanh::reduceMPI: MPI error", __FILE__, __LINE__ ) );
        }
    }

    virtual void
    printLocations (
        std::ostream & out
    ) const override
    {
        for ( int i = 0; i < m_locations_size; ++i )
            if ( m_locations[i] == INT_MAX )
            {
                out << "___ ";
            }
            else
            {
                out << std::setw ( 3 ) << m_locations[i] << " ";
            }
    }

    virtual void
    initializeData ()
    override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::initializeData " );

        int i = m_locations_size;
        for ( i = m_locations_size; i > 0 && m_locations[i - 1] == INT_MAX; --i );
        m_data_size = i--;

        if ( m_data_size )
        {
            for ( ; i > 1 && m_locations[i - 2] < INT_MAX && m_locations[i - 2] >= m_locations[i]; --i );

            m_data_n0 = i - 1;
            m_data_size -= m_data_n0;
        }
        else
        {
            m_data_n0 = -1;
        }

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
        View < ScalarField<const double> > const & psi,
        View< ScalarField<int> > * refine_flag
    ) override
    {
        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::setData: " << low << high << " " );

        // arm contour not here
        if ( !m_data_size || n_ind ( high ) < m_data_n0 - m_nnl )
        {
            return;
        }

        int ix0 = m_origin[0]; // symmetric
        int iy0 = m_origin[1]; // symmetric
        int ix1 = m_high[0]; // extend
        int iy1 = m_high[1]; // extend

        IntVector id, sym;
        int in0, in1;
        int it0, it1;
        int in, it, tr;

        int n_ = -1, t_; // tip candidates

        for ( in = std::max ( m_data_n0, n_ind ( low ) - m_n0h ); in < n_ind ( high ); ++in )
        {
            t_ = m_locations[ in - m_location_n0 ];
            if ( t_ < INT_MAX )
            {
                it = t_ - 2 * m_n0l;

                it0 = std::max (
                {
                    it,
                    in - 2 * ( high[1] - 1 ), 2 * low[0] - in,
                    in - 2 * ( iy1 - 1 ), 2 * ix0 - in,
                } );

                it1 = std::min (
                {
                    it + 2 * m_n0 - 1,
                    in - 2 * low[1], 2 * ( high[0] - 1 ) - in,
                    in - 2 * iy0, 2 * ( ix1 - 1 ) - in,
                } );

                id[0] = x_ind ( in, it );
                id[1] = y_ind ( in, it );
                id[2] = 0;

                n_ = in - m_data_n0;
                t_ = 0;

                // extend psi to -1 out computational boundary
                for ( ; id[1] >= iy1; it += 2, id += et, ++t_ )
                {
                    data_t ( n_, t_ ) = t_coord ( it );
                    data_z ( n_, t_ ) = -1.;
                }

                // symmetry on x axis
                sym = { m_origin[0] - id[0], id[1], 0 };
                if ( VAR == CC )
                {
                    --sym[0];
                }
                for ( ; id[0] < ix0; it += 2, id += et, --sym[0], --sym[1], ++t_ )
                    if ( sym[0] >= low[0] && sym[1] >= low[1] && sym[1] < high[1] )
                    {
                        data_t ( n_, t_ ) = t_coord ( it );
                        data_z ( n_, t_ ) = psi[sym];
                        if ( refine_flag ) ( *refine_flag ) [sym] = -2;
                    }

                // skip non patch region
                for ( ; it < it0; it += 2, id += et, ++t_ );

                // set patch data
                for ( ; it <= it1; it += 2, id += et, ++t_ )
                {
                    data_t ( n_, t_ ) = t_coord ( it );
                    data_z ( n_, t_ ) = psi[id];
                    if ( refine_flag ) ( *refine_flag ) [id] = 2;
                }

                // symmetry on y axis
                sym = { id[0], m_origin[1] - id[1], 0 };
                if ( VAR == CC )
                {
                    --sym[1];
                }
                for ( ; id[0] < ix1 && t_ < m_n0; it += 2, id += et, ++t_, ++sym[0], ++sym[1] )
                    if ( sym[0] < high[0] && low[1] <= sym[1] && sym[1] < high[1] )
                    {
                        data_t ( n_, t_ ) = t_coord ( it );
                        data_z ( n_, t_ ) = psi[sym];
                        if ( refine_flag ) ( *refine_flag ) [sym] = -2;
                    }

                // extend psi to -1 out computational boundary
                for ( ; t_ < m_n0; it += 2, id += et, ++t_, ++sym[0], ++sym[1] )
                {
                    data_t ( n_, t_ ) = t_coord ( it );
                    data_z ( n_, t_ ) = -1.;
                }
            }
        }

        // check for tip
        int dt = m_nt;
        if ( ( t_ind ( high[0] - 1, low[1] ) + dt ) * ( t_ind ( low[0], high[1] - 1 ) - dt ) <= 0 )
        {
            bool is_tip = false;
            if ( m_tip_n < 0 )
            {
                if ( n_ < 0 ) // move back from low[0] at most m_nnh steps
                {
                    in = n_ind ( low ) - m_location_n0;
                    int dn0 = std::max ( -m_nnh, -in );
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
                        {
                            n_ = i;
                        }
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
                in += tr = in % 2; // if lowest is odd increase n to be even
                it = 0;

                IntVector start { x_ind ( in, it ), y_ind ( in, it ), 0 };

                for ( t_ = 0; t_ < m_nt; ++t_, ++it )
                {
                    id = start;
                    in = n_ind ( start );
                    n_ = 0;

                    in0 = std::max (
                    {
                        in,
                        2 * low[0] - it, 2 * low[1] + it,
                        2 * ( 1 - ix1 ) - it, 2 * ( 1 - iy1 ) + it,
                    } );

                    in1 = std::min (
                    {
                        in + 2 * m_nn,
                        2 * high[0] - it, 2 * high[1] + it,
                        2 * ix1 - it, 2 * iy1 + it,
                    } );

                    // extend psi to -1 out computational boundary
                    for ( ; n_ < m_nn && id[0] < -ix1 + 1; in += 2, id += en, ++n_ )
                    {
                        tip_z ( t_, n_ ) = -1.;
                    }

                    // symmetry on x+y axes
                    sym = { m_origin[0] - id[0], m_origin[1] - id[1], 0 };
                    if ( VAR == CC )
                    {
                        --sym[0];
                        --sym[1];
                    }
                    for ( ; n_ < m_nn && id[0] < ix0; in += 2, id += en, ++n_, sym -= en )
                        if ( sym[0] >= low[0] && sym[1] >= low[1] && sym[1] < high[1] )
                        {
                            tip_z ( t_, n_ ) = psi[sym];
                            if ( refine_flag ) ( *refine_flag ) [sym] = -3;
                        }

                    // symmetry on y axis
                    sym = { id[0], m_origin[1] - id[1], 0 };
                    for ( ; n_ < m_nn && id[1] < iy0; in += 2, id += en, ++n_, sym += et )
                        if ( low[0] <= sym[0] && sym[0] < high[0] && low[1] <= sym[1] && sym[1] < high[1] )
                        {
                            tip_z ( t_, n_ ) = psi[sym];
                            if ( refine_flag ) ( *refine_flag ) [sym] = -3;
                        }

                    // skip non patch region
                    for ( ; n_ < m_nn && in < in0; in += 2, id += en, ++n_ );

                    // set patch data
                    for ( ; n_ < m_nn && in < in1; in += 2, id += en, ++n_ )
                        if ( low[0] <= id[0] && id[0] < high[0] && low[1] <= id[1] && id[1] < high[1] )
                        {
                            tip_z ( t_, n_ ) = psi[id];
                            if ( refine_flag ) ( *refine_flag ) [id] = 3;
                        }

                    // skip non patch region
                    for ( ; n_ < m_nn && id[0] < ix1; in += 2, id += en, ++n_ );

                    // skip non patch region
                    for ( ; n_ < m_nn && id[1] < iy1; in += 2, id += en, ++n_ );

                    // extend psi to -1 out computational boundary
                    for ( ; n_ < m_nn; in += 2, id += en, ++n_ )
                    {
                        tip_z ( t_, n_ ) = -1.;
                    }

                    // increment t
                    if ( ( it + tr ) % 2 )
                    {
                        --start[1];
                    }
                    else
                    {
                        ++start[0];
                    }
                }
            }
        }
    }

    virtual void
    reduceData (
        const ProcessorGroup * myworld
    ) override
    {
        if ( !m_data_size || myworld->nRanks() <= 1 )
            return;

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessorDiagonalTanh::reduceMPI " );

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

        DOUT ( g_mpi_dbg, "Rank-" << myworld->myRank() << " ArmPostProcessorDiagonalTanh::reduceMPI, done " );

        if ( error )
        {
            DOUT ( true, "ArmPostProcessorDiagonalTanh::reduceMPI: Uintah::MPI::Allreduce error: " << error );
            SCI_THROW ( InternalError ( "ArmPostProcessorDiagonalTanh::reduceMPI: MPI error", __FILE__, __LINE__ ) );
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
                {
                    out << "_________ ";
                }
                else
                {
                    out << std::setw ( 9 ) << tip_z ( j, i ) << " ";
                }
            out << "\n";
        }

        for ( int i = m_n0 - 1; i >= 0; --i )
        {
            out << "t: ";
            for ( int j = 0; j < m_data_size; ++j )
                if ( data_t ( j, i ) == -DBL_MAX )
                {
                    out << "___ ";
                }
                else
                {
                    out << std::setw ( 3 ) << data_t ( j, i ) << " ";
                }
            out << "\n";
        }
        for ( int i = m_n0 - 1; i >= 0; --i )
        {
            out << "z: ";
            for ( int j = 0; j < m_data_size; ++j )
                if ( data_z ( j, i ) == -DBL_MAX )
                {
                    out << "_________ ";
                }
                else
                {
                    out << std::setw ( 9 ) << data_z ( j, i ) << " ";
                }
            out << "\n";
        }
    }

    virtual void
    computeTipInfo (
        double & tip_position,
        double tip_curvatures[3]
    ) override
    {
        if ( !m_data_size )
            return;

        DOUTR ( m_dbg,  "ArmPostProcessorDiagonalTanh::computeTipInfo " );

        // Containers for 0-level
        double * arm_n = scinew double[m_data_size + m_nt];
        double * arm_t2 = scinew double[m_data_size + m_nt];

        // 1. Compute 0-level using data_t/data_z (first m_data_size points)

        Lapack::Tanh1 tanh1 ( m_psetup );
        double n, t;
        for ( int i = 0; i < m_data_size; ++i )
        {
            n = n_coord ( m_data_n0 + i );

            tanh1.fit ( m_n0, data_t ( i ), data_z ( i ) );
            t = tanh1.zero();

            arm_n[i] = n;
            arm_t2[i] = t * t;

            DOUT ( DBG_PRINT, "plot3(" << n << "+0*x,x,y,'ok')" );
            DOUT ( DBG_PRINT && i == 0, "hold on" );
            DOUT ( DBG_PRINT, "yy=linspace(x(1),x(end));" );
            DOUT ( DBG_PRINT, "zz=-tanh(beta(1)+beta(2)*yy+beta(3)*yy.^2);" );
            DOUT ( DBG_PRINT, "plot3(" << n << "+0*yy,yy,zz,'-k')\n" );
            DOUT ( DBG_PRINT, "plot3(" << n << "," << t << ",0,'*k')\n" );
        }

        // 2. 2D Interpolation at tip for position and curvature

        double * tip_n = scinew double[m_nn * m_nt];
        double * tip_t = scinew double[m_nn * m_nt];

        int n0 =  m_tip_n - m_nnl;

        for ( int i = 0; i < m_nt; ++i )
        {
            int j0 = ( n0 + i ) % 2;
            for ( int j = j0; j < 2 * m_nn; j += 2 )
            {
                tip_n[m_nn * i + j / 2] = n_coord ( n0 + j );
            }
            std::fill ( tip_t + m_nn * i, tip_t + m_nn * ( i + 1 ), t_coord ( i ) );
        }

        Lapack::Tanh2 tanh2 ( m_psetup );
        tanh2.fit ( m_nn * m_nt, tip_n, tip_t, tip_z ( 0 ) );

        DOUT ( DBG_PRINT, "plot3 (x,y,z,'bo');" );
        DOUT ( DBG_PRINT, "[X,Y]=meshgrid(linspace(x(1),x(end),10),linspace(y(1),y(end),10));" );
        DOUT ( DBG_PRINT, "Z=-tanh(beta(1)+beta(2)*X+beta(3)*X.^2+beta(4)*Y.^2);" );
        DOUT ( DBG_PRINT, "mesh(X,Y,Z);" );

        for ( int i = 0; i < m_nt; ++i )
        {
            t = t_coord ( i );
            n = tanh2.zero ( t );
            arm_n[m_data_size + m_nt - i - 1] = n;
            arm_t2[m_data_size + m_nt - i - 1] = t * t;

            DOUT ( DBG_PRINT, "plot3(" << n << "," << t << ",0,'ro')" );
        }

        tip_position = tanh2.zero ( 0. );
        tip_curvatures[0] = m_dn * std::abs ( ( tanh2.dxx0() + 3. * tanh2.dyy0() ) / ( 2. * M_SQRT2 * tanh2.dx0 ( tip_position ) ) );

        DOUT ( DBG_PRINT, "plot3(" << tip_position << ",0,0,'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10)" );

        delete[] tip_n;
        delete[] tip_t;

        // 3. Evaluate parabolic curvature using full 0-level

        int skip = 0;
        while ( arm_t2[skip] <= arm_t2[skip + 1] && skip < m_data_size - 1 )
        {
            ++skip;
        }

        Lapack::Poly parabola ( 1 );
        if ( m_data_size - skip < 2 )
        {
            tip_curvatures[1] = NAN;
        }
        else
        {
            parabola.fit ( m_data_size - skip, arm_t2 + skip, arm_n + skip );
            tip_curvatures[1] = 2.* std::abs ( parabola.cfx ( 1 ) );

            DOUT ( DBG_PRINT, "n=y;" );
            DOUT ( DBG_PRINT, "t1=sqrt(x);" );
            DOUT ( DBG_PRINT, "t2=-t1;" );
            DOUT ( DBG_PRINT, "tt=linspace(t2(1),t1(1));" );
            DOUT ( DBG_PRINT, "plot3(n,t1,0*n,'k.');" );
            DOUT ( DBG_PRINT, "plot3(n,t2,0*n,'k.');" )
            DOUT ( DBG_PRINT, "plot3(polyval(p,tt.^2),tt,0*tt,'r-');\n" );
        }

        // 4. Evaluate parabolic curvature excluding tip neighbor

        int arm_size = m_data_size - 2;
        skip = std::max ( skip, 1 );
        for ( ; arm_size >= skip; --arm_size )
            if ( arm_t2[arm_size - 1] - 2 * arm_t2[arm_size] + arm_t2[arm_size + 1] < m_alpha * m_dn )
                break;

        if ( arm_size - skip < 2 )
        {
            tip_curvatures[2] = NAN;
        }
        else
        {
            parabola.fit ( arm_size - skip, arm_t2 + skip, arm_n + skip );
            tip_curvatures[2] = 2.* std::abs ( parabola.cfx ( 1 ) );

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
}; // class ArmPostProcessorDiagonalTanh

template<VarType VAR> const IntVector ArmPostProcessorDiagonalTanh < VAR >::en { 1, 1, 0 };
template<VarType VAR> const IntVector ArmPostProcessorDiagonalTanh < VAR >::et { 1, -1, 0 };

} // namespace PhaseField
} // namespace Uintah

#undef DBG_PRINT

#endif // Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorDiagonalTanh_h
