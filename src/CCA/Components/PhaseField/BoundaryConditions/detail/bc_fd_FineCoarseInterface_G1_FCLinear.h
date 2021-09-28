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
 * @file CCA/Components/PhaseField/BoundaryConditions/detail/bc_fd_FineCoarseInterface_G1_FCLinear.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bc_fd_FineCoarseInterface_G1_FCLinear_h
#define Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bc_fd_FineCoarseInterface_G1_FCLinear_h

#include <CCA/Components/PhaseField/Views/detail/view.h>
#include <CCA/Components/PhaseField/Views/detail/basic_fd_view.h>
#include <CCA/Components/PhaseField/DataWarehouse/detail/dw_fd.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWInterface.h>
#include <CCA/Components/PhaseField/AMR/AMRInterface.h>

namespace Uintah
{
namespace PhaseField
{
namespace detail
{

/**
 * @brief Finite-differences scheme at fine/coarse interfaces (1-ghost stencils)
 * (ScalarField implementation)
 *
 * Implement FCLinear approximation of ghost value across the amr interface
 *
 * @image latex fclinear.eps "FCLinear"
 * @image html  fclinear.png "FCLinear"
 *
 * @todo generilize to DIM != D2
 *
 * @tparam T type of the field value at each point
 * @tparam STN finite-difference stencil
 * @tparam VAR type of variable representation
 * @tparam F patch face on which bc is to be applied
 */
template<typename T, StnType STN, VarType VAR, Patch::FaceType F>
class bc_fd < ScalarField<T>, STN, VAR, F, BC::FineCoarseInterface, 1, FC::FCLinear >
    : virtual public basic_fd_view < ScalarField<T>, STN >
{
private: // STATIC MEMBERS

    /// Number of ghosts required
    static constexpr int GN = 1;

    /// Problem dimension
    static constexpr DimType DIM = get_stn<STN>::dim;

    /// Boundary face normal vector direction
    static constexpr DirType D = get_face<F>::dir;

    /// Boundary face normal vector sign (int)
    static constexpr int SGN = get_face<F>::sgn;

    /// Boundary face normal vector sign (double)
    static constexpr double DSGN = get_face<F>::dsgn;

private: // TYPES

    /// Type of field
    using Field = ScalarField<T>;

    /// Non const type of the field value
    using V = typename std::remove_const<T>::type;

#ifdef HAVE_HYPRE
    /// Stencil entries type
    using S = typename get_stn<STN>::template type<T>;

    using A = HypreSStruct::AdditionalEntries;
#endif

private: // STATIC ASSERTIONS

    static_assert ( DIM == D2, "FCLinear must be used only for 2D problems" );

private: // MEMBERS

    /// View of the fine level for grid values
    const view<Field> * m_fine_view;

    /// View of the coarse level for computing fine ghost values
    view<Field> * m_coarse_view;

    /// Fine grid level
    const Level * m_level_fine;

    /// Coarse grid level
    const Level * m_level_coarse;

    /// Grid spacing of fine level
    Vector m_h;

    /// Region where the view is defined
    Support m_support;

private: // METHODS

    /**
     * @brief Get fine value at position
     *
     * @param id fine index
     * @return fine value at id
     */
    inline V
    fine_value (
        const IntVector & id
    ) const
    {
        return ( *m_fine_view ) [ id ];
    };

#ifdef HAVE_HYPRE
    inline Entries<V>
    fine_entries (
        const IntVector & id
    ) const
    {
        return m_fine_view->entries ( id );
    };
#endif

    /**
     * @brief Get coarse value at position
     *
     * @param id coarse index
     * @return fine value at id
     */
    inline V
    coarse_value (
        const IntVector & id
    ) const
    {
        const auto & coarse_view = *m_coarse_view;
        return coarse_view[id];
    };

#ifdef HAVE_HYPRE
    inline Entries<V>
    coarse_entries (
        const IntVector & id
    ) const
    {
        return m_coarse_view->entries ( id );
    };
#endif

    /**
     * @brief Get coarse interpolation at position
     *
     * @param id position index
     * @return coarse field interpolated value at id
     */
    inline V
    coarse_interp (
        const IntVector & id
    ) const
    {
        return operator[] ( id );
    };

#ifdef HAVE_HYPRE
    inline Entries<V>
    coarse_interp_entries (
        const IntVector & id
    ) const
    {
        return entries ( id );
    };
#endif

protected:

    /**
     * @brief Constructor
     *
     * Instantiate a copy of a given view
     *
     * @param copy source view for copying
     * @param deep if true inner grid variable is copied as well otherwise the
     * same grid variable is referenced
     */
    bc_fd (
        const bc_fd * copy,
        bool deep
    ) : m_fine_view ( copy->m_fine_view ),
        m_coarse_view ( copy->m_coarse_view->clone ( deep ) ),
        m_level_fine ( copy->m_level_fine ),
        m_level_coarse ( copy->m_level_coarse ),
        m_h ( copy->m_h ),
        m_support ( copy->m_support )
    {}


    /**
     * @brief Constructor
     *
     * Instantiate a copy of a given view using the specified fine view to
     * access the DataWarehouse.
     *
     * @remark When a bc_fd is instantiated from a finer level to enforce fine/coarse
     * interface conditions the copy fine view will not have access to the
     * DataWarehouse after set (bc_fd::set method is not retrieving data for
     * fine level)
     * @param fine_view view to be used for accessing the DataWarehouse
     * @param copy source view for copying
     * @param deep if true inner grid variable is copied as well otherwise the
     * same grid variable is referenced
     */
    bc_fd (
        const view<Field> * fine_view,
        const bc_fd * copy,
        bool deep
    ) : m_fine_view ( fine_view ),
        m_coarse_view ( copy->m_coarse_view->clone ( deep ) ),
        m_level_fine ( copy->m_level_fine ),
        m_level_coarse ( copy->m_level_coarse ),
        m_h ( copy->m_h ),
        m_support ( copy->m_support )
    {}

protected: // CONSTRUCTORS

    /**
     * @brief Constructor
     *
     * Instantiate a view without gathering info from the DataWarehouse
     *
     * @param fine_view view to access fine grid values
     * @param label unused since value is already in fine/coarse_view
     * @param material unused since value is already in fine/coarse_view
     * @param level fine grid level
     * @param coarse_view view to access coarse grid values
     */
    bc_fd (
        const view<Field> * fine_view,
        const VarLabel * _DOXYARG ( label ),
        int _DOXYARG ( material ),
        const Level * level,
        view<Field> * coarse_view
    ) : m_fine_view ( fine_view ),
        m_coarse_view ( coarse_view ),
        m_level_fine ( level ),
        m_level_coarse ( level->getCoarserLevel().get_rep() ),
        m_h ( level->dCell() )
    {}

public: // CONSTRUCTORS/DESTRUCTOR

    /// Destructor
    virtual ~bc_fd()
    {
        delete m_coarse_view; // created by bcs_basic_fd_view::get_value
    };

    /// Prevent copy (and move) constructor
    bc_fd ( const bc_fd & ) = delete;

    /// Prevent copy (and move) assignment
    /// @return deleted
    bc_fd & operator= ( const bc_fd & ) = delete;

public: // VIEW METHODS

    /**
     * @brief Retrieve value from the DataWarehouse for a given patch
     *
     * @remark bc_fd should not be set to range over a patch
     *
     * @param dw unused
     * @param patch unused
     * @param use_ghosts unused
     */
    virtual void
    set (
        DataWarehouse * _DOXYARG ( dw ),
        const Patch * _DOXYARG ( patch ),
        bool _DOXYARG ( use_ghosts )
    ) override
    {
        ASSERTFAIL ( "cannot set bc_fd over a patch" );
    };

    /**
     * @brief Retrieve value from the DataWarehouse for a given region
     *
     * retrieve also info from coarser level region by setting the coarser
     * view
     *
     * @param dw DataWarehouse from which data is retrieved
     * @param level grid level from which retrieve data
     * @param low lower bound of the region to retrieve
     * @param high higher bound of the region to retrieve
     * @param use_ghosts if ghosts value are to be retrieved
     */
    virtual void
    set (
        DataWarehouse * dw,
        const Level * level,
        const IntVector & low,
        const IntVector & high,
        bool use_ghosts
    ) override
    {
        m_support.clear();

        IntVector l ( low ), h ( high );

        if ( SGN > 0 ) // F=D+
        {
            l[D] = high[D];
            if ( use_ghosts ) h[D] += GN;
        }
        else // F=D-
        {
            h[D] = low[D];
            if ( use_ghosts ) l[D] -= GN;
        }

        m_support.emplace_back ( l, h );

        m_level_fine = level;
        m_level_coarse = level->getCoarserLevel().get_rep();
        m_h = level->dCell();

        if ( use_ghosts )
        {
            m_coarse_view->set ( dw, level, l, h );
        }
    };

    /**
     * @brief Get/Modify value at position with index id (virtual implementation)
     *
     * @remark value at boundary is computed at runtime thus doesn't exist in the
     * DataWarehouse
     *
     * @param id unused
     * @return nothing
     */
    virtual T &
    operator[] (
        const IntVector & _DOXYARG ( id )
    ) override VIRT;

    /**
     * @brief Get value at position with index id
     *
     * @param id position index (on fine level)
     * @return interpolated value from coarser level
     */
    virtual V
    operator[] (
        const IntVector & id
    ) const override
    {
        IntVector id0 = AMRInterface<VAR, DIM>::get_coarser ( m_level_fine, id );
        Point p_fine = DWInterface<VAR, DIM>::get_position ( m_level_fine, id );
        Point p0 = DWInterface<VAR, DIM>::get_position ( m_level_coarse, id0 );
        double w0, w1, w2;
        IntVector id1 { id }, id2 { id };
        id1[D] = id2[D] -= SGN;

        for ( size_t d = 0; d < DIM; ++d )
            if ( d != D )
            {
                if ( p0.asVector() [d] < p_fine.asVector() [d] )
                    id2[d] -= 1;
                else
                    id2[d] += 1;
            }
        w0 = 2. / 3.;
        w1 = 2. / 3.;
        w2 = -1. / 3.;

        V u0 = coarse_value ( id );
        V u1 = fine_value ( id1 );
        V u2 = fine_value ( id2 );

        return w0 * u0 + w1 * u1 + w2 * u2;
    };

#ifdef HAVE_HYPRE
    virtual Entries<V>
    entries (
        const IntVector & id_fine
    ) const override
    {
        IntVector id0 = AMRInterface<VAR, DIM>::get_coarser ( m_level_fine, id_fine );
        Point p_fine = DWInterface<VAR, DIM>::get_position ( m_level_fine, id_fine );
        Point p0 = DWInterface<VAR, DIM>::get_position ( m_level_coarse, id0 );
        double w0, w1, w2;
        IntVector id1 { id_fine }, id2 { id_fine };
        id1[D] = id2[D] -= SGN;

        for ( size_t d = 0; d < DIM; ++d )
            if ( d != D )
            {
                if ( p0.asVector() [d] < p_fine.asVector() [d] )
                    id2[d] -= 1;
                else
                    id2[d] += 1;
            }
        w0 = 2. / 3.;
        w1 = 2. / 3.;
        w2 = -1. / 3.;

        Entries<V> res;
        res.add ( coarse_entries ( id_fine ), w0 );
        res.add ( fine_entries ( id1 ), w1 );
        res.add ( fine_entries ( id2 ), w2 );

        res.simplify();
        return res;
    };
#endif

    /**
     * @brief Get the region on which the view is defined
     *
     * It is the ghost fine region across an amr interface
     *
     * @return support of the view
     */
    virtual
    Support
    get_support()
    const override
    {
        return m_support;
    };

    /**
     * @brief Check if the view has access to the fine position with index id
     *
     * @param id fine position index
     * @return check result
     */
    virtual
    bool
    is_defined_at (
        const IntVector & id
    ) const override
    {
        IntVector low { m_support.front().getLow() }, high { m_support.front().getHigh() };
        return ( low[X] <= id[X] && id[X] < high[X] ) &&
               ( low[Y] <= id[Y] && id[Y] < high[Y] ) &&
               ( low[Z] <= id[Z] && id[Z] < high[Z] );
    };

public: // BASIC FD VIEW METHODS

    /**
     * @brief Get base view (virtual implementation)
     *
     * @remark bc values are computed at runtime thus cannot create non const
     * view, there's nothing to modify in the DataWarehouse
     *
     * @return nothing
     */
    virtual inline view<Field> *
    get_view()
    override VIRT;

    /**
     * @brief Get base view
     *
     * @return const pointer to base view implementation
     */
    virtual inline const view<Field> *
    get_view()
    const override
    {
        return this;
    };

public: // BC FD MEMBERS

    /**
     * @brief First order derivative
     * (parallel direction implementation)
     *
     * First order derivative along DIR at index id at amr interface
     *
     * @remark derivative parallel to the face should not be computed by this bc_fd
     *
     * @tparam DIR Direction along with derivative is approximated
     * @param id unused
     * @return nothing
     */
    template < DirType DIR >
    inline typename std::enable_if < D != DIR, T >::type
    d (
        const IntVector & _DOXYARG ( id )
    ) const VIRT;

    /**
     * @brief First order derivative
     *
     * First order derivative along DIR at index id at amr interface
     *
     * @tparam DIR Direction along with derivative is approximated
     * @param id index where to evaluate the finite-difference
     * @return approximated value at id
     */
    template < DirType DIR >
    inline typename std::enable_if < D == DIR, T >::type
    d (
        const IntVector & id
    ) const
    {
        IntVector ip ( id ), im ( id );
        ip[D] += SGN;
        im[D] -= SGN;
        return 0.5 * DSGN * ( coarse_interp ( ip ) - fine_value ( im ) ) / m_h[D];
    }

    /**
     * @brief Second order derivative
     * (parallel direction implementation)
     *
     * Second order derivative along DIR at index id at amr interface
     *
     * @remark derivative parallel to the face should not be computed by this bc_fd
     *
     * @tparam DIR Direction along with derivative is approximated
     * @param id unused
     * @return nothing
     */
    template < DirType DIR >
    inline typename std::enable_if < D != DIR, T >::type
    d2 (
        const IntVector & _DOXYARG ( id )
    ) const VIRT;

    /**
     * @brief Second order derivative
     * (normal direction implementation)
     *
     * Second order derivative along DIR at index id at amr interface
     *
     * @tparam DIR Direction along with derivative is approximated
     * @param id index where to evaluate the finite-difference
     * @return approximated value at id
     */
    template < DirType DIR >
    inline typename std::enable_if < D == DIR, T >::type
    d2 (
        const IntVector & id
    ) const
    {
        IntVector ip ( id ), im ( id );
        ip[D] += SGN;
        im[D] -= SGN;
        return ( fine_value ( im ) + coarse_interp ( ip ) - 2. * fine_value ( id ) )
               / ( m_h[D] * m_h[D] );
    }

#ifdef HAVE_HYPRE
    template < DirType DIR >
    inline typename std::enable_if < D != DIR, void >::type
    add_d2_sys_hypre (
        const IntVector & _DOXYARG ( id ),
        S & _DOXYARG ( stencil_entries ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D == DIR, void >::type
    add_d2_sys_hypre (
        const IntVector &,
        S & _DOXYARG ( stencil_entries ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D != DIR, void >::type
    add_d2_rhs_hypre (
        const IntVector & _DOXYARG ( id ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D == DIR, void >::type
    add_d2_rhs_hypre (
        const IntVector & _DOXYARG ( id ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D != DIR, void >::type
    add_d2_sys_hypresstruct (
        const IntVector & _DOXYARG ( id ),
        S & _DOXYARG ( stencil_entries ),
        A & _DOXYARG ( additional_entries ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D == DIR, void >::type
    add_d2_sys_hypresstruct (
        const IntVector & id,
        S & stencil_entries,
        A & additional_entries,
        V & rhs
    ) const
    {
        if ( VAR == CC )
        {
            IntVector ip ( id );
            ip[D] += SGN;
            double h2 = m_h[D] * m_h[D];
            stencil_entries[F - SGN] += 1. / h2;
            stencil_entries.p += -2. / h2;

            Entries<V> extra = coarse_interp_entries ( ip );
            for ( auto & entry : extra.values )
            {
                if ( entry.level == m_level_fine->getIndex() ) // check if entry belongs to stencil
                {
                    IntVector tmp = entry.index - id;
                    if ( entry.index == id )
                    {
                        stencil_entries.p += entry.weight / h2;
                        goto next_entry;
                    }
                    for ( Patch::FaceType f = get_dim<DIM>::face_start; f <= get_dim<DIM>::face_end; f = Patch::nextFace ( f ) )
                        if ( tmp == Patch::getFaceDirection ( f ) )
                        {
                            stencil_entries[f] += entry.weight / h2;
                            goto next_entry;
                        }
                next_entry:
                    continue;
                }

                HypreSStruct::MatrixIndex index { id, entry.level, entry.index };
                auto it = additional_entries.find ( index );
                if ( it != additional_entries.end() ) it->second += entry.weight / h2;
                else additional_entries.emplace ( index, entry.weight / h2 );
            }
            rhs += extra.rhs;
        }
        else ASSERTFAIL ( "TODO" );
    }

    template < DirType DIR >
    inline typename std::enable_if < D != DIR, void >::type
    add_d2_rhs_hypresstruct (
        const IntVector & _DOXYARG ( id ),
        V & _DOXYARG ( rhs )
    ) const VIRT;

    template < DirType DIR >
    inline typename std::enable_if < D == DIR, void >::type
    add_d2_rhs_hypresstruct (
        const IntVector & id,
        V & rhs
    ) const
    {
        if ( VAR == CC )
        {
            IntVector ip ( id );
            ip[D] += SGN;
            Entries<V> extra = coarse_interp_entries ( ip );
            rhs += extra.rhs;
        }
        else ASSERTFAIL ( "TODO" );
    };
#endif

}; // class bc_fd

} // namespace detail
} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bc_fd_FineCoarseInterface_G1_FCLinear_h
