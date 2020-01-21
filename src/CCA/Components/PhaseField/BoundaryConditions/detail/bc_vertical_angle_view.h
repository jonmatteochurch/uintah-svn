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
 * FITNESS FOR stencil_entries PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 * @file CCA/Components/PhaseField/BoundaryConditions/detail/bc_vertical_angle_view.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2019/8
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bc_vertical_angle_view_h
#define Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bc_vertical_angle_view_h

#include <CCA/Components/PhaseField/Util/Expressions.h>
#include <CCA/Components/PhaseField/DataTypes/Support.h>
#include <CCA/Components/PhaseField/Views/detail/view.h>
#include <CCA/Components/PhaseField/Views/detail/basic_fd_view.h>
#include <CCA/Components/PhaseField/AMR/detail/amr_interpolator.h>

#include <memory>

namespace Uintah
{
namespace PhaseField
{
namespace detail
{

/**
 * @brief Detail implementation of variables wrapper for handling ghost regions
 * opposite to boundary edges and vertices
 *
 * Group together multiple views (one for each edge the boundary belongs to,
 * and one for accessing the DataWarehouse on internal indices) to allow
 * extrapolation in regions corresponding to vertical angle
 *
 * @tparam Field type of Field
 * @tparam STN finite-difference stencil
 * @tparam P list of BC, FC, and Patch::Face packs (FineCoarseInterface are not allowed)
 */
template<typename Field, StnType STN, BCF ... P > class bc_vertical_angle_view;

/**
 * @brief Detail implementation of variables wrapper for handling ghost regions
 * opposite to boundary edges and vertices (ScalarField implementation)
 *
 * Group together multiple views (one for each edge the boundary belongs to,
 * and one for accessing the DataWarehouse on internal indices) to allow
 * extrapolation in regions corresponding to vertical angle
 *
 * @tparam T type of Field
 * @tparam STN finite-difference stencil
 * @tparam P list of BC, FC, and Patch::Face packs (FineCoarseInterface are not allowed)
 */
template<typename T, StnType STN, BCF... P >
class bc_vertical_angle_view< ScalarField<T>, STN, P... >
    : virtual public view< ScalarField<T> >
{
private: // TYPES

    /// Type of field
    using Field = ScalarField<T>;

    /// Non const type of the field value
    using V = typename std::remove_const<T>::type;

    /// Boundary variety order (number of boundary faces incident to given region)
    static constexpr size_t ORD = sizeof... ( P );

    static_assert ( ORD >=2, "bc_vertical_angle_view must be used only over edges or vertices" );

private: // MEMBERS

    /// Inner view to grid variable (owned by parent bcs_basic_fd_view)
    const basic_fd_view<Field, STN> * m_dw_view;

    /// Inner view pointers to grid variable (indexed by direction - owned by parent bcs_basic_fd_view)
    const std::vector < basic_fd_view<Field, STN> * > & m_bc_views;

    /// Region where the view is defined
    Region m_support;

public: // CONSTRUCTOR

    /**
     * @brief Constructor
     *
     * Instantiate a view without gathering info from the DataWarehouse
     *
     * @param dw_view view to be used for accessing the DataWarehouse
     * @param bc_views views to be used for accessing ghost regions
     */
    bc_vertical_angle_view (
        const basic_fd_view<Field, STN> * dw_view,
        const std::vector < basic_fd_view<Field, STN> * > & bc_views
    ) : m_dw_view ( dw_view ),
        m_bc_views ( bc_views )
    {
    };

private: // SINGLE FACE METHODS

    /**
     * @brief Set the support to the region opposite to given one with respect
     * to the angle formed by P... (single face implementation)
     *
     * Adjust region bounding indices in the direction orthogonal to Q
     *
     * @tparam Q face boundary pack
     * @return unused value (to allow call in array initialization for variadic
     * template expansions)
     */
    template<BCF Q>
    bool
    set ()
    {
        static_assert ( get_bcf<Q>::bc != BC::FineCoarseInterface, "Unsupported BC" );
        constexpr Patch::FaceType F = get_bcf<Q>::face;
        constexpr DirType D = get_face<F>::dir;
        constexpr int SGN = get_face<F>::sgn;
        if ( SGN > 0 ) // F=D+
        {
            m_support.low() [D] = m_support.high() [D];
            m_support.high() [D] += 1;
        }
        else // F=D-
        {
            m_support.high() [D] = m_support.low() [D];
            m_support.low() [D] -= 1;
        }
        return true;
    }

    /**
     * @brief Compute the index on the grid boundary region corresponding to the
     * given ghost one in the region vertical to the angle formed by P...
     * (single face implementation)
     *
     * Adjust indices in the direction orthogonal to Q
     *
     * @tparam Q face boundary pack
     * @param[in,out] result grid index
     * @return unused value (to allow call in array initialization for variadic
     * template expansions)
     */
    template<BCF Q>
    bool
    opposite (
        IntVector & result
    ) const
    {
        static_assert ( get_bcf<Q>::bc != BC::FineCoarseInterface, "Unsupported BC" );

        constexpr Patch::FaceType F = get_bcf<Q>::face;
        constexpr DirType D = get_face<F>::dir;
        constexpr int SGN = get_face<F>::sgn;
        result[D] -= SGN;
        return true;
    }

    /**
     * @brief Add the value of the view at the index on the ghost boundary region
     * of Q corresponding to the given grid index on the angle formed by P...
     *
     * Adjust indices in the direction orthogonal to Q
     *
     * @tparam Q face boundary pack
     * @param opp grid index
     * @param[in,out] result value to add extrapolated value to
     * @return unused value (to allow call in array initialization for variadic
     * template expansions)
     */
    template<BCF Q>
    bool
    add_adj (
        const IntVector & opp,
        V & result
    ) const
    {
        static_assert ( get_bcf<Q>::bc != BC::FineCoarseInterface, "Unsupported BC" );

        constexpr Patch::FaceType F = get_bcf<Q>::face;
        constexpr DirType D = get_face<F>::dir;
        constexpr int SGN = get_face<F>::sgn;

        const basic_fd_view<Field, STN> & bc_view ( *m_bc_views[D] );

        IntVector adj = opp;
        adj[D] += SGN;
        result += bc_view[adj];
        return true;
    }

#ifdef HAVE_HYPRE
    template<BCF Q>
    bool
    add_adj_entries (
        const IntVector & opp,
        Entries<V> & res,
        const V & w
    ) const
    {
        static_assert ( get_bcf<Q>::bc != BC::FineCoarseInterface, "Unsupported BC" );

        constexpr Patch::FaceType F = get_bcf<Q>::face;
        constexpr DirType D = get_face<F>::dir;
        constexpr int SGN = get_face<F>::sgn;

        const basic_fd_view<Field, STN> & bc_view ( *m_bc_views[D] );

        IntVector adj = opp;
        adj[D] += SGN;
        res.add ( bc_view.entries ( adj ), w );
        return true;
    }
#endif

private: // METHODS

    /**
     * @brief Compute the index on the grid boundary region corresponding to the
     * given ghost one in the region vertical to the angle formed by P...
     *
     * @param id virtual angle index
     * @return grid index
     */
    IntVector
    opposite (
        const IntVector & id
    ) const
    {
        IntVector result ( id );
        std::array<bool, ORD> {{ opposite<P> ( result ) ... }};
        return result;
    }


public: // VIEW METHODS

    /**
     * @brief Retrieve values from the DataWarehouse for a given patch
     * (virtual implementation)
     *
     * @param dw unused
     * @param patch unused
     * @param use_ghosts unused
     */
    virtual void
    set (
        DataWarehouse * dw,
        const Patch * patch,
        bool use_ghosts
    ) override VIRT;

    /**
     * @brief Set the support to the region opposite to given one with respect
     * to the angle formed by P...
     *
     * @remark DataWarehouse is accessed through referenced views which are owned
     * by parent bcs_basic_fd_view so no data is retrieved here
     *
     * @param dw unused
     * @param level unused
     * @param low lower bound of the region to retrieve
     * @param high higher bound of the region to retrieve
     * @param use_ghosts unused
     */
    virtual void
    set (
        DataWarehouse * _DOXYARG ( dw ),
        const Level * _DOXYARG ( level ),
        const IntVector & low,
        const IntVector & high,
        bool _DOXYARG ( use_ghosts )
    ) override
    {
        m_support.low() = low;
        m_support.high() = high;
        std::array<bool, ORD> {{ set<P> () ... }};
    }

    /**
     * @brief Get a copy of the view (virtual implementation)
     *
     * @remark use constructor directly since this view only wraps references
     * to other views
     *
     * @param deep unused
     * @return nothing
     */
    virtual inline view<Field> *
    clone (
        bool _DOXYARG ( deep )
    )
    const override VIRT;

    /**
     * @brief Get a copy of the view and apply translate the support
     * (virtual implementation)
     *
     * @remark use constructor directly since this view only wraps references
     * to other views
     *
     * @param deep unused
     * @param offset unused
     * @return nothing
     */
    virtual inline view<Field> *
    clone (
        bool _DOXYARG ( deep ),
        const IntVector & _DOXYARG ( offset )
    )
    const override VIRT;

    /**
     * @brief Get the region for which the view has access to the DataWarehouse
     *
     * @return support of the view
     */
    virtual inline Support
    get_support()
    const override
    {
        return {{ m_support }};
    };

    /**
     * @brief Check if the view has access to the position with index id
     *
     * @param id position index
     * @return check result
     */
    virtual bool
    is_defined_at (
        const IntVector & id
    ) const override
    {
        IntVector low { m_support.getLow() };
        IntVector high { m_support.getHigh() };

        return ( low[X] <= id[X] && id[X] < high[X] ) &&
               ( low[Y] <= id[Y] && id[Y] < high[Y] ) &&
               ( low[Z] <= id[Z] && id[Z] < high[Z] );
    };

    /**
     * @brief Get/Modify value at position with index id
     * (virtual implementation)
     *
     * Value is extrapolate and cannot be modified
     *
     * @param id unused
     * @return nothing
     */
    virtual T &
    operator[] (
        const IntVector & _DOXYARG ( id )
    ) override VIRT;

    /**
     * @brief Get value at position
     *
     * @param id position index
     * @return field value at id
     */
    virtual V
    operator[] (
        const IntVector & id
    ) const override
    {
        IntVector opp = opposite ( id );
        V result = - ( *m_dw_view ) [opp];
        std::array<bool, ORD> {{ add_adj<P> ( opp, result ) ... }};
        result /= ( ORD - 1 );
        return result;
    };

#ifdef HAVE_HYPRE
    virtual Entries<V>
    entries (
        const IntVector & id
    ) const override
    {
        IntVector opp = opposite ( id );
        V w = 1./(ORD - 1);

        Entries<V> res;
        res.add ( m_dw_view->entries(opp), -w );
        std::array<bool, ORD> {{ add_adj_entries<P> ( opp, res, w ) ... }};
        return res;
    };
#endif

};

} // namespace detail
} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_detail_bcs_basic_fd_view_h
