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
 * @file CCA/Components/PhaseField/AMR/detail/amr_interpolator.h
 * @author Jon Matteo Church
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_AMR_detail_amr_interpolator_h
#define Packages_Uintah_CCA_Components_PhaseField_AMR_detail_amr_interpolator_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/Util/Expressions.h>

namespace Uintah
{
namespace PhaseField
{
namespace detail
{

/**
 * @brief Abstract wrapper of grid variables for interpolation from coarser to
 *        finer levels.
 *
 * Adds to view the possibility to compute multi-grid interpolation
 *
 * @remark All different interpolation strategies must specialize this class and
 *         implement the view< T > class
 *
 * @tparam Field type of Field (should be only ScalarField)
 * @tparam Problem type of PhaseField problem
 * @tparam Index index_sequence of Field within Problem (first element is variable index,
 * following ones, if present, are the component index within the variable)
 * @tparam FCI order of interpolation
 * @tparam DIM problem dimension
 */
template<typename Field, typename Problem, typename Index, FCIType FCI, DimType DIM> class amr_interpolator;

/**
 * @brief Abstract wrapper of grid variables for interpolation from coarser to
 * finer levels. (VectorField implementation)
 *
 * Adds to view the possibility to compute multi-grid interpolation
 *
 * @tparam T type of each component of the field at each point
 * @tparam N number of components
 * @tparam Problem type of PhaseField problem
 * @tparam Index index_sequence of Field within Problem (first element is variable index,
 * following ones, if present, are the component index within the variable)
 * @tparam FCI order of interpolation
 * @tparam DIM problem dimension
 */
template<typename T, size_t N, typename Problem, size_t... I, FCIType FCI, DimType DIM>
class amr_interpolator < VectorField<T, N>, Problem, index_sequence<I...>, FCI, DIM >
    : virtual public view < VectorField<T, N> >
    , virtual public view_array < view < ScalarField<T> >, ScalarField<T>, N >
{
private: // TYPES

    /// Type of field
    using Field = VectorField<T, N>;

    /// Type for component sub indices
    template <size_t J>
    using SubIndex = index_sequence<I..., J>;

    /// Type of View of each component
    template <size_t J>
    using View = amr_interpolator < ScalarField<T>, Problem, SubIndex<J>, FCI, DIM >;

private: // SINGLE INDEX METHODS

    /**
     * @brief Create amr_interpolator view for a given component
     *
     * Instantiate a view without gathering info from the DataWarehouse
     *
     * @tparam J component index
     * @param label label of variable in the DataWarehouse
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material index of material in the DataWarehouse
     * @return pointer to the newly created view
     */
    template<size_t J>
    void * create_element (
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material
    )
    {
        return this->m_view_ptr[J] = scinew View<J> ( label[J], subproblems_label, material );
    }

    /**
     * @brief Create amr_interpolator view for a given component
     *
     * Instantiate a view and gather info from dw
     *
     * @tparam J component index
     * @param dw DataWarehouse from which data is retrieved
     * @param label label of variable in the DataWarehouse
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material index of material in the DataWarehouse
     * @param patch grid patch on which data is retrieved
     * @param use_ghosts if ghosts value are to be retrieved
     * @return pointer to the newly created view
     */
    template<size_t J>
    void * create_element (
        DataWarehouse * dw,
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        bool use_ghosts
    )
    {
        return this->m_view_ptr[J] = scinew View<J> ( dw, label[J], subproblems_label, material, patch, use_ghosts );
    }

private: // INDEXED CONSTRUCTOR

    /**
     * @brief Indexed constructor
     *
     * Instantiate a view without gathering info from the DataWarehouse
     *
     * @tparam J indices for boundary views
     * @param unused to allow template argument deduction
     * @param label label of variable in the DataWarehouse
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material index of material in the DataWarehouse
     */
    template<size_t... J>
    amr_interpolator (
        index_sequence<J...> _DOXYARG ( unused ),
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material
    )
    {
        std::array<void *, N> {{ create_element<J> ( label, subproblems_label, material )... }};
    }

    /**
     * @brief Indexed constructor
     *
     * Instantiate a view and gather info from dw
     *
     * @tparam J indices for boundary views
     * @param unused to allow template argument deduction
     * @param dw DataWarehouse from which data is retrieved
     * @param label label of variable in the DataWarehouse
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material index of material in the DataWarehouse
     * @param patch grid patch
     * @param use_ghosts if ghosts value are to be retrieved
     */
    template<size_t... J>
    amr_interpolator (
        index_sequence<J...> _DOXYARG ( unused ),
        DataWarehouse * dw,
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        bool use_ghosts
    )
    {
        std::array<void *, N> {{ create_element<J> ( dw, label, subproblems_label, material, patch, use_ghosts )... }};
    }

private: // COPY CONSTRUCTOR

    /**
     * @brief Constructor
     *
     * Instantiate a copy of a given view
     *
     * @param copy source view for copying
     * @param deep if true inner grid variable is copied as well otherwise the
     * same grid variable is referenced
     */
    amr_interpolator (
        const amr_interpolator * copy,
        bool deep
    )
    {
        for ( size_t i = 0; i < N; ++i )
        {
            const auto & v = ( *copy ) [i];
            this->m_view_ptr[i] = v.clone ( deep );
        }
    }

public: // CONSTRUCTORS/DESTRUCTOR

    /**
     * @brief Constructor
     *
     * Instantiate amr_interpolator components without gathering info from the DataWarehouse
     *
     * @param label list of variable labels for each component
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material material index
     */
    amr_interpolator (
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material
    ) : amr_interpolator ( make_index_sequence<N> {}, label, subproblems_label, material )
    {}

    /**
     * @brief Constructor
     *
     * Instantiate amr_interpolator components and gather info from dw
     *
     * @param dw DataWarehouse from which data is retrieved
     * @param label list of variable labels for each component
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material material index
     * @param patch grid patch
     * @param use_ghosts if ghosts value are to be retrieved
     */
    amr_interpolator (
        DataWarehouse * dw,
        const typename Field::label_type & label,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        bool use_ghosts = View<0>::use_ghosts_dflt
    ) : amr_interpolator ( make_index_sequence<N> {}, dw, label, subproblems_label, material, patch, use_ghosts )
    {}

    /// Destructor
    virtual ~amr_interpolator()
    {
        for ( auto view : this->m_view_ptr )
            delete view;
    }

    /// Prevent copy (and move) constructor
    amr_interpolator ( const amr_interpolator & ) = delete;

    /// Prevent copy (and move) assignment
    /// @return deleted
    amr_interpolator & operator= ( const amr_interpolator & ) = delete;

public: // VIEW METHODS

    /**
     * @brief Get a copy of the view
     *
     * @param deep if true inner grid variable is copied as well otherwise the
     * same grid variable is referenced
     *
     * @return new view instance
     */
    virtual view<Field> *
    clone (
        bool deep
    )
    const override
    {
        return scinew amr_interpolator ( this, deep );
    };
};

} // namespace detail
} // namespace PhaseField
} // namespace Uintah

#include <CCA/Components/PhaseField/AMR/detail/amr_interpolator_I0.h>
#include <CCA/Components/PhaseField/AMR/detail/amr_interpolator_I1_D1.h>
#include <CCA/Components/PhaseField/AMR/detail/amr_interpolator_I1_D2.h>
#include <CCA/Components/PhaseField/AMR/detail/amr_interpolator_I1_D3.h>

#endif // Packages_Uintah_CCA_Components_PhaseField_AMR_detail_amr_interpolator_h
