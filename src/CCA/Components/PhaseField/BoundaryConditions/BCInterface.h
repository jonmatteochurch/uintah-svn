/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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
 * @file CCA/Components/PhaseField/BoundaryConditions/BCInterface.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_BCInterface_h
#define Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_BCInterface_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/BoundaryConditions/detail/get_bcs.h>
#include <CCA/Components/PhaseField/BoundaryConditions/detail/partition_range.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWInterface.h>

namespace Uintah
{
namespace PhaseField
{

/**
 * @brief Interface for boundary conditions
 *
 * groups together various methods to get info about boundary conditions which
 * depend on the different types of variable representation and fine-differences
 * stencils allowing to choose the relevant implementation at compile time
 * @tparam VAR type of variable representation
 * @tparam STN finite-difference stencil
 */
template < VarType VAR, StnType STN>
struct BCInterface
{
    /// Problem dimension
    static constexpr DimType DIM = get_stn<STN>::dim;

    /// Problem highest direction
    static constexpr DirType DIR = get_dim<DIM>::highest_dir;

    /// Number of ghost elements required by STN
    static constexpr int GN = get_stn<STN>::ghosts;

    /**
     * @brief Get boundary conditions on multiple faces static
     *
     * Execute the detail::get_bcs < DIM, Field > functor
     *
     * @tparam Field list of type of fields (ScalarField < T > or VectorField < T, N >) of the Problems
     * @param patch grid patch to be checked
     * @param child child index of face boundary condition (as per input file)
     * @param material problem material index
     * @param label variable label to check
     * @param c2f which fine/coarse interface conditions to use on each variable
     * @param[out] flags array of flags to check if any bc is applied to each one of faces
     * @return array of BCInfo
     */
    template <typename Field>
    static std::array< BCInfo<Field>, get_dim<DIM>::face_end >
    get_bcs (
        const Patch * patch,
        const int & child,
        const int & material,
        const typename Field::label_type & label,
        std::array<bool, 2 * get_stn<STN>::dim> & flags,
        const std::map < std::string, FC > * c2f = nullptr
    )
    {
        return detail::get_bcs<DIM, Field>::exec ( patch, child, material, label, c2f, flags );
    }

    /**
     * @brief Partition a patch into a list of Problems
     *
     * Partition a patch accordingly to given boundary conditions
     *
     * @tparam Field list of type of fields (ScalarField < T > or VectorField < T, N >) of the Problems
     * @param labels list of lables for each variable of the problem
     * @param subproblems_label label of subproblems in the DataWarehouse
     * @param material problem material index
     * @param patch grid patch to be partitioned
     * @param bcs list of arrays of BC info for each one of face of the patch for each one of the problem labels
     * @param flags array of flags to check if any bc is applied to each one of faces
     * @return list of partitioned Problems
     */
    template <typename... Field>
    static std::list < Problem<VAR, STN, Field...> >
    partition_patch (
        const typename Field::label_type & ...  labels,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        std::array< BCInfo<Field>, get_dim<DIM>::face_end >... bcs,
        const std::array<bool, 2 * DIM> & flags
    )
    {
        // output container
        std::list < Problem<VAR, STN, Field...> > problems;

        // start partitioning from the highest_dir detail::partition_range iterates on lower directions
        detail::partition_range < DIR, GN, 0, Field... >::template exec<VAR, STN> ( labels..., subproblems_label, material, patch->getLevel(), DWInterface<VAR, DIM>::get_low ( patch ), DWInterface<VAR, DIM>::get_high ( patch ), {}, bcs..., flags, problems );
        return problems;
    }

}; // struct BCInterface

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_BoundaryConditions_BCInterface_h
