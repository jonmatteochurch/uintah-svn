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
 * @file CCA/Components/PhaseField/DataTypes/BCInfo.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_BCInfo_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_BCInfo_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <CCA/Components/PhaseField/DataTypes/ScalarField.h>
#include <CCA/Components/PhaseField/DataTypes/VectorField.h>

#include <cstddef>

namespace Uintah
{
namespace PhaseField
{

/**
 * @brief Boundary Condition Information
 *
 * Stores the type of boundary conditions, o fine coarse interpolation and the
 * value to impose on a boundary face
 *
 * @tparam Field type of field of the variable on which the condition is applied
 */
template<typename Field> struct BCInfo;

/**
 * @brief Boundary Condition Information (ScalarField implementation)
 *
 * Stores the type of boundary conditions, o fine coarse interpolation and the
 * value to impose on a boundary face
 *
 * @tparam T type of the field value at each point
 */
template<typename T>
struct BCInfo< ScalarField<T> >
{
private: // STATIC MEMBERS

    /// Register
    /// @remark Cause the (instantiated) type of this class to be registered with the Core/Disclosure/TypeDescription class when this class.
    /// @remark It is not used for anything else in the
    // program.
    static TypeDescription::Register registerMe;

    /// Unique instance of TypeDescription
    static TypeDescription * td;

public: // MEMBERS

    /// value to impose on boundary
    typename std::remove_const<T>::type value;

    /// Type of boundary conditions
    BC bc;

    /// type of fine/coarse interface conditions
    FC c2f;

private: // STATIC METHODS

    /**
     * @brief Create and commit MPI_Datatype for BCInfo
     *
     * Define a new struct datatype for communicating BCInfo with MPI and commit it

     * @remark the reference to this private method is passed to the TypeDescription
     * constructor to enure a single instance of MPI_Datatype is created
     *
     * @return MPI_Datatype for the BCInfo struct instance
     */
    static MPI_Datatype
    make_mpitype()
    {
        static constexpr int count = 4;

        MPI_Datatype datatype;
        int blocklengths[count] = {1, 1, 1, 1};

        MPI_Aint displacements[count] =
        {
            offsetof ( BCInfo, value ),
            offsetof ( BCInfo, bc ),
            offsetof ( BCInfo, c2f ),
            sizeof ( BCInfo )
        };

        MPI_Datatype types[count] =
        {
            PhaseField::getTypeDescription<T>()->getMPIType(),
            PhaseField::getTypeDescription<size_t>()->getMPIType(),
            PhaseField::getTypeDescription<size_t>()->getMPIType(),
            MPI_UB
        };

        Uintah::MPI::Type_create_struct ( count, blocklengths, displacements, types, &datatype );
        Uintah::MPI::Type_commit ( &datatype );

        return datatype;
    }

public: // STATIC METHODS

    /**
     * @brief Get the TypeDescription for BCInfo
     *
     * Get the pointer to the single instance of TypeDescription for the BCInfo
     * type instance. If this is not yet initialized then a new instance is created.
     *
     * @return Pointer to the TypeDescription for the BCInfo struct instance
     */
    static const TypeDescription *
    getTypeDescription()
    {
        if ( !td )
        {
            td = scinew TypeDescription ( TypeDescription::Other, "BCInfo", true, BCInfo::make_mpitype );
        }
        return td;
    }
};

/**
 * @brief Boundary Condition Information (VectorField implementation)
 *
 * Stores the type of boundary conditions, o fine coarse interpolation and the
 * value to impose on a boundary face
 *
 * @tparam T type of each component of the field at each point
 * @tparam N number of components
 */
template<typename T, size_t N>
struct BCInfo< VectorField<T, N> >
{
private: // STATIC MEMBERS

    /// Register
    /// @remark Cause the (instantiated) type of this class to be registered with the Core/Disclosure/TypeDescription class when this class.
    /// @remark It is not used for anything else in the
    // program.
    static TypeDescription::Register registerMe;

    /// Unique instance of TypeDescription
    static TypeDescription * td;

public: // MEMBERS

    /// array of values to impose on boundary
    std::array < typename std::remove_const<T>::type, N > value;

    /// Type of boundary conditions
    BC bc;

    /// type of fine/coarse interface conditions
    FC c2f;

public: // CONSTRUCTORS

    /**
     * @brief Constructor
     *
     * From a list of BCInfo
     *
     * @tparam I0 first BCInfo type
     * @tparam I following BCInfo types
     * @param i0 first BCInfo
     * @param i following BCInfo's
     */
    template<typename I0, typename... I>
    BCInfo ( I0 && i0, I && ... i ) :
        value { std::forward<T> ( i0.value ), std::forward<T> ( i.value )... },
          bc ( i0.bc ),
          c2f ( i0.c2f )
    {};

    /**
     * @brief Default Constructor
     */
    BCInfo() = default;

private: // STATIC METHODS

    /**
     * @brief Create and commit MPI_Datatype for BCInfo
     *
     * Define a new struct datatype for communicating BCInfo with MPI and commit it

     * @remark the reference to this private method is passed to the TypeDescription 
     * constructor to enure a single instance of MPI_Datatype is created
     *
     * @return MPI_Datatype for the BCInfo struct instance
     */
    static MPI_Datatype
    make_mpitype()
    {
        static constexpr int count = 4;

        MPI_Datatype datatype;
        int blocklengths[count] = {N, 1, 1, 1};

        MPI_Aint displacements[count] =
        {
            offsetof ( BCInfo, value ),
            offsetof ( BCInfo, bc ),
            offsetof ( BCInfo, c2f ),
            sizeof ( BCInfo )
        };

        MPI_Datatype types[count] =
        {
            getMPIType<T>(),
            getMPIType<size_t>(),
            getMPIType<size_t>(),
            MPI_UB
        };

        Uintah::MPI::Type_create_struct ( count, blocklengths, displacements, types, &datatype );
        Uintah::MPI::Type_commit ( &datatype );

        return datatype;
    }

public: // STATIC METHODS

    /**
     * @brief Get the TypeDescription for BCInfo
     *
     * Get the pointer to the single instance of TypeDescription for the BCInfo
     * type instance. If this is not yet initialized then a new instance is created.
     *
     * @return Pointer to the TypeDescription for the BCInfo struct instance
     */
    static const TypeDescription *
    getTypeDescription()
    {
        if ( !td )
        {
            td = scinew TypeDescription ( TypeDescription::Other, "SubProblemsMpiData", true, BCInfo::make_mpitype );
        }
        return td;
    }
};

template<typename T> TypeDescription * BCInfo < ScalarField<T> >::td = nullptr;

template<typename T> TypeDescription::Register BCInfo< ScalarField<T> >::registerMe ( getTypeDescription() );

template<typename T, size_t N> TypeDescription * BCInfo < VectorField<T, N> >::td = nullptr;

template<typename T, size_t N> TypeDescription::Register BCInfo< VectorField<T, N> >::registerMe ( getTypeDescription() );

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_BCInfo_h
