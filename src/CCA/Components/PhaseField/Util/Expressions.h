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
 * @file CCA/Components/PhaseField/Util/Expressions.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Util_Expressions_h
#define Packages_Uintah_CCA_Components_PhaseField_Util_Expressions_h

namespace Uintah
{
namespace PhaseField
{
#if __cplusplus < 201402L
/**
 * @brief C++11 implementation of C++14 integer sequence
 *
 * The class template integer_sequence represents a compile-time sequence of
 * integers. When used as an argument to a function template, the parameter
 * pack I can be deduced and used in pack expansion.
 *
 * @tparam T an integer type to use for the elements of the sequence
 * @tparam I a non-type parameter pack representing the sequence
 */
template<typename T, T... I>
struct integer_sequence
{
    /// value type
    typedef T value_type;

    /**
     * @brief Returns the number of elements in I
     *
     * Returns the number of elements in Ints. Equivalent to sizeof...(I)
     *
     * @return The number of elements in I
     */
    static constexpr size_t size() noexcept
    {
        return sizeof... ( I );
    }
};

/**
 * A helper alias template index_sequence is defined for the common case where T
 * is std::size_t.
 *
 * @tparam I a non-type parameter pack representing the sequence
 */
template<size_t... I>
using index_sequence = integer_sequence<size_t, I...>;

/**
 * A helper alias template make_integer_sequence is defined to simplify
 * creation of integer_sequence types with 0, 1, 2, ..., N-1 as I.
 *
 * @tparam T an integer type to use for the elements of the sequence
 * @tparam N size of the sequence to create
 * @tparam I a non-type parameter pack representing the sequence
 */
template<typename T, size_t N, T... I>
struct make_integer_sequence
/// @cond DOXYIGNORE
        : make_integer_sequence < T, N - 1, N - 1, I... >
/// @endcond
{};

/**
 * A helper alias template make_integer_sequence is defined to simplify
 * creation of integer_sequence types with 0, 1, 2, ..., N-1 as I.
 * (N=0 implementation)
 *
 * @tparam T an integer type to use for the elements of the sequence
 * @tparam I a non-type parameter pack representing the sequence
 */
template<typename T, T... I>
struct make_integer_sequence<T, 0, I...> : integer_sequence<T, I...> {};

/**
 * A helper alias template make_index_sequence is defined to simplify
 * creation of index_sequence types with 0, 1, 2, ..., N-1 as I.
 * @tparam N size of the sequence to create
 */
template<size_t N>
using make_index_sequence = make_integer_sequence<size_t, N>;

/**
 * A helper alias template index_sequence_for is defined to convert any type
 * parameter pack into an index sequence of the same length
 *
 * @tparam T a type parameter pack for witch to create the sequence
 */
template<typename... T>
using index_sequence_for = make_index_sequence<sizeof... ( T ) >;
#else
using std::integer_sequence;
using std::index_sequence;
using std::make_integer_sequence;
using std::make_index_sequence;
using std::index_sequence_for;
#endif

/**
 * @brief Helper struct template filtering out FineCoarseInterface from a BCF
 * integer_sequence
 *
 * This class template represents a compile-time integer_sequence of BCF s.
 * When given an integer_sequence of any BCF s it result in an integer_sequence
 * where all BCF corresponding to FineCoarseInterface  are filtered out.
 *
 * @tparam Q an integer_sequence of BCF
 * @tparam N the number of BCF from the end of the sequence to not filter
 */
template<typename Q, size_t N, typename = void> struct filtering_bcf;

/**
 * @brief Helper struct template filtering out FineCoarseInterface from a BCF
 * integer_sequence ( base case implementation )
 *
 * When N is equal to the number of BCF there is nothing left to filter then
 * inherit final integer_sequence
 *
 * @tparam P a variadic list of BCF
 * @tparam N the number of BCF from the end of the sequence to not filter
 */
template<size_t N, BCF... P>
struct filtering_bcf< integer_sequence<BCF, P...>, N, typename std::enable_if< ( sizeof... ( P ) == N ) >::type >
: public integer_sequence<BCF, P...>
{};

/**
 * @brief Helper struct template filtering out FineCoarseInterface from a BCF
 * integer_sequence ( not filtering implementation )
 *
 * When P0 is not FineCoarseInterface put it at the end of the list and
 * increment N then recursively inherit
 *
 * @tparam P0 first BCF in the list
 * @tparam P a variadic list of remaining BCF
 * @tparam N the number of BCF from the end of the sequence to not filter
 */
template<size_t N, BCF P0, BCF... P>
struct filtering_bcf < integer_sequence<BCF, P0, P...>, N, typename std::enable_if < ( N <= sizeof... ( P ) && get_bcf<P0>::bc != BC::FineCoarseInterface ) >::type >
: public filtering_bcf < integer_sequence<BCF, P..., P0>, N + 1 >
{};

/**
 * @brief Helper struct template filtering out FineCoarseInterface from a BCF
 * integer_sequence ( filtering implementation )
 *
 * When P0 is FineCoarseInterface remove it from the list then recursively inherit
 *
 * @tparam P0 first BCF in the list
 * @tparam P a variadic list of remaining BCF
 * @tparam N the number of BCF from the end of the sequence to not filter
 */
template<size_t N, BCF P0, BCF... P>
struct filtering_bcf < integer_sequence<BCF, P0, P...>, N, typename std::enable_if < ( N <= sizeof... ( P ) && get_bcf<P0>::bc == BC::FineCoarseInterface ) >::type >
: public filtering_bcf<integer_sequence<BCF, P...>, N>
{};

/**
 * A helper alias template filtering out FineCoarseInterface from a variadic
 * list of BCF s
 *
 * @tparam P a variadic list of BCF
 */
template<BCF... P> using filter_bcf = filtering_bcf < integer_sequence<BCF, P...>, 0 >;

/**
 * @brief Compile time power
 *
 * expression template for computing the power \f$ N^M \f$ at compile time
 * @tparam N integer base \f$ N \f$
 * @tparam M integer exponent \f$ M \f$
 */
template<size_t N, size_t M>
struct power
{
    /// value of \f$ N^M \f$
    static constexpr size_t value = power<N,M/2>::value * power<N,M-M/2>::value;
};

/**
 * @brief Compile time factorial (M = O case)
 *
 * expression template for computing the power \f$ N^0 \f$ at compile time
 * @tparam N integer base \f$ N \f$
 */
template<size_t N>
struct power<N,0>
{
    /// value of \f$ N^0 \f$
    static constexpr size_t value = 1;
};

/**
 * @brief Compile time factorial (M = 1 case)
 *
 * expression template for computing the power \f$ N^1 \f$ at compile time
 * @tparam N integer base \f$ N \f$
 */
template<size_t N>
struct power<N,1>
{
    /// value of \f$ N^1 \f$
    static constexpr size_t value = N;
};

/**
 * @brief Compile time factorial
 *
 * expression template for computing the factorial \f$ N! \f$ at compile time
 * @tparam N integer argument \f$ N \f$
 */
template<size_t N>
struct factorial
{
    /// value of \f$ N! \f$
    static constexpr size_t value = N * factorial < N - 1 >::value;
};

/**
 * @brief Compile time factorial (O case)
 *
 * expression template for computing the factorial \f$ N! \f$ at compile time
 * for \f$ N=0 \f$
 */
template<>
struct factorial<0>
{
    /// value of \f$ 0! \f$
    static constexpr size_t value = 1;
};

/**
 * @brief Compile time combinations
 *
 * Expression template for computing the binomial coefficient \f$ N\choose R \f$
 * at compile time
 * @tparam N set size
 * @tparam R combinations size
 */
template<size_t N, size_t R>
struct combinations
{
    /// value of \f$ N\choose R \f$
    static constexpr size_t value = factorial<N>::value / ( factorial<R>::value * factorial < N - R >::value );
};

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_Util_Expressions_h
