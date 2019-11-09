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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 * @file CCA/Components/Util/Factory/Factory.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_Util_Factory_Factory_h
#define Packages_Uintah_CCA_Components_Util_Factory_Factory_h

#include <functional>
#include <map>

namespace Uintah
{

/**
 * @brief Generic factory creator class
 *
 * Make possible to create a new instance of a virtual class choosing which
 * derived class implementation to create according to a string
 *
 * @tparam B virtual base class type
 * @tparam Args types of the arguments for the derived classes constructors
 */
template<typename B, typename ... Args>
class Factory
{
public:
    /// Pointer type of the derived classes constructors
    using FactoryMethod = std::function< Base<B> * ( Args... ) >;

    /// Type of the map between strings and derived classes constructors
    using FactoryMap = std::map<std::string, FactoryMethod>;

// protected:
    /// Mapping between strings and derived classes constructors
    static FactoryMap RegisteredNames;

public:
    /**
     * @brief Register to the factory
     *
     * Register a string with a particular derived class constructor
     *
     * @param name string to be registered
     * @param constructor pointer to the derived class constructor
     * @return whether it was added or updated
     */
    static bool
    Register (
        std::string name,
        FactoryMethod constructor
    )
    {
        // add the pair to the map
        auto registeredPair = Factory::RegisteredNames.insert ( std::make_pair ( name.c_str(), constructor ) );
        // return whether it was added or updated
        return registeredPair.second;
    }

    /**
     * @brief Factory create
     *
     * Create a derived class given a string
     *
     * @param name string identifying a registered derived class
     * @param args parameters forwarded to the constructor
     * @return new instance of derived class
     */
    static Base<B> *
    Create (
        std::string name,
        Args ... args
    )
    {
        // attempt to get the pair from the map
        auto registeredPair = Factory::RegisteredNames.find ( name );
        // did we find one?
        if ( registeredPair == Factory::RegisteredNames.end() )
            return nullptr; // return NULL
        // return a new instance of derived class
        return registeredPair->second ( args... );
    }

    /**
     * @brief Factory creator
     *
     * Get the creator for a derived class given a string
     *
     * @param name string identifying a registered derived class
     * @param args parameters forwarded to the constructor
     * @return creator method for derived class
     */
    static FactoryMethod
    Creator (
        std::string name
    )
    {
        // attempt to get the pair from the map
        auto registeredPair = Factory::RegisteredNames.find ( name );
        // did we find one?
        if ( registeredPair == Factory::RegisteredNames.end() )
            return FactoryMethod(); // return empty function
        // return registered creator of derived class
        return registeredPair->second;
    }
}; // class typename

} // namespace Uintah

#endif
