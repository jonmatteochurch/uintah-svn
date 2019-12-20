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
#include <iostream>

// #define FACTORY_DBG

namespace Uintah
{

#ifdef FACTORY_DBG
struct FactoryString
        : public std::string
{
    FactoryString ( const char str[] )
        : std::string ( str )
    {
        std::cout << "Created string `" << std::string ( str ) << "`" << std::endl;
    }
};
#else
using FactoryString = std::string;
#endif

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
    using FactoryMethod = std::function< B * ( Args... ) >;

    /// Type of the map between strings and derived classes constructors
#ifdef FACTORY_DBG
    struct FactoryMap
            : public std::map<std::string, FactoryMethod>
    {
        FactoryMap() : std::map<std::string, FactoryMethod>()
        {
            std::cout << "Created factory map for " << typeid ( B ).name() << std::endl;
        }
    };
#else
    using FactoryMap = std::map<std::string, FactoryMethod>;
#endif

protected:
    /// Mapping between strings and derived classes constructors
    static FactoryMap m_registry;

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
#ifdef FACTORY_DBG
        std::cout << "Registering `" << name << "` to `" << typeid ( B ).name() << "` factory" << std::endl;
#endif
        // add the pair to the map
        auto it = m_registry.emplace ( name, constructor );
        // return whether it was added or updated
        return it.second;
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
    static B *
    Create (
        std::string name,
        Args ... args
    )
    {
#ifdef FACTORY_DBG
        std::cout << "Factory creating `" << name << "` instance of `" << typeid ( B ).name()  << "`... ";
#endif
        // attempt to get the pair from the map
        auto it = m_registry.find ( name );
        // did we find one?
        if ( it == m_registry.end() )
        {
#ifdef FACTORY_DBG
            std::cout << "FAILED" << std::endl;
#endif
            return nullptr; // return NULL
        }
        // return a new instance of derived class
#ifdef FACTORY_DBG
        std::cout << "SUCCESS" << std::endl;
#endif
        return it->second ( args... );
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
        auto it = m_registry.find ( name );
        // did we find one?
        if ( it == m_registry.end() )
            return FactoryMethod(); // return empty function
        // return registered creator of derived class
        return it->second;
    }
}; // class typename

// WARN: still needs to be referenced to be initialized
template<typename B, typename ... Args>
typename Factory<B, Args...>::FactoryMap Factory<B, Args...>::m_registry = {};

} // namespace Uintah

#endif
