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
 * @file CCA/Components/Solver/HypreSStruct/SolverFactory.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2019/10
 */

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverFactory_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverFactory_h

#include <Core/Util/Factory/Base.h>
#include <Core/Util/Factory/Factory.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Components/Solvers/HypreSStruct/GlobalDataP.h>

namespace Uintah
{

class ProcessorGroup;
class SolverCommon;

namespace HypreSStruct
{

class SStructInterface;

/// Base class for SStructInterface
using SStructInterfaceBase = Base< SStructInterface >;

/// Factory class for SStructInterface
using SStructInterfaceFactory = Factory< SStructInterface, /*const int &,*/ const GlobalDataP & >;

/**
 * @brief Factory class for different HpreSStruct::Solver implementations
 *
 * Factory class for creating new instances of SolverInterface with the 
 * HpreSStruct interface
 */
class SolverFactory
{
public:
    /**
    * @brief factory create method
    *
    * Factory method for creating new instances of SolverInterface
    * with the HpreSStruct interface
    *
    * @param myWorld data structure to manage mpi processes
    * @param solverSpec specifications parsed from ups input file
    */
    static SolverCommon * create (
        const ProcessorGroup * myWorld,
        ProblemSpecP solverSpec
    );
}; // class SolverFactory

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverFactory_h
