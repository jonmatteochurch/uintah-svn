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
 * @file CCA/Components/PhaseField/Solver/HypreSStruct/SolverFactory.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#include <CCA/Components/Solvers/HypreSStruct/SolverFactory.h>
#include <CCA/Components/Solvers/HypreSStruct/GlobalData.h>
#include <CCA/Components/Solvers/HypreSStruct/Solver.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;
using namespace HypreSStruct;

SolverCommon *
SolverFactory::create (
    const ProcessorGroup * myWorld,
    ProblemSpecP solvSpec
)
{
    std::string ndim, c2f;
    solvSpec->getAttribute ( "ndim", ndim );
    solvSpec->getAttribute ( "c2f", c2f );

    std::string solver, precond;
    ProblemSpecP solvParam = solvSpec->getFirstChild();
    solvParam->getWithDefault ( "solver", solver, "fac" );
    solvParam->getWithDefault ( "preconditioner", precond, "none" );
    std::string impl = solver + "|"  + ( precond == "none" ? "" : precond + "|" ) + ndim + "|" + c2f;

    SStructInterfaceFactory::FactoryMethod interface_creator;
    if ( solver == "fac" )
    {
        int csolver_type;
        solvParam->getWithDefault ( "csolver_type", csolver_type, 1 );
        std::string cimpl;
        switch ( csolver_type )
        {
        case 1: // SysPFMG-PCG
            cimpl = "pcg|syspfmg|" + ndim + "|" + c2f;
            break;
        case 2: // SysPFMG
            cimpl = "syspfmg|" + ndim + "|" + c2f;
            break;
        default:
            SCI_THROW ( ProblemSetupException ( "Cannot Create HypreSStruct FAC coarse solver with csolver_type='" + std::to_string ( csolver_type ) + "'", __FILE__, __LINE__ ) );
        }

        GlobalDataP gdata;
        auto creator = SStructInterfaceFactory::Creator ( impl );
        if ( !creator ) SCI_THROW ( ProblemSetupException ( "Cannot Create HypreSStruct Interface impl='" + impl + "'", __FILE__, __LINE__ ) );
        auto ccreator = SStructInterfaceFactory::Creator ( cimpl );
        if ( !ccreator ) SCI_THROW ( ProblemSetupException ( "Cannot Create HypreSStruct Interface impl='" + cimpl + "'", __FILE__, __LINE__ ) );
        interface_creator = [creator, ccreator] ( const GlobalDataP & gdata )
        {
            if ( gdata->nParts() > 1 )
                return creator ( gdata );
            else
                return ccreator ( gdata );
        };
    }
    else
    {
        interface_creator = SStructInterfaceFactory::Creator ( impl );
    }
    if ( !interface_creator ) SCI_THROW ( ProblemSetupException ( "Cannot Create HypreSStruct Interface impl='" + impl + "'", __FILE__, __LINE__ ) );

    if ( ndim == "2" )
        return scinew Solver<2> ( myWorld, std::move ( interface_creator ) );

    if ( ndim == "3" )
        return scinew Solver<3> ( myWorld, std::move ( interface_creator ) );

    SCI_THROW ( ProblemSetupException ( "Cannot Create HypreSStruct Solver with ndim='" + ndim + "'", __FILE__, __LINE__ ) );
    return nullptr;
}

