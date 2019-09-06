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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSolverFACParams_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSolverFACParams_h

#include <CCA/Ports/SolverInterface.h>

namespace Uintah
{
namespace HypreFAC
{

enum CoarseSolverType : int
{
    DefaultCoarseSolverType = -1,
    SysPFMG_PCG = 1,
    SysPFMG     = 2
};

enum RelaxType : int
{
    DefaultRelaxType    = -1,
    Jacobi              =  0,
    WeightedJacobi      =  1,
    RedBlackGaussSeidel =  2  // symmetrix: RB pre-relaxation, BR post-relaxation
};

class SolverParams : public SolverParameters
{
public:
    int                 solveFrequency; // Frequency for solving the linear system. timestep % solveFrequency

    int                 max_levels;     // maximum number of FAC levels
    HYPRE_Real          tol;            // convergence tolerance
    int                 max_iter;       // maximum number of iterations
    int                 rel_change;     // require that the relative difference in successive iterates be small
    RelaxType           relax_type;     // relaxation type
    HYPRE_Real          weight;         // set Jacobi weight if WeightedJacobi is used
    int                 num_pre_relax;  // number of relaxation sweeps before coarse-grid correction
    int                 num_post_relax; // number of relaxation sweeps after coarse-grid correction
    CoarseSolverType    csolver_type;   // coarsest solver type
    int                 logging;        // amount of logging to do

    int                 C2F; // restriction operation type: 0 use closest fine (should be NC default), 1 interpolate fine (should be CC default)

    SolverParams()
        : solveFrequency ( -1 )
        , max_levels ( -1 )
        , tol ( -1. )
        , max_iter ( -1 )
        , relax_type ( DefaultRelaxType )
        , weight ( -1 )
        , num_pre_relax ( -1 )
        , num_post_relax ( -1 )
        , csolver_type ( DefaultCoarseSolverType )
        , logging ( -1 )
        , C2F ( -1 )

    {}

    ~SolverParams()
    {}
};

} // namespace HypreFAC
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSolverFAC_h


