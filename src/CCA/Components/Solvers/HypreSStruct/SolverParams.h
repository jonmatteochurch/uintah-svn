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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverParams_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverParams_h

#include <CCA/Ports/SolverInterface.h>

#include <HYPRE_utilities.h>

#define IGNOREPARAM(solver,param,name) \
DOUT( \
  MPI::Impl::prank( MPI_COMM_WORLD ) == 0 && param->name != -1, \
  "" #solver ": Ignoring '" #name "' parameter." \
)

namespace Uintah
{
namespace HypreSStruct
{

class SolverParams
    : public SolverParameters
{
public:
    int solveFrequency; // Frequency for solving the linear system. timestep % solveFrequency

    int max_levels; // maximum number of FAC levels
    /**
     * (Optional) Set the relative convergence tolerance.
     **/
    HYPRE_Real tol; // convergence tolerance
    /**
     * (Optional) Set the absolute convergence tolerance (default is
     * 0). If one desires the convergence test to check the absolute
     * convergence tolerance {\it only}, then set the relative convergence
     * tolerance to 0.0. (The default convergence test is $ <C*r,r> \leq$
     * max(relative$\_$tolerance$^{2} \ast <C*b, b>$, absolute$\_$tolerance$^2$).)
     **/
    HYPRE_Real abs_tol; // convergence tolerance
    /**
     * (Optional) Set a residual-based convergence tolerance which checks if
     * $\|r_{old}-r_{new}\| < rtol \|b\|$. This is useful when trying to converge to
     * very low relative and/or absolute tolerances, in order to bail-out before
     * roundoff errors affect the approximation.
     **/
    HYPRE_Real res_tol; // convergence tolerance
    HYPRE_Real abstolf; // convergence tolerance
    HYPRE_Real cf_tol; // convergence tolerance

    int min_iter; // minimum number of iterations
    int max_iter; // maximum number of iterations
    /**
     * (Optional) Use the two-norm in stopping criteria.
     **/
    int two_norm;
    int rel_change; // require that the relative difference in successive iterates be small
    int k_dim; // require that the relative difference in successive iterates be small
    int aug_dim; // require that the relative difference in successive iterates be small
    int skip_real_r_check; // require that the relative difference in successive iterates be small
    /**
     * (Optional) Recompute the residual at the end to double-check convergence.
     **/
    int recompute_residual;
    /**
     * (Optional) Periodically recompute the residual while iterating.
     **/
    int recompute_residual_p;
    RelaxType relax_type; // relaxation type
    HYPRE_Real weight; // set Jacobi weight if WeightedJacobi is used
    int num_pre_relax; // number of relaxation sweeps before coarse-grid correction
    int num_post_relax; // number of relaxation sweeps after coarse-grid correction
    int skip_relax; // Skip relaxation on certain grids for isotropic problems. This can greatly improve efficiency by eliminating unnecessary relaxations when the underlying problem is isotropic
    CoarseSolverType csolver_type; // coarsest solver type
    StructSolverType ssolver;
    int logging; // amount of logging to do

    SolverParams() :
        solveFrequency ( 1 ),
        max_levels ( -1 ),
        tol ( -1. ),
        abs_tol ( -1. ),
        abstolf ( -1. ),
        cf_tol ( -1. ),
        min_iter ( -1 ),
        max_iter ( -1 ),
        two_norm ( -1 ),
        rel_change ( -1 ),
        k_dim ( -1 ),
        aug_dim ( -1 ),
        skip_real_r_check ( -1 ),
        recompute_residual ( -1 ),
        recompute_residual_p ( -1 ),
        relax_type ( DefaultRelaxType ),
        weight ( -1. ),
        num_pre_relax ( -1 ),
        num_post_relax ( -1 ),
        skip_relax ( -1 ),
        csolver_type ( DefaultCoarseSolverType ),
        ssolver ( DefaultStructSolverType ),
        logging ( -1 )
    {
        this->setSetupFrequency ( 0 );
        this->setUpdateCoefFrequency ( 0 );
    }

    ~SolverParams()
    {}
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SolverParams_h


