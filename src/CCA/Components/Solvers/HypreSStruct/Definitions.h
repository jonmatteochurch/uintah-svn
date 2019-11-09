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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Definitions_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Definitions_h

#define PRINTSYSTEM 1

#include <HYPRE.h>

namespace Uintah
{
namespace HypreSStruct
{

// positive also precond, GMRES default
enum class S : int
{
    GMRES = 0, 
    PCG = -1, 
    FlexGMRES = -2, 
    LGMRES = -3, 
    BiCGSTAB = -4, 
    SysPFMG = 1, 
    Split = 2, 
    FAC = 3, 
    Maxwell = 4
};

// negative non solver, None default
enum class P : int
{
    None = 0, 
    Diag = -1, 
    SysPFMG = ( int ) S::SysPFMG, 
    Split = ( int ) S::Split, 
    FAC = ( int ) S::FAC, 
    Maxwell = ( int ) S::Maxwell
};

enum CoarseSolverType : int
{
    DefaultCoarseSolverType = -1,
    SysPFMG_PCG = 1,
    SysPFMG     = 2
};

enum StructSolverType : int
{
    DefaultStructSolverType = -1,
    SMG = HYPRE_SMG,
    PFMG = HYPRE_PFMG
};

enum RelaxType : int
{
    DefaultRelaxType    = -1,
    Jacobi              =  0,
    WeightedJacobi      =  1,
    RedBlackGaussSeidel =  2  // symmetrix: RB pre-relaxation, BR post-relaxation
};

template <int DIM>
struct SStructStencil
{
    static constexpr int size = 2 * DIM + 1;
    static int offsets[size][DIM];
    static const int entry[size];
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_Definitions_h

