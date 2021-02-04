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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_solver_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_detail_sstruct_solver_h

#include <CCA/Components/Solvers/HypreSStruct/Definitions.h>

namespace Uintah
{
namespace HypreSStruct
{
namespace detail
{

template <int SLV, int DIM, int C2F, bool precond> class sstruct_solver;

template<S SLV, int DIM, int C2F> using sstruct_ssolver = sstruct_solver<(int)SLV, DIM, C2F, false>;
template<P PCN, int DIM, int C2F> using sstruct_precond = sstruct_solver<(int)PCN, DIM, C2F, true>;

} // namespace detail
} // namespace HypreSStruct
} // namespace Uintah

#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_Diagonal.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_SysPFMG.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_Split.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_FAC.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_PCG.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_GMRES.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_FlexGMRES.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_LGMRES.h>
#include <CCA/Components/Solvers/HypreSStruct/detail/sstruct_BiCGSTAB.h>

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_PartDataP_h


