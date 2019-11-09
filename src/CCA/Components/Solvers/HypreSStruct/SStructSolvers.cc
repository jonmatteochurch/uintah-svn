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
 * @file CCA/Components/Solvers/HypreSStruct/SStructSolvers.cc
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2019/10
 *
 * In this the names used by ApplicationFactory for the differnt implementations
 * of the Heat application are defined as well as their explicit instantiation.
 */

#include <CCA/Components/Solvers/HypreSStruct/SStructSolver.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructImplementation.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructPCG.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructSysPFMG.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructFAC.h>

namespace Uintah
{
namespace HypreSStruct
{

/// @cond DOXYIGNORE
template<> const std::string SStructSolver<S::PCG, 2, 0>::Name = "pcg|2|0";
template<> const std::string SStructSolver<S::PCG, 2, 1>::Name = "pcg|2|1";
template<> const std::string SStructSolver<S::PCG, 3, 0>::Name = "pcg|3|0";
template<> const std::string SStructSolver<S::PCG, 3, 1>::Name = "pcg|3|1";

template<> const std::string SStructSolver<S::SysPFMG, 2, 0>::Name = "sys_pfmg|2|0";
template<> const std::string SStructSolver<S::SysPFMG, 2, 1>::Name = "sys_pfmg|2|1";
template<> const std::string SStructSolver<S::SysPFMG, 3, 0>::Name = "sys_pfmg|3|0";
template<> const std::string SStructSolver<S::SysPFMG, 3, 1>::Name = "sys_pfmg|3|1";

template<> const std::string SStructSolver<S::FAC, 2, 0>::Name = "fac|2|0";
template<> const std::string SStructSolver<S::FAC, 2, 1>::Name = "fac|2|1";
template<> const std::string SStructSolver<S::FAC, 3, 0>::Name = "fac|3|0";
template<> const std::string SStructSolver<S::FAC, 3, 1>::Name = "fac|3|1";

template<> const std::string SStructSolver<S::PCG, 2, 0, P::SysPFMG>::Name = "pcg|sys_pfmg|2|0";
template<> const std::string SStructSolver<S::PCG, 2, 1, P::SysPFMG>::Name = "pcg|sys_pfmg|2|1";
template<> const std::string SStructSolver<S::PCG, 3, 0, P::SysPFMG>::Name = "pcg|sys_pfmg|3|0";
template<> const std::string SStructSolver<S::PCG, 3, 1, P::SysPFMG>::Name = "pcg|sys_pfmg|3|1";

template class SStructSolver<S::PCG, 2, 0>;
template class SStructSolver<S::PCG, 2, 1>;
template class SStructSolver<S::PCG, 3, 0>;
template class SStructSolver<S::PCG, 3, 1>;

template class SStructSolver<S::FAC, 2, 0>;
template class SStructSolver<S::FAC, 2, 1>;
template class SStructSolver<S::FAC, 3, 0>;
template class SStructSolver<S::FAC, 3, 1>;

template class SStructSolver<S::SysPFMG, 2, 0>;
template class SStructSolver<S::SysPFMG, 2, 1>;
template class SStructSolver<S::SysPFMG, 3, 0>;
template class SStructSolver<S::SysPFMG, 3, 1>;

template class SStructSolver<S::PCG, 2, 0, P::SysPFMG>;
template class SStructSolver<S::PCG, 2, 1, P::SysPFMG>;
template class SStructSolver<S::PCG, 3, 0, P::SysPFMG>;
template class SStructSolver<S::PCG, 3, 1, P::SysPFMG>;
/// @endcond

} // namespace HypreSStruct
} // namespace PhaseField
