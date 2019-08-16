/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#include <CCA/Components/Solvers/HypreFAC/Solver.h>

namespace Uintah
{
namespace HypreFAC
{

DebugStream cout_doing ( "SOLVER ASSEMBLING", false );
DebugStream cout_assembling ( "SOLVER_DOING_COUT", false );

template<> std::string Solver<2>::AdditionalEntriesSuffix = "_extra";

template<> int Solver<2>::offsets[5][2] = {
    {0,0},
    {-1,0},
    {1,0},
    {0,-1},
    {0,1}
};

template<> int Solver<2>::stn0 = 0;
template<> int Solver<2>::stnW = 1;
template<> int Solver<2>::stnE = 2; 
template<> int Solver<2>::stnS = 3;
template<> int Solver<2>::stnN = 4;
template<> int Solver<2>::stnB = -1;
template<> int Solver<2>::stnT = -1;

template<> std::string Solver<3>::AdditionalEntriesSuffix = "_extra";

template<> int Solver<3>::offsets[7][3] = {
    {0,0,0},
    {-1,0,0},
    {1,0,0},
    {0,-1,0},
    {0,1,0},
    {0,0,-1},
    {0,0,1}
};

template<> int Solver<3>::stn0 = 0; 
template<> int Solver<3>::stnW = 1;
template<> int Solver<3>::stnE = 2; 
template<> int Solver<3>::stnS = 3;
template<> int Solver<3>::stnN = 4;
template<> int Solver<3>::stnB = 5;
template<> int Solver<3>::stnT = 6;

} // namespace HypreFAC
} // namespace Uintah



